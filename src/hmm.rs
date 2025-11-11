use serde::Serialize;
use std::collections::BTreeMap;
use std::collections::VecDeque;

#[derive(Debug, Serialize)]
pub struct HeaderField {
    pub key: String,
    pub value: String,
}

#[derive(Debug, Serialize)]
pub struct StatsLine {
    pub scope: String,
    pub method: String,
    pub values: Vec<Option<f64>>,
    pub raw_values: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct HmmComposition {
    pub match_emissions: Vec<Option<f64>>,
    pub insert_emissions: Vec<Option<f64>>,
    pub transitions: Vec<Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct HmmNode {
    pub index: usize,
    pub match_emissions: Vec<Option<f64>>,
    pub insert_emissions: Vec<Option<f64>>,
    pub transitions: Vec<Option<f64>>,
    pub annotations: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct HmmProfile {
    pub header: Vec<HeaderField>,
    pub header_map: BTreeMap<String, String>,
    pub stats: Vec<StatsLine>,
    pub alphabet: Vec<String>,
    pub transition_labels: Vec<String>,
    pub composition: Option<HmmComposition>,
    pub nodes: Vec<HmmNode>,
}

pub fn parse_hmm(input: &str) -> Result<HmmProfile, String> {
    let mut header_fields = Vec::new();
    let mut header_map = BTreeMap::new();
    let mut stats = Vec::new();
    let mut alphabet = Vec::new();
    let mut transition_labels = Vec::new();

    let mut lines = input.lines();

    while let Some(line) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "//" {
            return Err("Unexpected end of HMM before body".into());
        }
        if trimmed.starts_with("HMM") {
            let mut tokens = trimmed.split_whitespace();
            tokens.next(); // skip HMM keyword
            alphabet = tokens.map(|t| t.to_string()).collect();
            let transition_line = lines
                .next()
                .ok_or_else(|| "Missing transition label line".to_string())?;
            transition_labels = transition_line
                .trim()
                .split_whitespace()
                .map(|t| t.to_string())
                .collect();
            break;
        }
        if trimmed.starts_with("STATS") {
            stats.push(parse_stats_line(trimmed)?);
            continue;
        }

        let key = trimmed
            .split_whitespace()
            .next()
            .ok_or_else(|| "Invalid header line".to_string())?;
        let value = trimmed[key.len()..].trim().to_string();
        header_fields.push(HeaderField {
            key: key.to_string(),
            value: value.clone(),
        });
        header_map.insert(key.to_string(), value);
    }

    if alphabet.is_empty() {
        return Err("Failed to locate HMM alphabet".into());
    }

    let mut remaining: VecDeque<&str> = lines.collect();

    let mut composition = None;
    let mut nodes = Vec::new();

    while let Some(line) = remaining.pop_front() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "//" {
            break;
        }
        if trimmed.starts_with("STATS") {
            stats.push(parse_stats_line(trimmed)?);
            continue;
        }
        if trimmed.starts_with("COMPO") {
            if composition.is_some() {
                return Err("Multiple COMPO sections found".into());
            }
            composition = Some(parse_composition(
                trimmed,
                &mut remaining,
                alphabet.len(),
                transition_labels.len(),
            )?);
            continue;
        }

        if trimmed.as_bytes()[0].is_ascii_digit() {
            let node = parse_node(
                trimmed,
                &mut remaining,
                alphabet.len(),
                transition_labels.len(),
            )?;
            nodes.push(node);
            continue;
        }

        return Err(format!("Unrecognized line in HMM body: {trimmed}"));
    }

    Ok(HmmProfile {
        header: header_fields,
        header_map,
        stats,
        alphabet,
        transition_labels,
        composition,
        nodes,
    })
}

fn parse_stats_line(line: &str) -> Result<StatsLine, String> {
    let mut tokens = line.split_whitespace();
    tokens.next(); // STATS
    let scope = tokens
        .next()
        .ok_or_else(|| "Missing STATS scope".to_string())?
        .to_string();
    let method = tokens
        .next()
        .ok_or_else(|| "Missing STATS method".to_string())?
        .to_string();

    let mut values = Vec::new();
    let mut raw_values = Vec::new();
    for token in tokens {
        raw_values.push(token.to_string());
        values.push(parse_optional_number(token)?);
    }

    Ok(StatsLine {
        scope,
        method,
        values,
        raw_values,
    })
}

fn parse_composition(
    line: &str,
    remaining: &mut VecDeque<&str>,
    alphabet_len: usize,
    transition_len: usize,
) -> Result<HmmComposition, String> {
    let mut tokens = line.split_whitespace();
    let label = tokens
        .next()
        .ok_or_else(|| "Malformed COMPO line".to_string())?;
    if label != "COMPO" {
        return Err("Expected COMPO label".into());
    }

    let match_emissions = collect_scores(tokens, alphabet_len)?;

    let insert_line =
        pop_nonempty(remaining).ok_or_else(|| "Missing COMPO insert line".to_string())?;
    let insert_tokens = insert_line.split_whitespace();
    let insert_emissions = collect_scores(insert_tokens, alphabet_len)?;

    let transition_line =
        pop_nonempty(remaining).ok_or_else(|| "Missing COMPO transition line".to_string())?;
    let transition_tokens = transition_line.split_whitespace();
    let transitions = collect_scores(transition_tokens, transition_len)?;

    Ok(HmmComposition {
        match_emissions,
        insert_emissions,
        transitions,
    })
}

fn parse_node(
    line: &str,
    remaining: &mut VecDeque<&str>,
    alphabet_len: usize,
    transition_len: usize,
) -> Result<HmmNode, String> {
    let mut tokens = line.split_whitespace();
    let index_token = tokens
        .next()
        .ok_or_else(|| "Missing node index".to_string())?;
    let index = index_token
        .parse::<usize>()
        .map_err(|_| format!("Invalid node index: {index_token}"))?;

    let match_emissions = collect_scores(&mut tokens, alphabet_len)?;
    let annotations = tokens.map(|t| t.to_string()).collect();

    let insert_line =
        pop_nonempty(remaining).ok_or_else(|| format!("Missing insert line for node {index}"))?;
    if insert_line.trim() == "//" {
        return Err(format!("Unexpected end of HMM while parsing node {index}"));
    }
    let insert_tokens = insert_line.split_whitespace();
    let insert_emissions = collect_scores(insert_tokens, alphabet_len)?;

    let transition_line = pop_nonempty(remaining)
        .ok_or_else(|| format!("Missing transition line for node {index}"))?;
    if transition_line.trim() == "//" {
        return Err(format!("Unexpected end of HMM while parsing node {index}"));
    }
    let transition_tokens = transition_line.split_whitespace();
    let transitions = collect_scores(transition_tokens, transition_len)?;

    Ok(HmmNode {
        index,
        match_emissions,
        insert_emissions,
        transitions,
        annotations,
    })
}

fn collect_scores<'a, I>(tokens: I, expected: usize) -> Result<Vec<Option<f64>>, String>
where
    I: IntoIterator<Item = &'a str>,
{
    let mut result = Vec::with_capacity(expected);
    for token in tokens.into_iter().take(expected) {
        result.push(parse_optional_number(token)?);
    }
    if result.len() != expected {
        return Err(format!(
            "Expected {expected} score values, found {}",
            result.len()
        ));
    }
    Ok(result)
}

fn parse_optional_number(token: &str) -> Result<Option<f64>, String> {
    if token == "*" {
        return Ok(None);
    }
    token
        .parse::<f64>()
        .map(Some)
        .map_err(|_| format!("Invalid numeric value: {token}"))
}

fn pop_nonempty<'a>(lines: &mut VecDeque<&'a str>) -> Option<&'a str> {
    while let Some(line) = lines.pop_front() {
        if !line.trim().is_empty() {
            return Some(line);
        }
    }
    None
}
