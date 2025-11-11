use serde::Serialize;
use std::iter::Peekable;

#[derive(Debug, Serialize)]
pub struct HeaderField {
    pub key: String,
    pub value: String,
}

#[derive(Debug, Serialize)]
pub struct HmmState {
    pub label: String,
    pub match_emissions: Vec<Option<f64>>,
    pub annotation: Vec<String>,
    pub insert_emissions: Vec<Option<f64>>,
    pub transitions: Vec<Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct Hmm {
    pub format_line: String,
    pub metadata: Vec<HeaderField>,
    pub alphabet: Vec<String>,
    pub transition_order: Vec<String>,
    pub states: Vec<HmmState>,
}

pub fn parse_hmm(input: &str) -> Result<Hmm, String> {
    let mut lines = input
        .lines()
        .map(|line| line.trim_end_matches('\r'))
        .peekable();

    let format_line = next_nonempty_line(&mut lines)
        .ok_or_else(|| "HMM text is empty".to_string())?
        .trim()
        .to_string();

    let mut metadata = Vec::new();
    let header_line = loop {
        let line = lines
            .next()
            .ok_or_else(|| "Unexpected end of input before HMM section".to_string())?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("HMM") {
            break trimmed.to_string();
        }
        let mut parts = trimmed.split_whitespace();
        let key = parts
            .next()
            .ok_or_else(|| "Malformed metadata line".to_string())?
            .to_string();
        let value = parts.collect::<Vec<_>>().join(" ");
        metadata.push(HeaderField { key, value });
    };

    let header_tokens: Vec<_> = header_line.split_whitespace().collect();
    if header_tokens.len() < 2 {
        return Err("HMM header line is missing the alphabet".into());
    }
    if header_tokens[0] != "HMM" {
        return Err("HMM header line must start with 'HMM'".into());
    }
    let alphabet = header_tokens[1..]
        .iter()
        .map(|sym| sym.to_string())
        .collect::<Vec<_>>();
    if alphabet.is_empty() {
        return Err("Alphabet cannot be empty".into());
    }

    let transition_header = loop {
        let line = lines
            .next()
            .ok_or_else(|| "Missing transition column header".to_string())?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "//" {
            return Err("Encountered end of record before transition header".into());
        }
        break trimmed.to_string();
    };
    let transition_order = transition_header
        .split_whitespace()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    if transition_order.is_empty() {
        return Err("Transition header must list at least one column".into());
    }

    let mut states = Vec::new();
    let mut saw_terminator = false;
    while let Some(line) = next_nonempty_line(&mut lines) {
        let trimmed = line.trim();
        if trimmed == "//" {
            saw_terminator = true;
            break;
        }

        let match_tokens: Vec<_> = trimmed.split_whitespace().collect();
        if match_tokens.len() < alphabet.len() + 1 {
            return Err(format!(
                "Match emission line for state '{}' does not have enough values",
                match_tokens.get(0).unwrap_or(&"?")
            ));
        }
        let label = match_tokens[0].to_string();
        let match_emissions =
            parse_value_slice(&match_tokens[1..=alphabet.len()], &label, "match")?;
        let annotation = match_tokens[alphabet.len() + 1..]
            .iter()
            .map(|s| (*s).to_string())
            .collect::<Vec<_>>();

        let insert_line = next_nonempty_line(&mut lines)
            .ok_or_else(|| format!("Missing insert emission line for state '{}'", label))?;
        let insert_trimmed = insert_line.trim();
        if insert_trimmed == "//" {
            return Err(format!(
                "Unexpected end of record after match emissions for state '{}'",
                label
            ));
        }
        let insert_tokens: Vec<_> = insert_trimmed.split_whitespace().collect();
        if insert_tokens.len() < alphabet.len() {
            return Err(format!(
                "Insert emission line for state '{}' does not have enough values",
                label
            ));
        }
        let insert_emissions =
            parse_value_slice(&insert_tokens[..alphabet.len()], &label, "insert")?;

        let transition_line = next_nonempty_line(&mut lines)
            .ok_or_else(|| format!("Missing transition line for state '{}'", label))?;
        let transition_trimmed = transition_line.trim();
        if transition_trimmed == "//" {
            return Err(format!(
                "Unexpected end of record after insert emissions for state '{}'",
                label
            ));
        }
        let transition_tokens: Vec<_> = transition_trimmed.split_whitespace().collect();
        if transition_tokens.len() < transition_order.len() {
            return Err(format!(
                "Transition line for state '{}' does not have enough values",
                label
            ));
        }
        let transitions = parse_value_slice(
            &transition_tokens[..transition_order.len()],
            &label,
            "transition",
        )?;

        states.push(HmmState {
            label,
            match_emissions,
            annotation,
            insert_emissions,
            transitions,
        });
    }

    if !saw_terminator {
        return Err("HMM record is missing the '//' terminator".into());
    }
    if states.is_empty() {
        return Err("No HMM states were parsed".into());
    }

    Ok(Hmm {
        format_line,
        metadata,
        alphabet,
        transition_order,
        states,
    })
}

fn next_nonempty_line<'a, I>(lines: &mut Peekable<I>) -> Option<&'a str>
where
    I: Iterator<Item = &'a str>,
{
    while let Some(line) = lines.next() {
        if !line.trim().is_empty() {
            return Some(line);
        }
    }
    None
}

fn parse_value_slice(
    tokens: &[&str],
    label: &str,
    field: &str,
) -> Result<Vec<Option<f64>>, String> {
    tokens
        .iter()
        .map(|token| {
            parse_value(token)
                .map_err(|err| format!("{} line for state '{}': {}", field, label, err))
        })
        .collect()
}

fn parse_value(token: &str) -> Result<Option<f64>, String> {
    if token == "*" {
        return Ok(None);
    }
    token
        .parse::<f64>()
        .map(Some)
        .map_err(|_| format!("invalid numeric value '{}'", token))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_example_hmm() {
        let text = include_str!("../static/hmm-example.hmm");
        let hmm = parse_hmm(text).expect("failed to parse HMM");
        assert_eq!(hmm.format_line, "HMMER3/f [3.3.2 | Nov 2020]");
        assert!(hmm
            .metadata
            .iter()
            .any(|field| field.key == "NAME" && field.value == "SPHERE-III.008_786"));
        assert_eq!(hmm.alphabet.len(), 20);
        assert_eq!(hmm.transition_order.len(), 7);
        assert_eq!(hmm.states.len(), 41);
        assert_eq!(hmm.states[0].label, "COMPO");
        assert_eq!(hmm.states[1].label, "1");
        assert_eq!(hmm.states[1].annotation.get(0), Some(&"1".to_string()));
    }
}
