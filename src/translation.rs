const VALID_NUCLEOTIDES: &[char] = &[
    'A', 'C', 'G', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N',
];

fn sanitize_sequence(seq: &str) -> Result<String, String> {
    let mut cleaned = String::with_capacity(seq.len());
    for ch in seq.chars() {
        if ch.is_whitespace() {
            continue;
        }
        if ch == '-' {
            // Ignore gap symbols that might be present in FASTA inputs.
            continue;
        }
        if ch.is_ascii_digit() {
            return Err(format!(
                "Invalid character '{}' found in nucleotide sequence",
                ch
            ));
        }
        if ch.is_ascii_alphabetic() {
            let upper = ch.to_ascii_uppercase();
            if VALID_NUCLEOTIDES.contains(&upper) {
                cleaned.push(match upper {
                    'U' => 'T',
                    _ => upper,
                });
            } else {
                return Err(format!(
                    "Unsupported nucleotide '{}' encountered in sequence",
                    ch
                ));
            }
            continue;
        }
        return Err(format!(
            "Invalid character '{}' found in nucleotide sequence",
            ch
        ));
    }
    Ok(cleaned)
}

fn translate_codon(codon: &str) -> char {
    if codon.chars().any(|c| !matches!(c, 'A' | 'C' | 'G' | 'T')) {
        return 'X';
    }
    match codon {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        "TAA" | "TAG" | "TGA" => '*',
        _ => 'X',
    }
}

fn translate_frame_internal(
    sanitized: &str,
    frame: usize,
    stop_at_first_stop: bool,
) -> Result<String, String> {
    if !(1..=3).contains(&frame) {
        return Err(format!(
            "Frame must be between 1 and 3 (inclusive). Got {}",
            frame
        ));
    }
    if sanitized.is_empty() {
        return Ok(String::new());
    }

    let mut translated = String::new();
    let mut index = frame - 1;

    while index + 3 <= sanitized.len() {
        let codon = &sanitized[index..index + 3];
        let aa = translate_codon(codon);
        if aa == '*' {
            if stop_at_first_stop {
                break;
            }
            translated.push('*');
        } else {
            translated.push(aa);
        }
        index += 3;
    }

    Ok(translated)
}

pub fn translate_frame(
    seq: &str,
    frame: usize,
    stop_at_first_stop: bool,
) -> Result<String, String> {
    let sanitized = sanitize_sequence(seq)?;
    translate_frame_internal(&sanitized, frame, stop_at_first_stop)
}

pub fn translate_all_frames(seq: &str, stop_at_first_stop: bool) -> Result<Vec<String>, String> {
    let sanitized = sanitize_sequence(seq)?;
    let mut out = Vec::with_capacity(3);
    for frame in 1..=3 {
        out.push(translate_frame_internal(
            &sanitized,
            frame,
            stop_at_first_stop,
        )?);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_simple_sequence() {
        let result = translate_frame("ATGGCC", 1, false).unwrap();
        assert_eq!(result, "MA");
    }

    #[test]
    fn respects_frame_offset() {
        let result = translate_frame("AATGGCC", 2, false).unwrap();
        assert_eq!(result, "MA");
    }

    #[test]
    fn stops_at_first_stop() {
        let result = translate_frame("ATGTAAATG", 1, true).unwrap();
        assert_eq!(result, "M");
    }

    #[test]
    fn handles_all_frames() {
        let result = translate_all_frames("ATGGCC", false).unwrap();
        assert_eq!(
            result,
            vec!["MA".to_string(), "W".to_string(), "G".to_string()]
        );
    }

    #[test]
    fn rejects_invalid_characters() {
        let err = translate_frame("ATG1CC", 1, false).unwrap_err();
        assert!(err.contains("Invalid character"));
    }

    #[test]
    fn rejects_invalid_frame() {
        let err = translate_frame("ATGGCC", 0, false).unwrap_err();
        assert!(err.contains("Frame"));
    }
}
