use crate::fna2faa::{self, CodonEncoder, FrameTranslation, TranslationSummary};

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

fn normalize_frame(frame: i8) -> Result<i8, String> {
    match frame {
        1..=3 => Ok(frame - 1),
        -3..=-1 => Ok(frame),
        _ => Err(format!(
            "Frame must be between -3 and 3 (excluding 0). Got {}",
            frame
        )),
    }
}

fn translate_frame_with_encoder(
    encoder: &CodonEncoder,
    sanitized: &str,
    frame: i8,
    stop_at_first_stop: bool,
) -> Result<FrameTranslation, String> {
    let internal_frame = normalize_frame(frame)?;
    fna2faa::translate_frame_internal(encoder, sanitized, internal_frame, stop_at_first_stop)
}

pub fn translate_frame(seq: &str, frame: i8, stop_at_first_stop: bool) -> Result<String, String> {
    let sanitized = sanitize_sequence(seq)?;
    let encoder = CodonEncoder::mk_encoder();
    let translation =
        translate_frame_with_encoder(&encoder, &sanitized, frame, stop_at_first_stop)?;
    Ok(translation.amino_acids)
}

pub fn translate_all_frames(
    seq: &str,
    stop_at_first_stop: bool,
) -> Result<TranslationSummary, String> {
    let sanitized = sanitize_sequence(seq)?;
    let encoder = CodonEncoder::mk_encoder();
    Ok(fna2faa::translate_all_frames_internal(
        &encoder,
        &sanitized,
        stop_at_first_stop,
    ))
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
        assert_eq!(result, "M*");
    }

    #[test]
    fn handles_all_frames() {
        let summary = translate_all_frames("ATGGCC", false).unwrap();
        assert_eq!(summary.frames.len(), 6);
        assert_eq!(summary.frames[0].frame, 1);
        assert_eq!(summary.frames[0].amino_acids, "MA");
        assert_eq!(summary.frames[3].frame, -1);
        assert_eq!(summary.frames[3].amino_acids, "GH");
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

    #[test]
    fn translates_reverse_frame() {
        let result = translate_frame("ATGGCC", -1, false).unwrap();
        assert_eq!(result, "GH");
    }
}
