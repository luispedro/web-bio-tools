use serde::Serialize;

/// Mapping for nucleotide characters to bitmasks of canonical bases.
/// The canonical bases use the order A, C, G, T (U maps to T).
fn base_mask(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0b0001),
        b'C' | b'c' => Some(0b0010),
        b'G' | b'g' => Some(0b0100),
        b'T' | b't' | b'U' | b'u' => Some(0b1000),
        b'R' | b'r' => Some(0b0101),
        b'Y' | b'y' => Some(0b1010),
        b'S' | b's' => Some(0b0110),
        b'W' | b'w' => Some(0b1001),
        b'K' | b'k' => Some(0b1100),
        b'M' | b'm' => Some(0b0011),
        b'B' | b'b' => Some(0b1110),
        b'D' | b'd' => Some(0b1101),
        b'H' | b'h' => Some(0b1011),
        b'V' | b'v' => Some(0b0111),
        b'N' | b'n' | b'X' | b'x' => Some(0b1111),
        b'-' => Some(0),
        _ => None,
    }
}

fn canonical_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

fn codon_index(first: usize, second: usize, third: usize) -> usize {
    (first << 4) | (second << 2) | third
}

/// Reverse complement nucleotide text into the provided buffer.
pub fn rev_compl_to(src: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(src.len());
    for &b in src.iter().rev() {
        dst.push(match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'R' | b'r' => b'Y',
            b'Y' | b'y' => b'R',
            b'S' | b's' => b'S',
            b'W' | b'w' => b'W',
            b'K' | b'k' => b'M',
            b'M' | b'm' => b'K',
            b'B' | b'b' => b'V',
            b'D' | b'd' => b'H',
            b'H' | b'h' => b'D',
            b'V' | b'v' => b'B',
            b'N' | b'n' | b'X' | b'x' => b'N',
            b'-' => b'-',
            other => other,
        });
    }
}

/// Encapsulates an amino acid lookup table for canonical codons.
#[derive(Debug, Clone)]
pub struct CodonEncoder {
    table: [u8; 64],
}

impl CodonEncoder {
    pub fn mk_encoder() -> Self {
        CodonEncoder {
            table: build_table(),
        }
    }

    pub fn translate_triplet(&self, codon: &[u8]) -> (u8, bool) {
        if codon.len() != 3 {
            return (b'X', true);
        }

        let mut masks = [0u8; 3];
        for (idx, &nt) in codon.iter().enumerate() {
            match base_mask(nt) {
                Some(mask) => masks[idx] = mask,
                None => return (b'X', true),
            }
        }

        if masks.iter().any(|&m| m == 0) {
            return (b'X', true);
        }

        let mut aa: Option<u8> = None;

        for first in 0..4 {
            if masks[0] & (1 << first) == 0 {
                continue;
            }
            for second in 0..4 {
                if masks[1] & (1 << second) == 0 {
                    continue;
                }
                for third in 0..4 {
                    if masks[2] & (1 << third) == 0 {
                        continue;
                    }
                    let idx = codon_index(first, second, third);
                    let current = self.table[idx];
                    match aa {
                        None => aa = Some(current),
                        Some(existing) if existing == current => {}
                        _ => {
                            return (b'X', true);
                        }
                    }
                }
            }
        }

        let aa = aa.unwrap_or(b'X');
        (aa, aa == b'X')
    }
}

pub fn build_table() -> [u8; 64] {
    let mut table = [b'X'; 64];
    for (codon, aa) in STANDARD_TABLE.iter() {
        let chars: Vec<u8> = codon.bytes().collect();
        let idx = codon_index(
            canonical_index(chars[0]).unwrap(),
            canonical_index(chars[1]).unwrap(),
            canonical_index(chars[2]).unwrap(),
        );
        table[idx] = *aa as u8;
    }
    table
}

const STANDARD_TABLE: [(&str, char); 64] = [
    ("TTT", 'F'),
    ("TTC", 'F'),
    ("TTA", 'L'),
    ("TTG", 'L'),
    ("TCT", 'S'),
    ("TCC", 'S'),
    ("TCA", 'S'),
    ("TCG", 'S'),
    ("TAT", 'Y'),
    ("TAC", 'Y'),
    ("TAA", '*'),
    ("TAG", '*'),
    ("TGT", 'C'),
    ("TGC", 'C'),
    ("TGA", '*'),
    ("TGG", 'W'),
    ("CTT", 'L'),
    ("CTC", 'L'),
    ("CTA", 'L'),
    ("CTG", 'L'),
    ("CCT", 'P'),
    ("CCC", 'P'),
    ("CCA", 'P'),
    ("CCG", 'P'),
    ("CAT", 'H'),
    ("CAC", 'H'),
    ("CAA", 'Q'),
    ("CAG", 'Q'),
    ("CGT", 'R'),
    ("CGC", 'R'),
    ("CGA", 'R'),
    ("CGG", 'R'),
    ("ATT", 'I'),
    ("ATC", 'I'),
    ("ATA", 'I'),
    ("ATG", 'M'),
    ("ACT", 'T'),
    ("ACC", 'T'),
    ("ACA", 'T'),
    ("ACG", 'T'),
    ("AAT", 'N'),
    ("AAC", 'N'),
    ("AAA", 'K'),
    ("AAG", 'K'),
    ("AGT", 'S'),
    ("AGC", 'S'),
    ("AGA", 'R'),
    ("AGG", 'R'),
    ("GTT", 'V'),
    ("GTC", 'V'),
    ("GTA", 'V'),
    ("GTG", 'V'),
    ("GCT", 'A'),
    ("GCC", 'A'),
    ("GCA", 'A'),
    ("GCG", 'A'),
    ("GAT", 'D'),
    ("GAC", 'D'),
    ("GAA", 'E'),
    ("GAG", 'E'),
    ("GGT", 'G'),
    ("GGC", 'G'),
    ("GGA", 'G'),
    ("GGG", 'G'),
];

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct FrameTranslation {
    pub frame: i8,
    pub amino_acids: String,
    pub stops: Vec<usize>,
    pub ambiguous: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct TranslationSummary {
    pub frames: Vec<FrameTranslation>,
}

pub fn translate_frame_internal(
    encoder: &CodonEncoder,
    sequence: &str,
    frame: i8,
    stop_at_first: bool,
) -> Result<FrameTranslation, String> {
    let (is_reverse, offset) = match frame {
        0 | 1 | 2 => (false, frame as usize),
        -1 | -2 | -3 => (true, (-frame - 1) as usize),
        _ => return Err(format!("Invalid frame: {}", frame)),
    };

    let seq_bytes = sequence.as_bytes();
    let mut buffer;
    let working = if is_reverse {
        buffer = Vec::with_capacity(seq_bytes.len());
        rev_compl_to(seq_bytes, &mut buffer);
        &buffer[offset..]
    } else {
        if seq_bytes.len() < offset {
            return Ok(FrameTranslation {
                frame,
                amino_acids: String::new(),
                stops: Vec::new(),
                ambiguous: Vec::new(),
            });
        }
        &seq_bytes[offset..]
    };

    let mut amino_acids = Vec::with_capacity(working.len() / 3);
    let mut stops = Vec::new();
    let mut ambiguous_positions = Vec::new();

    let mut index = 0;
    while index + 3 <= working.len() {
        let codon = &working[index..index + 3];
        let (aa, ambiguous) = encoder.translate_triplet(codon);
        let aa_char = aa as char;
        let aa_index = amino_acids.len();

        if ambiguous {
            ambiguous_positions.push(aa_index);
        }

        amino_acids.push(aa);

        if aa_char == '*' {
            stops.push(aa_index);
            if stop_at_first {
                break;
            }
        }

        index += 3;
    }

    let amino_acids = String::from_utf8(amino_acids).unwrap_or_default();

    Ok(FrameTranslation {
        frame,
        amino_acids,
        stops,
        ambiguous: ambiguous_positions,
    })
}

pub fn translate_all_frames_internal(
    encoder: &CodonEncoder,
    sequence: &str,
    stop_at_first: bool,
) -> TranslationSummary {
    let mut frames = Vec::with_capacity(6);
    for &frame in &[0, 1, 2, -1, -2, -3] {
        match translate_frame_internal(encoder, sequence, frame, stop_at_first) {
            Ok(result) => frames.push(result),
            Err(_) => {}
        }
    }
    TranslationSummary { frames }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_codon_translation() {
        let encoder = CodonEncoder::mk_encoder();
        assert_eq!(encoder.translate_triplet(b"ATG"), (b'M', false));
        assert_eq!(encoder.translate_triplet(b"TAA"), (b'*', false));
    }

    #[test]
    fn test_ambiguous_codon_translation() {
        let encoder = CodonEncoder::mk_encoder();
        assert_eq!(encoder.translate_triplet(b"ATN"), (b'X', true));
        assert_eq!(encoder.translate_triplet(b"NNN"), (b'X', true));
    }

    #[test]
    fn test_rev_complement() {
        let mut buffer = Vec::new();
        rev_compl_to(b"ATGC", &mut buffer);
        assert_eq!(buffer, b"GCAT");
    }

    #[test]
    fn test_translate_frame_forward() {
        let encoder = CodonEncoder::mk_encoder();
        let result = translate_frame_internal(&encoder, "ATGGCC", 0, false).unwrap();
        assert_eq!(result.amino_acids, "MA");
        assert!(result.stops.is_empty());
        assert!(result.ambiguous.is_empty());
    }

    #[test]
    fn test_translate_frame_reverse() {
        let encoder = CodonEncoder::mk_encoder();
        let result = translate_frame_internal(&encoder, "ATGGCC", -1, false).unwrap();
        assert_eq!(result.amino_acids, "GH");
    }

    #[test]
    fn test_translate_frame_stop_at_first() {
        let encoder = CodonEncoder::mk_encoder();
        let result = translate_frame_internal(&encoder, "ATGTAATTT", 0, true).unwrap();
        assert_eq!(result.amino_acids, "M*");
        assert_eq!(result.stops, vec![1]);
    }

    #[test]
    fn test_translate_frame_with_ambiguity_marks_position() {
        let encoder = CodonEncoder::mk_encoder();
        let result = translate_frame_internal(&encoder, "ATNCCC", 0, false).unwrap();
        assert_eq!(result.amino_acids, "XP");
        assert_eq!(result.ambiguous, vec![0]);
    }

    #[test]
    fn test_translate_all_frames() {
        let encoder = CodonEncoder::mk_encoder();
        let summary = translate_all_frames_internal(&encoder, "ATGGCC", false);
        assert_eq!(summary.frames.len(), 6);
        assert_eq!(summary.frames[0].frame, 0);
    }
}
