import web_bio_tools
from Bio.Align import PairwiseAligner, substitution_matrices


def test_sw_blosum62_score():
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    result = web_bio_tools.smith_waterman_blosum62(seq1, seq2, -10.0, -0.5)
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    ref = aligner.align(seq1, seq2)[0]
    assert abs(result.score - ref.score) < 1e-6


def test_nw_globalms():
    seq1 = "GATTACA"
    seq2 = "GATTACA"
    result = web_bio_tools.needleman_wunsch(seq1, seq2)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    ref = aligner.align(seq1, seq2)[0]
    assert result.aligned_seq1 == ref[0]
    assert result.aligned_seq2 == ref[1]
    assert abs(result.score - ref.score) < 1e-6
