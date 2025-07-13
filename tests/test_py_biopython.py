import web_bio_tools
from Bio import pairwise2
from Bio.Align import substitution_matrices


def test_sw_blosum62_score():
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    result = web_bio_tools.smith_waterman_blosum62(seq1, seq2, -10.0, -0.5)
    matrix = substitution_matrices.load("BLOSUM62")
    ref = pairwise2.align.localds(seq1, seq2, matrix, -10.0, -0.5)[0]
    assert abs(result.score - ref.score) < 1e-6


def test_nw_globalms():
    seq1 = "GATTACA"
    seq2 = "GATTACA"
    result = web_bio_tools.needleman_wunsch(seq1, seq2)
    ref = pairwise2.align.globalms(seq1, seq2, 2, -1, -1, -0.5)[0]
    assert result.aligned_seq1 == ref.seqA
    assert result.aligned_seq2 == ref.seqB
    assert abs(result.score - ref.score) < 1e-6
