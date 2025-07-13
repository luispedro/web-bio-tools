import web_bio_tools
from Bio.Align import PairwiseAligner, substitution_matrices
from hypothesis import given, settings, assume
from hypothesis import strategies as st


def test_sw_blosum62_score():
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    result = web_bio_tools.smith_waterman_blosum62(seq1, seq2, -10.0, -0.5)
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    ref = aligner.align(seq1, seq2)
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


aa = "ARNDCQEGHILKMFPSTWYV"


@given(
    seq1=st.text(alphabet=aa, min_size=5, max_size=30),
    seq2=st.text(alphabet=aa, min_size=5, max_size=30),
    gap_open=st.floats(min_value=-20, max_value=-1),
    gap_extend=st.floats(min_value=-2, max_value=-0.1),
)
@settings(max_examples=20, deadline=None)
def test_sw_blosum62_hypothesis(seq1, seq2, gap_open, gap_extend):
    result = web_bio_tools.smith_waterman_blosum62(seq1, seq2, gap_open, gap_extend)
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    ref = aligner.align(seq1, seq2)
    assert abs(result.score - ref.score) < 1e-6


dna = "ACGT"


@given(
    seq1=st.text(alphabet=dna, min_size=5, max_size=40),
    seq2=st.text(alphabet=dna, min_size=5, max_size=40),
    match_score=st.integers(min_value=1, max_value=5),
    mismatch_penalty=st.integers(min_value=-5, max_value=-1),
    gap_open=st.floats(min_value=-5, max_value=-1),
    gap_extend=st.floats(min_value=-2, max_value=-0.1),
)
@settings(max_examples=20, deadline=None)
def test_nw_globalms_hypothesis(seq1, seq2, match_score, mismatch_penalty, gap_open, gap_extend):
    assume(gap_open <= gap_extend)
    result = web_bio_tools.needleman_wunsch_custom(
        seq1, seq2, match_score, mismatch_penalty, gap_open, gap_extend
    )
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_penalty
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    ref = aligner.align(seq1, seq2)
    assert abs(result.score - ref.score) < 1e-6
