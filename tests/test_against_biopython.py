import pytest
import web_bio_tools
from Bio import pairwise2
from Bio.Align import substitution_matrices


def test_sw_blosum62_vs_biopython():
    mat = substitution_matrices.load("BLOSUM62")
    biopy = pairwise2.align.localds("ACACACTA", "AGCACACA", mat, -1.0, -0.5)[0]

    r = web_bio_tools.py_smith_waterman_blosum62(
        "ACACACTA", "AGCACACA", gap_open=-1.0, gap_extend=-0.5
    )
    assert r["aligned_seq1"] == biopy.seqA
    assert r["aligned_seq2"] == biopy.seqB
    assert r["score"] == pytest.approx(biopy.score)


def test_nw_blosum62_vs_biopython():
    mat = substitution_matrices.load("BLOSUM62")
    biopy = pairwise2.align.globalds("GATTACA", "GATTACA", mat, -1.0, -0.5)[0]

    r = web_bio_tools.py_needleman_wunsch_blosum62(
        "GATTACA", "GATTACA", gap_open=-1.0, gap_extend=-0.5
    )
    assert r["aligned_seq1"] == biopy.seqA
    assert r["aligned_seq2"] == biopy.seqB
    assert r["score"] == pytest.approx(biopy.score)
