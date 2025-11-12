import web_bio_tools

import warnings
from Bio import BiopythonWarning
from Bio.Seq import Seq
from hypothesis import given, settings
import hypothesis as hp
from hypothesis import strategies as st

dna = "ACGT"


def translate_with_biopython(seq: str, frame: int, stop_at_first_stop: bool) -> str:
    normalized = seq.replace("u", "t").replace("U", "T").upper()
    seq_obj = Seq(normalized)
    if frame > 0:
        working = seq_obj[frame - 1:]
    else:
        working = seq_obj.reverse_complement()[(-frame) - 1 :]

    working_str = str(working)
    trimmed_len = (len(working_str) // 3) * 3
    if trimmed_len == 0:
        return ""

    trimmed = Seq(working_str[:trimmed_len])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        translated = str(trimmed.translate(table=1, stop_symbol="*"))

    if stop_at_first_stop:
        stop_index = translated.find("*")
        if stop_index != -1:
            translated = translated[: stop_index + 1]

    return translated


@given(
    seq=st.text(alphabet=dna, min_size=3, max_size=120),
    frame=st.sampled_from([-3, -2, -1, 1, 2, 3]),
    stop_at_first_stop=st.booleans(),
)
@settings(
    max_examples=20,
    suppress_health_check=[hp.HealthCheck.data_too_large],
    deadline=None)
def test_translate_dna_frame_matches_biopython(seq, frame, stop_at_first_stop):
    ours = web_bio_tools.translate_dna_frame(seq, frame, stop_at_first_stop)
    expected = translate_with_biopython(seq, frame, stop_at_first_stop)
    assert ours == expected


def test_fna2faa():
    # https://gmsc.big-data-biology.org/cluster/GMSC10.90AA.283_000_000
    dna_seq = "ATGCACGGACACTCCCCGGACGTCACGACCACCACGGTGGACGTGGTCGCCCACGCGGGTTACCGCATCGGGGACCGCGTCCTGCGGGCCGCGAAGGTGACCGTGCTGGATCCTGAGAGCTGA"
    expected_aa_seq = "MHGHSPDVTTTTVDVVAHAGYRIGDRVLRAAKVTVLDPES*"
    assert web_bio_tools.translate_dna_frame(dna_seq, 1, True) == expected_aa_seq
