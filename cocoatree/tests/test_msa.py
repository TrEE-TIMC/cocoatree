from cocoatree import msa
from cocoatree.datasets import load_S1A_serine_proteases


def test_filter_seq_id():
    sequences = load_S1A_serine_proteases()
    sequence_ids = sequences["sequence_ids"]
    sequences = sequences["alignment"]

    filtered_seq = msa.filter_seq_id(
        sequences,
        sequence_ids,
        sequence_ids[:100])
    assert len(filtered_seq) <= 100
