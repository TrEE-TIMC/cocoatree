from cocoatree.datasets import load_S1A_serine_proteases


def test_load_S1A_serine_proteases():
    dataset = load_S1A_serine_proteases()
    assert "sequence_ids" in dataset.keys()
    assert "alignment" in dataset.keys()
