import cocoatree
import itertools


def test_ghost():
    print(cocoatree.__file__)
    assert True


def test_perform_sca():
    data = cocoatree.datasets.load_rhomboid_proteases()
    sequences_id = data["sequence_ids"][:300]
    sequences = data["alignment"][:300]

    coevol_metrics = ["SCA", "MI", "NMI"]
    corrections = [None, "entropy", "APC"]

    for coevol_met, corr in itertools.product(coevol_metrics, corrections):
        cocoatree.perform_sca(
            sequences_id, sequences, n_components=2,
            coevolution_metric=coevol_met,
            correction=corr
            )
    assert True
