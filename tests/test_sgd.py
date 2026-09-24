from ptm_lollipop.sgd import normalize_ptm_type


def test_normalize_ptm_type():
    assert normalize_ptm_type("phosphorylated residue") == "phosphorylation"
    assert normalize_ptm_type("SUMOylation site") == "sumoylation"
    assert normalize_ptm_type("unknown modification") == "other"

