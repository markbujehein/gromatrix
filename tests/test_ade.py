"""Steps A, D, E: offline tests (metrics, cache). Network-free by design."""

import numpy as np
import pytest

from gromatrix.afdb import align, kabsch, lddt_ca, rmsd_after_fit, tm_score
from gromatrix.net import Cache


def _helix(n=60):
    t = np.arange(n) * np.radians(100)
    return np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)]


def _rot(deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def test_metrics_invariant_to_rigid_motion():
    ref = _helix()
    mod = ref @ _rot(73).T + [5, -3, 12]
    r, t = kabsch(mod, ref)
    assert np.allclose(mod @ r.T + t, ref, atol=1e-6)
    assert rmsd_after_fit(mod, ref) < 1e-6
    assert tm_score(mod, ref, len(ref))[0] == pytest.approx(1.0, abs=1e-6)
    assert lddt_ca(mod, ref) == pytest.approx(1.0)


def test_metrics_detect_a_displaced_domain():
    ref = _helix()
    mod = ref.copy()
    mod[30:] = mod[30:] @ _rot(60).T + [8, 0, 0]  # second half swings away
    tm = tm_score(mod, ref, len(ref))[0]
    assert 0.3 < tm < 0.7  # global: one half still superposes
    # local (15 Å) lDDT barely notices a rigid hinge motion - why both are reported
    assert 0.9 < lddt_ca(mod, ref) < 0.99
    assert rmsd_after_fit(mod[:30], ref[:30]) < 1e-6  # the rigid half is intact


def test_tm_score_normalisation_penalises_partial_coverage():
    ref = _helix()
    assert tm_score(ref[:30], ref[:30], 60)[0] == pytest.approx(0.5, abs=1e-6)


def test_alignment_handles_terminal_truncation():
    pairs = align("MKTAYIAKQR", "GSHMMKTAYIAKQRQISFVK")
    assert pairs[0] == (0, 4) and len(pairs) == 10


def test_cache_reuses_and_repairs(tmp_path):
    src = tmp_path / "src.cif"
    src.write_text("data_TEST\n")
    cache = Cache(tmp_path / "raw")
    url = src.as_uri()
    assert cache.fetch("T.cif", url)[1] == "downloaded"
    assert cache.fetch("T.cif", url)[1] == "cached"
    (tmp_path / "raw" / "T.cif").write_text("truncat")
    assert cache.fetch("T.cif", url)[1] == "repaired"
    assert (tmp_path / "raw" / "T.cif").read_text() == "data_TEST\n"
    assert Cache(tmp_path / "raw").valid("T.cif")  # persisted in cache.json
