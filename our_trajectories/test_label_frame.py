"""Orchestration checks for label_frame.py's process(), run without a real
chargemol call: `label_vasp_dir` (the chargemol/DDEC6 mechanics, now living
in volumetric_tools) is monkeypatched to a fake that records its kwargs
instead of running anything real, so these exercise the part of the
pipeline this refactor actually touches -- formal-charges/config_type
threading and the CHARGEMOL_OK/CHARGEMOL_FAILED/gzip_all bookkeeping.

    conda run -n dft python test_label_frame.py
"""
from __future__ import annotations

import gzip
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

import label_frame
from formal_charges import formal_charges_for


def _make_case_dir(tmp_dir: Path, compound: str, system_tag: str) -> Path:
    # label_frame.py's process() derives `compound` from
    # case_dir.resolve().parents[2].name -- reproduce that depth here.
    case_dir = tmp_dir / compound / "1000K" / "seg00-eos1.00" / "f00000"
    case_dir.mkdir(parents=True)
    (case_dir / "INCAR").write_text(f"SYSTEM = {system_tag}\nENCUT = 520\n")
    # some other file, to check gzip_all() still runs over the rest of the dir
    (case_dir / "OUTCAR").write_text("dummy outcar contents\n")
    return case_dir


class _FakePreviewAtoms:
    def get_chemical_symbols(self):
        return ["Li", "P"]


def test_process_success_threads_formal_charges_and_config_type():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case_dir(Path(tmp), "Li-P", "Li-P_1000K_seg00-eos1.00_f0")

        recorded = {}

        def fake_label_vasp_dir(vasp_dir, **kwargs):
            recorded["vasp_dir"] = vasp_dir
            recorded.update(kwargs)
            return SimpleNamespace()

        orig_ase_read = label_frame.ase_read
        orig_label_vasp_dir = label_frame.label_vasp_dir
        label_frame.ase_read = lambda p: _FakePreviewAtoms()
        label_frame.label_vasp_dir = fake_label_vasp_dir
        try:
            ok = label_frame.process(case_dir)
        finally:
            label_frame.ase_read = orig_ase_read
            label_frame.label_vasp_dir = orig_label_vasp_dir

        assert ok is True
        assert recorded["extra_arrays"]["REF_formal_charges"].tolist() == [
            formal_charges_for("Li-P")["Li"],
            formal_charges_for("Li-P")["P"],
        ]
        assert recorded["extra_info"]["config_type"] == "Li-P_1000K_seg00-eos1.00_f0"
        assert recorded["delete_originals"] is True


def test_process_success_writes_ok_marker_and_gzips():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case_dir(Path(tmp), "Li-P", "Li-P_tag")

        orig_ase_read = label_frame.ase_read
        orig_label_vasp_dir = label_frame.label_vasp_dir
        label_frame.ase_read = lambda p: _FakePreviewAtoms()
        label_frame.label_vasp_dir = lambda vasp_dir, **kwargs: SimpleNamespace()
        try:
            ok = label_frame.process(case_dir)
        finally:
            label_frame.ase_read = orig_ase_read
            label_frame.label_vasp_dir = orig_label_vasp_dir

        assert ok is True
        # gzip_all() runs after the OK marker is touched, so it ends up
        # CHARGEMOL_OK.gz like every other non-.h5 file in the directory.
        assert (case_dir / "CHARGEMOL_OK.gz").exists()
        assert not (case_dir / "CHARGEMOL_FAILED.gz").exists()
        remaining = [f for f in case_dir.iterdir() if f.is_file()]
        assert remaining, "expected some files to remain"
        assert all(f.suffix == ".gz" for f in remaining), remaining


def test_process_failure_writes_failed_marker_with_traceback():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case_dir(Path(tmp), "Li-P", "Li-P_tag")

        def raising_label_vasp_dir(vasp_dir, **kwargs):
            raise RuntimeError("simulated chargemol failure")

        orig_ase_read = label_frame.ase_read
        orig_label_vasp_dir = label_frame.label_vasp_dir
        label_frame.ase_read = lambda p: _FakePreviewAtoms()
        label_frame.label_vasp_dir = raising_label_vasp_dir
        try:
            ok = label_frame.process(case_dir)
        finally:
            label_frame.ase_read = orig_ase_read
            label_frame.label_vasp_dir = orig_label_vasp_dir

        assert ok is False
        assert not (case_dir / "CHARGEMOL_OK.gz").exists()
        # gzip_all() runs after the traceback is written, so it ends up
        # CHARGEMOL_FAILED.gz like every other non-.h5 file in the directory.
        failed = case_dir / "CHARGEMOL_FAILED.gz"
        assert failed.exists()
        with gzip.open(failed, "rt") as fh:
            assert "simulated chargemol failure" in fh.read()
        remaining = [f for f in case_dir.iterdir() if f.is_file()]
        assert remaining
        assert all(f.suffix == ".gz" for f in remaining), remaining


def main() -> None:
    tests = [
        test_process_success_threads_formal_charges_and_config_type,
        test_process_success_writes_ok_marker_and_gzips,
        test_process_failure_writes_failed_marker_with_traceback,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"\n{len(tests)}/{len(tests)} passed")


if __name__ == "__main__":
    main()
