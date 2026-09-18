"""Run chargemol/DDEC6 on a finished single-point VASP directory (CHGCAR/AECCAR0/
AECCAR2/POTCAR present), then clean up: on success, delete AECCAR0 and CHG (AECCAR0
is redundant once DDEC charges exist; CHG is VASP's non-augmented pseudo density,
redundant with CHGCAR), lossily compress the remaining volumetric files
(CHGCAR/AECCAR1/AECCAR2) to zfp-compressed .h5 via `voltools compress
--delete-originals` (see ../../volumetric-tools; verified error is far below
DFT/AIMD noise, see compression_test/RESULTS.md and volumetric-tools/benchmarks/
vasp.md), and mark the directory OK; on failure, leave everything in place
(including the raw volumetric files, uncompressed) and mark it failed, since we
can't be sure they won't be needed to retry. Either way, gzip everything else
remaining (voltools' own .h5 output is already compressed, so it's skipped).

    conda run -n dft python run_chargemol.py <case_dir> [<case_dir> ...]
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

# ~/bin (where Chargemol_09_26_2017_linux_parallel lives) isn't always on PATH in
# non-interactive/batch shells -- make sure it is before pymatgen looks for it.
os.environ["PATH"] = f"{Path.home() / 'bin'}:{os.environ.get('PATH', '')}"

from pymatgen.command_line.chargemol_caller import ChargemolAnalysis  # noqa: E402

ATOMIC_DENSITIES = Path.home() / "chargemol" / "chargemol_09_26_2017" / "atomic_densities"
VOLTOOLS = Path.home() / "volumetric-tools" / ".venv" / "bin" / "voltools"


def compress_volumetric(case_dir: Path) -> None:
    """Compress CHGCAR/AECCAR1/AECCAR2 to zfp .h5 (default accuracy) and delete
    the raw files, logging the CLI's per-file ratio/error report."""
    with open(case_dir / "voltools.log", "wb") as log:
        subprocess.run(
            [str(VOLTOOLS), "compress", str(case_dir), "--delete-originals"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def gzip_all(case_dir: Path) -> None:
    for f in case_dir.iterdir():
        if f.is_file() and f.suffix not in (".gz", ".h5"):
            with open(f, "rb") as fin, gzip.open(f"{f}.gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)
            f.unlink()


def process(case_dir: Path) -> bool:
    os.environ.setdefault("OMP_NUM_THREADS", str(os.cpu_count() or 1))
    try:
        ChargemolAnalysis(path=str(case_dir), atomic_densities_path=str(ATOMIC_DENSITIES), run_chargemol=True)
    except Exception:
        (case_dir / "CHARGEMOL_FAILED").write_text(traceback.format_exc())
        gzip_all(case_dir)
        return False

    for name in ("AECCAR0", "CHG"):
        f = case_dir / name
        if f.exists():
            f.unlink()
    compress_volumetric(case_dir)
    (case_dir / "CHARGEMOL_OK").touch()
    gzip_all(case_dir)
    return True


def main() -> None:
    ok = 0
    for arg in sys.argv[1:]:
        case_dir = Path(arg)
        success = process(case_dir)
        ok += success
        print(f"{'OK  ' if success else 'FAIL'} {case_dir}")
    print(f"{ok}/{len(sys.argv) - 1} succeeded")


if __name__ == "__main__":
    main()
