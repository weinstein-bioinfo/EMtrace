"""Smoke test of EMtrace's external tools on this machine.

Runs the tool wrappers and shell pipelines end to end on a tiny synthetic
data set, in a folder whose name contains spaces: chopper (pigz pipeline),
minimap2, samtools (sort | fixmate | sort, markdup, index), bcftools
(mpileup | call | norm) and bgzip/tabix. A known SNP is planted in the reads
and must be called.

On Windows this exercises the WSL route (see docs/windows_install.md): each
tool is reported as found natively or in WSL.

Usage (from an environment where emtrace is installed):

    python scripts/smoke_test_tools.py [work_dir]

Exits 0 when every check passes. The work folder is kept for inspection.
"""

from __future__ import annotations

import gzip
import random
import sys
import tempfile
import traceback
from pathlib import Path

from emtrace.utils import wsl
from emtrace.utils.shell import run_shell, shell_path
from emtrace.wrappers.bcftools import Bcftools
from emtrace.wrappers.chopper import Chopper
from emtrace.wrappers.minimap2 import Minimap2
from emtrace.wrappers.samtools import Samtools

TOOLS = ("minimap2", "samtools", "bcftools", "chopper", "pigz", "bgzip", "tabix")
SNP_POS = 2501  # 1-based position of the planted SNP


def report_tools() -> bool:
    print(f"Platform: {sys.platform}")
    if wsl.IS_WINDOWS:
        print(f"WSL available: {wsl.is_available()}")
        print(f"WSL distribution: {wsl.distro() or '(default)'}")
        print(f"WSL tool folders: {wsl.search_path() or '(none found)'}")
    ok = True
    for tool in TOOLS:
        found = wsl.locate_tool(tool)
        if found is None:
            print(f"  {tool:9s} NOT FOUND")
            ok = False
        else:
            path, via_wsl = found
            print(f"  {tool:9s} {'WSL   ' if via_wsl else 'native'} {path}")
    return ok


def make_data(root: Path) -> tuple[Path, Path]:
    random.seed(1)
    ref = "".join(random.choice("ACGT") for _ in range(5000))
    i = SNP_POS - 1
    alt = ref[:i] + ("A" if ref[i] != "A" else "C") + ref[i + 1:]
    ref_fa = root / "ref seq.fa"
    ref_fa.write_text(f">ref1\n{ref}\n", newline="\n")
    reads_fq = root / "reads in.fq.gz"
    with gzip.open(reads_fq, "wt", newline="\n") as fh:
        for n in range(300):
            s = random.randint(0, 4000)
            seq = alt[s:s + 900]
            fh.write(f"@r{n}\n{seq}\n+\n{'I' * len(seq)}\n")
    return ref_fa, reads_fq


def run_checks(root: Path) -> None:
    ref_fa, reads_fq = make_data(root)

    print("chopper ...", flush=True)
    fq = Chopper(num_threads=2).filter(
        reads_fq, root / "filtered reads.fq.gz",
        quality=10, min_length=500, log_file=root / "chopper log.txt",
    )
    assert fq.exists() and fq.stat().st_size > 0, "chopper produced no reads"

    print("minimap2 ...", flush=True)
    sam = Path(Minimap2(num_threads=2).align(
        ref_fa, fq, root / "aln out.sam", sam_output=True, n=-1,
    ))
    assert sam.stat().st_size > 0, "minimap2 produced an empty SAM"

    print("samtools sort | fixmate | sort, markdup, index ...", flush=True)
    st = Samtools(num_threads=2)
    bam = st.sort_fixmate_markdup(sam, root / "aln sorted.bam")
    mapped = st.count_mapped(bam)
    assert mapped > 250, f"only {mapped} reads mapped"
    run_shell(f"samtools faidx {shell_path(ref_fa)}", tool_name="samtools")

    print("bcftools mpileup | call | norm ...", flush=True)
    vcf = Bcftools(num_threads=2).mpileup_call(bam, ref_fa, root / "calls out.vcf")
    calls = [ln.split("\t") for ln in vcf.read_text().splitlines() if not ln.startswith("#")]
    assert any(c[1] == str(SNP_POS) for c in calls), f"SNP at {SNP_POS} not called: {calls}"

    print("bgzip, tabix ...", flush=True)
    run_shell(f"bgzip -f {shell_path(vcf)}", tool_name="bgzip")
    run_shell(f"tabix -f {shell_path(str(vcf) + '.gz')}", tool_name="tabix")
    assert Path(str(vcf) + ".gz.tbi").exists(), "tabix index missing"


def main() -> int:
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if base is not None:
        base.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="emtrace-smoke-", dir=base)) / "run dir with spaces"
    root.mkdir()

    if not report_tools():
        print("FAILED: missing tools (see docs/windows_install.md on Windows)")
        return 1
    print(f"Work folder: {root}")
    try:
        run_checks(root)
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        print("FAILED")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
