#!/usr/bin/env python3
"""Create local compatibility links after placing ROOT data in inputs/."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OTHER = ROOT / "other"

LINKS = {
    OTHER / "File": Path("../inputs/File"),
    OTHER / "OutputFiles": Path("../inputs/OutputFiles"),
    OTHER / "EPOS_data": Path("../inputs/EPOS_data"),
    OTHER / "farm_output": Path("../farm/fit/farm_output"),
    OTHER / "nuclearStructureDatafit": Path("../fit/nuclearStructureDatafit"),
    ROOT / "inputs/epos/EPOS_data": Path("../EPOS_data"),
    ROOT / "inputs/epos/resonance_PR.root": Path("../EPOS_data/resonance_PR.root"),
    ROOT / "inputs/epos/resonance_RR.root": Path("../EPOS_data/resonance_RR.root"),
}


def main() -> None:
    for link, target in LINKS.items():
        link.parent.mkdir(parents=True, exist_ok=True)
        if link.is_symlink():
            if link.readlink() != target:
                raise RuntimeError(f"Unexpected link target: {link} -> {link.readlink()}")
        elif link.exists():
            raise RuntimeError(f"Path exists and is not a link: {link}")
        else:
            link.symlink_to(target)
    kernels = ROOT / "inputs/kernels"
    kernels.mkdir(exist_ok=True)
    for name in (
        "sourceKernel_q10_r31_rmax20_coulomb_byQ.root",
        "adaptive_pipi_coulomb_q19_r34_byQ.root",
        "dataKernel_qData_r31_coulomb_byQ.root",
    ):
        link = kernels / name
        target = Path("../OutputFiles") / name
        if link.is_symlink():
            if link.readlink() != target:
                raise RuntimeError(f"Unexpected link target: {link}")
        elif link.exists():
            raise RuntimeError(f"Path exists and is not a link: {link}")
        else:
            link.symlink_to(target)
    print("Local data links ready")


if __name__ == "__main__":
    main()
