#!/usr/bin/env bash
#SBATCH --job-name=cats3d_kernel
#SBATCH --output=/scratch8/Wenya/cats3d_kernel_%j.log
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

dlm_root=${DLM_ROOT:-/home/ktas/go54bav/DLM}
kernel_out=${KERNEL_OUT:-/scratch8/Wenya/pipi_coulomb_q24_r72_nonuniform.root}

cd "$dlm_root"
test -d "$(dirname "$kernel_out")"
if [[ -e "$kernel_out" ]]; then
    echo "Refusing to overwrite existing kernel: $kernel_out" >&2
    exit 1
fi

export OMP_NUM_THREADS=1
make -C CATS_3D/farm/kernel run KERNEL_OUT="$kernel_out"
