"""Fold the nonuniform source-independent pi+pi+ kernel with a Gaussian source."""
import argparse
from array import array

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("kernel")
parser.add_argument("output")
parser.add_argument("--radii", nargs=3, type=float, metavar=("ROUT", "RSIDE", "RLONG"), required=True)
parser.add_argument("--lambda-effective", type=float, default=1.0)
args = parser.parse_args()

source = ROOT.TFile.Open(args.kernel)
if not source or source.IsZombie():
    raise RuntimeError("Cannot open kernel")
meta = source.Get("KernelMetadata")
tree = source.Get("KernelByQ")
if not meta or not tree or meta.GetEntries() != 1:
    raise RuntimeError("Missing kernel metadata or cells")
meta.GetEntry(0)
q_edges = np.asarray(meta.qOutEdgesMeV, dtype=np.float64)
r_edges = np.asarray(meta.rOutEdgesFm, dtype=np.float64)
if len(q_edges) != 20 or len(r_edges) != 35:
    raise RuntimeError("Unexpected grid size")

mid = 0.5 * (r_edges[:-1] + r_edges[1:])
width = np.diff(r_edges)
axes = [np.exp(-mid*mid/(4*radius*radius))*width for radius in args.radii]
weight = np.einsum("i,j,k->ijk", *axes).reshape(-1).astype(np.float64)
weight /= weight.sum()

out = ROOT.TFile(args.output, "RECREATE")
edges = array("d", q_edges)
h = ROOT.TH3D("hCkPositiveOctant", ";q_{out} (MeV/c);q_{side} (MeV/c);q_{long} (MeV/c)",
              len(q_edges)-1, edges, len(q_edges)-1, edges, len(q_edges)-1, edges)
for entry in tree:
    values = np.asarray(entry.wavefunction2, dtype=np.float64)
    if len(values) != len(weight):
        raise RuntimeError("Spatial kernel size mismatch")
    corr = 1.0 + args.lambda_effective * (np.dot(weight, values) - 1.0)
    h.SetBinContent(entry.iqOut+1, entry.iqSide+1, entry.iqLong+1, corr)
h.Write()
ROOT.TNamed("kernelPath", args.kernel).Write()
ROOT.TNamed("GaussianRadiiFm", ",".join(map(str, args.radii))).Write()
ROOT.TNamed("lambdaEffective", str(args.lambda_effective)).Write()
out.Close()
print(f"Folded {tree.GetEntries()} q cells into {args.output}")
