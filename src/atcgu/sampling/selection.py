"""Conformer selection from a CREST ensemble: energy window or radius-of-gyration bins.
No jobs are submitted.
"""
import numpy as np

def read_multi_xyz(p):
    """CREST ensemble file -> [(energy_hartree, [sym], Nx3 array), ...]"""
    lines = open(p).read().splitlines(); out = []; i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1; continue
        n = int(lines[i].strip())
        if n < 1 or len(lines) < i + n + 2 or not lines[i+1].split():
            raise ValueError("Incomplete CREST frame or missing energy")
        e = float(lines[i+1].split()[0])
        sym = []; X = []
        for ln in lines[i+2:i+2+n]:
            q = ln.split()
            if len(q) != 4:
                raise ValueError("Malformed CREST coordinate row")
            sym.append(q[0]); X.append([float(v) for v in q[1:4]])
        out.append((e, sym, np.array(X))); i += 2 + n
    return out

def kabsch_rmsd(P, Q, mask):
    P = P[mask] - P[mask].mean(0); Q = Q[mask] - Q[mask].mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    return float(np.sqrt(((P @ (V @ np.diag([1, 1, d]) @ W) - Q) ** 2).sum() / len(P)))

def radius_of_gyration(X, mask):
    Y = X[mask]
    return float(np.sqrt(((Y - Y.mean(0)) ** 2).sum(1).mean()))

def select_energy(confs, rel, mask, pol, rmsd_cutoff):
    """Frames in file order until the first one above the window; RMSD-deduplicated, capped at n_keep."""
    keep = []
    for k in range(len(confs)):
        if rel[k] > pol["window"]:
            break
        if any(kabsch_rmsd(confs[k][2], confs[j][2], mask) < rmsd_cutoff for j in keep):
            continue
        keep.append(k)
        if len(keep) >= pol["n_keep"]:
            break
    return keep

def select_rg_stratified(confs, rel, mask, pol, rg_seed, rmsd_cutoff):
    """Bin the ensemble by radius of gyration, take the lowest per_bin in each bin.

    The upper bin edge is extended to the seed's R_g so the most extended CREST geometries
    -- which are the ones that can rebut the ALPB compaction -- always fall inside a bin
    rather than off the end of the histogram.
    """
    idx_ok = np.where(rel <= pol["window"])[0]
    if len(idx_ok) == 0:
        return []
    R = np.array([radius_of_gyration(confs[k][2], mask) for k in idx_ok])
    lo, hi = R.min(), max(R.max(), rg_seed)
    edges = np.linspace(lo, hi, pol["n_bins"] + 1)
    edges[-1] += 1e-9
    keep = []
    for b in range(pol["n_bins"]):
        inb = idx_ok[(R >= edges[b]) & (R < edges[b+1])]
        inb = inb[np.argsort(rel[inb])]
        got = 0
        for k in inb:
            if any(kabsch_rmsd(confs[k][2], confs[j][2], mask) < rmsd_cutoff for j in keep):
                continue
            keep.append(k); got += 1
            if got >= pol["per_bin"]:
                break
    return sorted(keep, key=lambda k: rel[k])
