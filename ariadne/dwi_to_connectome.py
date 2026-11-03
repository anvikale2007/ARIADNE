"""Diffusion MRI (DWI) -> structural connectome, using DIPY.

    python dwi_to_connectome.py --dwi scan.nii.gz --bval scan.bval --bvec scan.bvec --out connectomes/subject01

Steps: brain mask -> diffusion tensor fit (FA, principal directions) -> deterministic tractography seeded in
white matter -> parcellate the brain into N regions -> count streamlines between every pair of regions.

!! LIMITATION !! Without a T1 image + anatomical atlas, regions here are *geometric* parcels (k-means on voxel
positions), not named anatomical structures (no "hippocampus"). With a registered atlas (e.g. AAL, Desikan) you
can pass it with --atlas (NIfTI of integer labels in DWI space) and the same code builds a labelled connectome.
"""
import argparse, json, os, warnings
import numpy as np
import nibabel as nib

warnings.filterwarnings("ignore")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dwi", required=True); ap.add_argument("--bval", required=True); ap.add_argument("--bvec", required=True)
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--parcels", type=int, default=70, help="number of geometric parcels (ignored with --atlas)")
    ap.add_argument("--atlas", default=None, help="optional integer-label NIfTI in DWI space")
    ap.add_argument("--max-b", type=float, default=1500, help="use volumes up to this b-value for the tensor fit")
    ap.add_argument("--fa-seed", type=float, default=0.20); ap.add_argument("--fa-stop", type=float, default=0.12)
    ap.add_argument("--min-len", type=float, default=20.0, help="drop streamlines shorter than this (mm)")
    ap.add_argument("--density", type=int, default=2, help="seeds per white-matter voxel (per axis)")
    a = ap.parse_args()

    from dipy.core.gradients import gradient_table
    from dipy.segment.mask import median_otsu
    from dipy.reconst.dti import TensorModel
    from dipy.data import get_sphere
    from dipy.direction import DeterministicMaximumDirectionGetter
    from dipy.tracking.local_tracking import LocalTracking
    from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
    from dipy.tracking import utils
    from dipy.tracking.streamline import Streamlines
    from sklearn.cluster import KMeans

    os.makedirs(a.out, exist_ok=True)
    img = nib.load(a.dwi); data = np.asarray(img.get_fdata(), dtype=np.float32); affine = img.affine
    bvals = np.loadtxt(a.bval); bvecs = np.loadtxt(a.bvec)
    if bvecs.shape[0] == 3 and bvecs.shape[1] != 3: bvecs = bvecs.T
    print(f"DWI {data.shape}, voxel {img.header.get_zooms()[:3]} mm, b-values {sorted(set(np.round(bvals,-2)))}")

    # 1. mask + tensor fit on low/medium b-values
    b0s = np.where(bvals < 50)[0]
    _, mask = median_otsu(data, vol_idx=b0s, median_radius=2, numpass=1, autocrop=False)
    sel = np.where(bvals <= a.max_b)[0]
    if (bvals[sel] > 50).sum() < 12:   # single high-b shell (e.g. Stanford b=2000): use all volumes
        sel = np.arange(len(bvals))
    gtab = gradient_table(bvals[sel], bvecs=bvecs[sel])
    tenfit = TensorModel(gtab).fit(data[..., sel], mask=mask)
    fa = np.nan_to_num(tenfit.fa)
    print(f"mask voxels: {int(mask.sum())} | white-matter-like voxels (FA>{a.fa_seed}): {int((fa > a.fa_seed).sum())}")

    # 2. tractography (deterministic, tensor-derived direction field)
    sphere = get_sphere(name="repulsion100")
    odf = tenfit.odf(sphere)
    pmf = np.clip(odf, 0, None)
    pmf[~mask] = 0
    dg = DeterministicMaximumDirectionGetter.from_pmf(pmf, max_angle=40.0, sphere=sphere)
    stop = ThresholdStoppingCriterion(fa, a.fa_stop)
    seeds = utils.seeds_from_mask(fa > a.fa_seed, affine, density=a.density)
    print(f"tracking from {len(seeds)} seeds ...")
    sl = Streamlines(LocalTracking(dg, stop, seeds, affine, step_size=1.0, max_cross=1))
    lens = np.array(list(utils.length(sl)))
    sl = Streamlines([s for s, L in zip(sl, lens) if L >= a.min_len])
    print(f"streamlines kept (>= {a.min_len} mm): {len(sl)}")
    if len(sl) < 200:
        raise SystemExit("Too few streamlines: check b-vectors orientation / FA thresholds.")

    # 3. parcellation
    if a.atlas:
        labels = np.asarray(nib.load(a.atlas).get_fdata()).astype(int)
    else:
        idx = np.argwhere(mask)
        km = KMeans(n_clusters=a.parcels, n_init=4, random_state=0).fit(idx * np.array(img.header.get_zooms()[:3]))
        labels = np.zeros(mask.shape, dtype=int); labels[tuple(idx.T)] = km.labels_ + 1
    ids = np.array(sorted(i for i in np.unique(labels) if i > 0)); n = len(ids)
    remap = np.zeros(labels.max() + 1, dtype=int); remap[ids] = np.arange(1, n + 1)
    labels = remap[labels]

    # 4. connectome: streamline endpoint counts between parcels
    # endpoints often stop just inside white matter, so dilate parcels over the brain mask for endpoint lookup
    from scipy import ndimage
    near = ndimage.distance_transform_edt(labels == 0, return_indices=True)
    filled = labels[tuple(near[1])]
    M = utils.connectivity_matrix(sl, affine, filled, symmetric=True, return_mapping=False)
    C = np.asarray(M, dtype=float)[1:n + 1, 1:n + 1]
    np.fill_diagonal(C, 0)
    centres = np.array([np.argwhere(labels == k).mean(0) for k in range(1, n + 1)])
    centres_mm = nib.affines.apply_affine(affine, centres)
    names = [f"P{k:02d}" for k in range(1, n + 1)]

    np.savetxt(os.path.join(a.out, "weights.csv"), C, delimiter=",", fmt="%.0f")
    np.savetxt(os.path.join(a.out, "labels.csv"), names, fmt="%s")
    np.savetxt(os.path.join(a.out, "centres_mm.csv"), centres_mm, delimiter=",", fmt="%.1f")
    nib.save(nib.Nifti1Image(labels.astype(np.int16), affine), os.path.join(a.out, "parcellation.nii.gz"))
    nib.save(nib.Nifti1Image(fa.astype(np.float32), affine), os.path.join(a.out, "fa.nii.gz"))

    deg = (C > 0).sum(1)
    info = {"regions": n, "streamlines": len(sl), "edges": int((np.triu(C, 1) > 0).sum()),
            "density": float((C > 0).sum() / (n * (n - 1))), "isolated_regions": int((deg == 0).sum()),
            "median_strength": float(np.median(C.sum(1))), "parcellation": "atlas" if a.atlas else "geometric k-means"}
    json.dump(info, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))

    # QC figure
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    z = fa.shape[2] // 2
    ax[0].imshow(np.rot90(fa[:, :, z]), cmap="gray"); ax[0].set_title(f"FA (axial slice {z})"); ax[0].axis("off")
    ax[1].imshow(np.log1p(C), cmap="magma"); ax[1].set_title("Connectome (log streamline count)")
    ax[1].set_xlabel("region"); ax[1].set_ylabel("region")
    d = np.linalg.norm(centres_mm[:, None] - centres_mm[None], axis=-1); iu = np.triu_indices(n, 1)
    ax[2].scatter(d[iu], np.log1p(C[iu]), s=4, alpha=.4); ax[2].set_xlabel("distance between regions (mm)")
    ax[2].set_ylabel("log(1+streamlines)"); ax[2].set_title("Sanity: connection strength falls with distance")
    plt.tight_layout(); plt.savefig(os.path.join(a.out, "qc.png"), dpi=110)
    print("Saved to", a.out)


if __name__ == "__main__":
    main()
