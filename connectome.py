"""Load the TVB 76-region template connectome and clean it."""
import os, zipfile
import numpy as np

def load_template():
    import tvb_data
    path = os.path.join(os.path.dirname(tvb_data.__file__), "connectivity", "connectivity_76.zip")
    z = zipfile.ZipFile(path)
    W = np.loadtxt(z.open("weights.txt"))
    labels = [l.split()[0] for l in z.open("centres.txt").read().decode().splitlines() if l.strip()]
    W = (W + W.T) / 2          # symmetrize
    np.fill_diagonal(W, 0)     # no self-connections
    W = W / W.max()            # scale to 0-1
    return W, labels

if __name__ == "__main__":
    C, labels = load_template()
    print("shape:", C.shape, "| symmetric:", np.allclose(C, C.T), "| diagonal zero:", np.allclose(np.diag(C), 0), "| max:", C.max())
    print("labels:", labels)
    hc = labels.index("rHC")
    top = np.argsort(C[hc])[::-1][:6]
    print("rHC strongest links:", [(labels[i], round(float(C[hc, i]), 3)) for i in top])
