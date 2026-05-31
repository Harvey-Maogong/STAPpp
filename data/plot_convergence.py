import os
import numpy as np
import matplotlib.pyplot as plt


def find_outfile(basename):
    for ext in ['.dat.out', '.out', '_out.txt']:
        f = f"data/{basename}{ext}"
        if os.path.exists(f): return f
    import glob
    files = glob.glob(f"data/{basename}*.out")
    return files[0] if files else None


def extract_disp(filename, node_id, dof=2):
    if not filename or not os.path.exists(filename): return None
    with open(filename, 'r') as f:
        lines = f.readlines()
    in_disp = False
    for line in lines:
        if "D I S P L A C E M E N T S" in line:
            in_disp = True
            continue
        if in_disp and ("S T R E S S" in line or "S O L U T I O N" in line):
            break
        if in_disp:
            p = line.split()
            if len(p) >= 4 and p[0].isdigit() and int(p[0]) == node_id:
                return float(p[dof])
    return None


Ns = np.array([1, 2, 4, 8])
h = 2.0 / Ns
nodes = {1: 2, 2: 3, 4: 5, 8: 9}
exact = 0.004  # PL^3/(3EI) = 1*8/(3*1000*2/3) = 0.004

disps = []
for N in Ns:
    d = extract_disp(find_outfile(f"beam_short_N{N}"), nodes[N], dof=2)
    disps.append(abs(d) if d else np.nan)
    if d:
        print(
            f"N={N}: |disp|={abs(d):.6f}, error={abs(abs(d)-exact)/exact*100:.1f}%"
        )

disps = np.array(disps)
valid = ~np.isnan(disps)
h_v, d_v = h[valid], disps[valid]
err = np.abs(d_v - exact) / exact

p = np.polyfit(np.log(h_v), np.log(err), 1)[0] if len(h_v) >= 2 else 0

plt.figure(figsize=(8, 6))
plt.loglog(h_v,
           err,
           'ro-',
           linewidth=2,
           markersize=10,
           label=f'T3 Short Beam (p≈{p:.2f})')
plt.loglog(h_v, (h_v / h_v[0])**1 * err[0], 'k--', label='O(h)')
plt.loglog(h_v, (h_v / h_v[0])**2 * err[0], 'k:', label='O(h²)')
plt.xlabel('Element size $h$')
plt.ylabel('Relative error')
plt.title(
    'Convergence: Short Cantilever Beam ($L=2, H=2$)\nT3 Element, Plane Stress'
)
plt.legend()
plt.grid(True, which='both', ls='--', alpha=0.6)
plt.tight_layout()
plt.savefig('convergence_short.png', dpi=300)
plt.show()
print(f"\nObserved rate: {p:.2f}")
