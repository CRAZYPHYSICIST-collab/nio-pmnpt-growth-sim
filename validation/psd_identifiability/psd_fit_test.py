"""
Identifiability test: can the surface-diffusion coefficient (b = kappa/F) and its
activation energy (Ed) be recovered from noisy AFM-like PSDs? No PINN is used; the
analytic solution of the linear model is fitted by least squares.

Model (lengths in nm):
  omega(k)  = a k^2 + b k^4
  S(k,tau)  = Omega/(2 omega) (1 - exp(-2 omega tau)) + S0(k) exp(-2 omega tau)
  S0(k)     = pi xi^2 W0^2 exp(-(k xi)^2 / 4)
Noise: radially averaged PSD ~ scaled chi-square, m = max(2, int(k L)) * 3
"""
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

kB = 8.617e-5
W0, XI, OMEGA, A_TRUE, ED_TRUE = 0.2, 20.0, 0.01822, 0.04, 0.50
T = np.array([300.0, 400.0, 500.0, 573.0])
TAU = np.array([10, 30, 50, 70, 90, 120], float)
K = np.logspace(np.log10(0.006), np.log10(0.6), 60)
L = 1024.0
M = np.maximum(2, (K * L).astype(int)) * 3          # degrees of freedom for each k
B0 = 0.01 / np.exp(-ED_TRUE / (kB * 300.0))          # b(300 K) = 0.01 nm^3
B_TRUE = B0 * np.exp(-ED_TRUE / (kB * T))
N_SEEDS = 20


def S0(k):
    return np.pi * XI**2 * W0**2 * np.exp(-(k * XI) ** 2 / 4)


def S_model(k, tau, a, b, Om, s0):
    w = a * k**2 + b * k**4
    e = np.exp(-2 * w * tau)
    return Om / (2 * w) * (1 - e) + s0 * e


def measure(S, rng):
    return S if rng is None else S * rng.gamma(M, 1.0 / M)


def fit_one(b_true, rng):
    """One temperature: measure at 6 thicknesses, fit a, b, Omega together."""
    s0m = measure(S0(K), rng)                         # measured (noisy) substrate PSD
    data = np.array([measure(S_model(K, t, A_TRUE, b_true, OMEGA, S0(K)), rng) for t in TAU])
    w8 = np.sqrt(M)

    def res(p):
        a, b, Om = np.exp(p)
        return np.concatenate([w8 * (np.log(data[i]) - np.log(S_model(K, t, a, b, Om, s0m)))
                               for i, t in enumerate(TAU)])
    # several starting points to avoid local minima
    best = None
    for b_start in [1e-3, 1.0, 1e2]:
        r = least_squares(res, np.log([0.01, b_start, 0.01]), x_scale="jac", max_nfev=5000)
        if best is None or r.cost < best.cost:
            best = r
    return np.exp(best.x), data, s0m


def arrhenius(b):
    x = 1 / (kB * T)
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, np.log(b), rcond=None)
    dof = len(x) - 2
    s2 = res[0] / dof if res.size else 0.0
    se = np.sqrt(s2 / np.sum((x - x.mean()) ** 2))
    return -coef[0], se


# ---- 1) noise-free check ----
b_clean = np.array([fit_one(bt, None)[0][1] for bt in B_TRUE])
a_clean = np.array([fit_one(bt, None)[0][0] for bt in B_TRUE])
assert np.all(np.abs(b_clean / B_TRUE - 1) < 1e-4), b_clean / B_TRUE
assert np.all(np.abs(a_clean / A_TRUE - 1) < 1e-4), a_clean / A_TRUE
print("Noise-free check passed: b and a recovered exactly.\n")

# ---- 2) 20 noise seeds ----
bf = np.zeros((N_SEEDS, len(T)))
af = np.zeros_like(bf)
Ed = np.zeros(N_SEEDS)
keep = {}
for s in range(N_SEEDS):
    rng = np.random.default_rng(1000 + s)
    for j, bt in enumerate(B_TRUE):
        p, data, s0m = fit_one(bt, rng)
        af[s, j], bf[s, j] = p[0], p[1]
        if s == 0:
            keep[j] = (p, data, s0m)
    Ed[s] = arrhenius(bf[s])[0]

print(f"{'T (K)':>6} {'b_true':>10} {'b_fit (mean ± SD)':>24} {'error %':>8} | "
      f"{'a_fit (mean ± SD)':>22} {'error %':>8}")
for j in range(len(T)):
    eb = 100 * (bf[:, j].mean() / B_TRUE[j] - 1)
    ea = 100 * (af[:, j].mean() / A_TRUE - 1)
    print(f"{T[j]:6.0f} {B_TRUE[j]:10.4g} {bf[:, j].mean():11.4g} ± {bf[:, j].std(ddof=1):<10.3g} "
          f"{eb:8.1f} | {af[:, j].mean():9.4g} ± {af[:, j].std(ddof=1):<10.3g} {ea:8.1f}")
    rel_sd = 100 * bf[:, j].std(ddof=1) / B_TRUE[j]
    if abs(eb) > 20 or rel_sd > 20:
        print(f"   WARNING: at T = {T[j]:.0f} K the bias/scatter of b exceeds 20 % "
              f"(mean error {eb:.1f} %, SD {rel_sd:.1f} %)")

print(f"\nEd true = {ED_TRUE:.3f} eV   Ed fit = {Ed.mean():.3f} ± {Ed.std(ddof=1):.3f} eV "
      f"(bias {100 * (Ed.mean() / ED_TRUE - 1):+.1f} %)")

# ---- 3) example figure ----
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
for ax, j in zip(axs, [0, 3]):
    p, data, s0m = keep[j]
    for i, t in enumerate(TAU):
        c = plt.cm.viridis(i / (len(TAU) - 1))
        ax.loglog(K, data[i], ".", color=c, ms=3)
        ax.loglog(K, S_model(K, t, *p, s0m), "-", color=c, lw=1.2, label=f"{t:.0f} nm")
    ax.set_title(f"T = {T[j]:.0f} K  (b true {B_TRUE[j]:.3g}, fit {p[1]:.3g} nm³)")
    ax.set_xlabel("k (nm⁻¹)")
axs[0].set_ylabel("PSD (nm⁴)")
axs[1].legend(title="thickness", fontsize=8)
fig.tight_layout()
fig.savefig("fit_example.png", dpi=150)
print("\nFigure: fit_example.png")
