"""
Scenario: "drive-in" anneal of an implanted boron (B) layer in Si.
  Physics: dC/dt = D(T) d2C/dx2,   D(T) = D0 exp(-Ea / kT)
  Literature (intrinsic B in Si): D0 = 0.76 cm^2/s, Ea = 3.46 eV -> "true" values, used only to generate data
  Experiment: 1 h anneal at 4 temperatures, then a SIMS depth profile (5 % multiplicative noise)

  Dimensionless form: xi = x/L, tau = t/t_max, u = C/C_peak
      u_tau = kappa u_xixi,   kappa = D t_max / L^2
      IC: u(xi,0) = exp(-xi^2/a^2)  (a = w/L, implant profile)
      BC: u_xi(0,tau) = 0 (no out-diffusion),  u(1,tau) = 0 (deep region)
      Analytic: u = a/sqrt(a^2+4 kappa tau) * exp(-xi^2/(a^2+4 kappa tau))

Steps:
  1) Generate synthetic SIMS data
  2) Tune inverse-PINN hyperparameters with Bayesian optimization
     (objective: error on held-out points; the true D is NOT used)
  3) PINN with the best setting at each temperature -> D(T)
  4) Arrhenius fit -> Ea, D0; compare with literature
Requires: numpy, scipy, scikit-learn, torch (CPU is enough; runtime ~tens of minutes).
"""
import math, time
import numpy as np
import torch
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel

torch.manual_seed(0); rng = np.random.default_rng(0)
kB = 8.617e-5                      # eV/K
D0_TRUE, EA_TRUE = 0.76, 3.46      # cm²/s, eV
L, W, T_MAX = 600e-7, 60e-7, 3600  # 600 nm, 60 nm, 1 h
A = W / L
TEMPS_C = [925, 975, 1025, 1075]
STEPS = 4000                       # PINN Adam steps (reduce if slow)

def D_of_T(T_C):  return D0_TRUE * math.exp(-EA_TRUE / (kB * (T_C + 273.15)))
def kappa_of_D(D): return D * T_MAX / L**2
def D_of_kappa(k): return k * L**2 / T_MAX

def exact_u(xi, tau, kappa):
    s2 = A**2 + 4 * kappa * tau
    return A / np.sqrt(s2) * np.exp(-xi**2 / s2)

# ---------- 1) Synthetic SIMS data ----------
def make_sims(T_C, n=40, noise=0.05):
    xi = np.sort(rng.uniform(0, 0.8, n))                     # 0-480 nm depth
    u = exact_u(xi, 1.0, kappa_of_D(D_of_T(T_C))) * (1 + noise * rng.standard_normal(n))
    idx = rng.permutation(n)
    tr, va = idx[:30], idx[30:]                               # 30 training, 10 validation
    return (xi[tr], u[tr]), (xi[va], u[va])

data = {T: make_sims(T) for T in TEMPS_C}

# ---------- 2) Inverse PINN ----------
def grad(y, x):
    return torch.autograd.grad(y, x, torch.ones_like(y), create_graph=True)[0]

def train_pinn(train, lr, lam_data, lam_ic=1.0, width=32, steps=STEPS, seed=0):
    torch.manual_seed(seed)
    net = torch.nn.Sequential(torch.nn.Linear(2, width), torch.nn.Tanh(),
                              torch.nn.Linear(width, width), torch.nn.Tanh(),
                              torch.nn.Linear(width, width), torch.nn.Tanh(),
                              torch.nn.Linear(width, 1))
    log_k = torch.nn.Parameter(torch.tensor(math.log(0.01)))   # rough initial guess for kappa
    opt = torch.optim.Adam(list(net.parameters()) + [log_k], lr=lr)
    f = lambda x, t: net(torch.cat([x, t], 1))

    xf = torch.rand(2000, 1, requires_grad=True); tf = torch.rand(2000, 1, requires_grad=True)
    x0 = torch.rand(300, 1); t0 = torch.zeros(300, 1)
    tb = torch.rand(200, 1, requires_grad=True)
    xs = torch.zeros(200, 1, requires_grad=True); xe = torch.ones(200, 1)
    xd = torch.tensor(train[0], dtype=torch.float32)[:, None]
    ud = torch.tensor(train[1], dtype=torch.float32)[:, None]
    td = torch.ones_like(xd)

    for _ in range(steps):
        k = torch.exp(log_k)
        u = f(xf, tf)
        u_x = grad(u, xf)
        r = grad(u, tf) - k * grad(u_x, xf)
        loss = ((r**2).mean()
                + lam_ic * ((f(x0, t0) - torch.exp(-x0**2 / A**2))**2).mean()   # IC
                + (grad(f(xs, tb), xs)**2).mean()                      # u_ξ(0,τ)=0
                + (f(xe, tb)**2).mean()                                # u(1,τ)=0
                + lam_data * ((f(xd, td) - ud)**2).mean())             # SIMS data
        opt.zero_grad(); loss.backward(); opt.step()
    return net, torch.exp(log_k).item()

def val_error(kappa, val):
    """Physics-consistent validation: solve the PDE with the learned kappa and compare with held-out SIMS points.
    (Not with the network output itself: with a large lambda_data the network can fit the data while violating the physics.)"""
    pred = exact_u(val[0], 1.0, kappa)   # shortcut: analytic solver; use a finite-difference solver for general IC/BC
    return float(np.sqrt(np.mean((pred - val[1])**2)))

# ---------- 3) BO: separate inverse-PINN hyperparameter tuning at each temperature ----------
# Note: a single setting tuned at 975 C did not transfer to the other temperatures (kappa changes 40-fold),
# so BO runs separately at each temperature. The kappa of the best trial is taken as the result.
bounds = np.array([[-3.5, -2.0],    # log10(lr)
                   [0.0, 2.0],      # log10(λ_IC)
                   [0.0, 3.0]])     # log10(λ_data)
N_INIT, N_ITER = 4, 10

def ei(Xc, gp, best, xi=0.01):
    mu, s = gp.predict(Xc, return_std=True); s = np.maximum(s, 1e-9)
    z = (mu - best - xi) / s
    return (mu - best - xi) * norm.cdf(z) + s * norm.pdf(z)

def bo_tune(T):
    tr, va = data[T]
    k_true = kappa_of_D(D_of_T(T))
    Xo, Yo, Ko = [], [], []
    for it in range(N_ITER):
        if it < N_INIT:
            x = rng.uniform(bounds[:, 0], bounds[:, 1])
        else:
            gp = GaussianProcessRegressor(
                ConstantKernel() * Matern([1.0, 1.0, 1.0], length_scale_bounds=(1e-2, 1e2), nu=2.5)
                + WhiteKernel(1e-3, noise_level_bounds=(1e-8, 1e-1)),
                normalize_y=True, n_restarts_optimizer=5, random_state=0).fit(np.array(Xo), Yo)
            cand = rng.uniform(bounds[:, 0], bounds[:, 1], (3000, 3))
            x = cand[np.argmax(ei(cand, gp, max(Yo)))]
        _, k = train_pinn(tr, lr=10**x[0], lam_ic=10**x[1], lam_data=10**x[2])
        y = -math.log10(val_error(k, va))
        Xo.append(x); Yo.append(y); Ko.append(k)
        print(f"  {T}°C {it+1:2d}: lr={10**x[0]:.1e} λ_IC={10**x[1]:5.1f} λ_data={10**x[2]:6.1f} "
              f"val RMSE={10**-y:.4f}  κ={k:.4f} (true {k_true:.4f})", flush=True)
    b = int(np.argmax(Yo))
    return Ko[b], 10**-Yo[b]

print("== BO + inverse PINN: D at each temperature ==")
Ds = []
for T in TEMPS_C:
    k, rmse = bo_tune(T)
    D, Dt = D_of_kappa(k), D_of_T(T)
    Ds.append(D)
    print(f"  -> {T}°C: D_PINN={D:.3e}  D_true={Dt:.3e}  error={100*abs(D-Dt)/Dt:5.1f}%  val RMSE={rmse:.4f}\n", flush=True)

invkT = 1 / (kB * (np.array(TEMPS_C) + 273.15))
slope, intercept = np.polyfit(invkT, np.log(Ds), 1)        # ln D = ln D0 - Ea·(1/kT)
Ea, D0 = -slope, math.exp(intercept)
print(f"\n== Arrhenius ==\n  Ea = {Ea:.3f} eV (literature {EA_TRUE})   D0 = {D0:.3g} cm²/s (literature {D0_TRUE})")

# Self-check: Ea must be within 5 %
assert abs(Ea - EA_TRUE) / EA_TRUE < 0.05, "Ea bias too large: increase STEPS or widen the BO bounds"
