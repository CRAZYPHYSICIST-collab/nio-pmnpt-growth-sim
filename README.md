# Surface growth simulation for NiO wedge films on PMN-PT

A 2D surface-growth simulation (linear exact solution + nonlinear pseudo-spectral solver) with a virtual AFM measurement, for NiO wedge films grown by pulsed laser deposition (PLD) on PMN-PT substrates.

**Question addressed:** can weak and strong surface diffusion (low vs. high substrate temperature) be told apart in quantities that an AFM actually measures: RMS roughness, growth exponent β and power spectral density (PSD)?

The code is a single C# file with no dependencies beyond the .NET / Mono standard library. Random seeds are fixed, so every run gives identical numbers.

## Quick start

| System | Command | Requirements |
|---|---|---|
| Windows | `powershell -ExecutionPolicy Bypass -File run.ps1 -Mode all` | PowerShell 5.1 (.NET Framework) |
| Linux / macOS | `./run.sh all` | Mono (`mcs`, `mono`) and `libgdiplus` |

Modes: `test` (a few seconds), `exact`, `virtual`, `nonlinear`, `all`. The `all` mode takes roughly 5–15 minutes and prints a summary table at the end. Results are written to `output/`.

## Model

The surface height h(x, y) evolves with the mean film thickness τ as

```
∂h/∂τ = a ∇²h − b ∇⁴h + (c/2) |∇h|² + ξ ,     ⟨ξ ξ⟩ = Ω δ(x) δ(τ)
```

| Symbol | Meaning | Unit |
|---|---|---|
| a = ν/F | surface tension coefficient / growth rate | nm |
| b = κ/F | surface diffusion coefficient / growth rate | nm³ |
| c = λ/F | slope-dependent (KPZ) term / growth rate | – |
| Ω = a₀³/4 = 0.0182 nm³ | volume of one NiO formula unit (a₀ = 0.4177 nm); strength of the shot noise | nm³ |

Because τ = F·t, all coefficients are effective coefficients scaled by the growth rate F. Surface tension and diffusion are specified through characteristic lengths at 100 nm thickness: b = ℓ_D⁴ / 100 nm, a = ℓ_ν² / 100 nm.

**Linear part (c = 0)** is solved analytically. With the discrete Laplacian eigenvalue L(k) and ω = aL + bL²:

```
S(k, τ) = Ω (1 − e^(−2ωτ)) / (2ω)  +  S₀(k) e^(−2ωτ)
```

S₀(k) is a Gaussian-correlated substrate roughness (W₀ = 0.2 nm, correlation length 20 nm). Surfaces are sampled from this PSD on a 512×512 grid (2 nm pixels, i.e. a 1×1 µm² AFM image). The AFM tip is modeled as morphological dilation with a parabolic tip of 10 nm radius (Villarrubia 1997).

**Nonlinear part (c ≠ 0)** is solved on a 256×256 grid with a pseudo-spectral method (2/3 de-aliasing) and an ETD2RK time stepper (Cox & Matthews 2002). The |∇h|² term is regularized as 1 − exp(−|∇h|²) (Dasgupta, Das Sarma & Kim 1996).

## Validation (`test` mode)

| Check | Result |
|---|---|
| FFT round-trip error | 1.8 × 10⁻¹⁵ |
| Random deposition, W(100 nm): sampled / analytic √(τΩ)/dx | 0.676 / 0.675 nm |
| Surface diffusion only (ℓ_D = 20 nm): discrete solution / continuum analytic | 0.021 / 0.021 nm |
| β in the boundary cases (random deposition / surface tension only / surface diffusion only) | 0.50 / 0.06 / 0.25 (theory 1/2, 0 (logarithmic), 1/4; Barabási & Stanley 1995) |
| Time step (strong diffusion, c = 3): W(100 nm) for dt = 0.1 / 0.05 | 0.041 / 0.040 nm |

## Summary of results

The "SUMMARY TABLE" block at the end of `output/summary_all.txt`:

| Quantity | Result |
|---|---|
| β in the boundary cases (random deposition / surface tension only / surface diffusion only) | 0.50 / 0.06 / 0.25 |
| Effective β, weak (ℓ_D = 1 nm) / strong (ℓ_D = 10 nm) diffusion, 10–120 nm, 8 repeats | 0.37 ± 0.01 / 0.250 ± 0.004 |
| Same, with 0.2 nm substrate roughness and AFM tip | 0.18 ± 0.01 / −0.106 ± 0.004 |
| Roughness ratio at 100 nm, weak / strong | 11.2 (ideal) – 2.5 (substrate + tip) |
| Spatial frequency where the PSDs separate (ratio > 3 at 120 nm) | 0.10 nm⁻¹ (tip limit ≈ 2π/R = 0.63 nm⁻¹) |
| Effect of the nonlinear term, strong diffusion (λ/F = −1 … 3) | β = 0.23–0.25 |
| Pixel-size sensitivity of β, weak diffusion (1 / 2 / 4 nm) | 0.23 / 0.36 / 0.48 |

## Output files (`output/`)

| File | Content |
|---|---|
| `summary_all.txt`, `summary_test.txt` | All numerical results and the summary table |
| `scan.csv`, `fig3_beta_scan.png`, `fig4_W100_scan.png` | β and W(100 nm) scanned over ℓ_D and ℓ_ν (exact solution) |
| `virtual_experiment_W.csv`, `fig1_W_thickness.png` | Virtual AFM measurement at 6 thicknesses along the wedge (10–120 nm), 3 regions × 8 repeats |
| `psd_120nm.csv`, `fig2_PSD.png` | Radial PSDs, measured with and without the tip |
| `fig5_surfaces.png` | Example surface maps at 100 nm |
| `nonlinear_W.csv`, `fig6_nonlinear.png` | Effect of the λ term |

## Additional validation (`validation/`)

| Folder | What it tests | Run | Result |
|---|---|---|---|
| `psd_identifiability/` | Can κ/F and its activation energy E_d be recovered from noisy AFM-like PSDs? The analytic linear model is fitted by least squares to synthetic PSDs at 4 temperatures × 6 thicknesses, 20 noise seeds. | `python psd_fit_test.py` (numpy, scipy, matplotlib; seconds) | κ/F recovered within 3 % at 400–573 K; E_d = 0.53 ± 0.06 eV (true 0.50 eV). At 300 K (diffusion length ≈ 1 nm) κ/F is not identifiable (48 % scatter) because its signature lies above the measurable k range. |
| `pinn_bo_boron/` | Generic test of the PINN + Bayesian optimization + Arrhenius chain on synthetic SIMS profiles of boron diffusion in Si (inverse PINN, hyperparameters tuned by BO at each temperature). | `python sims_pinn_bo.py` (numpy, scipy, scikit-learn, torch; tens of minutes on CPU) | E_a = 3.45 eV (true 3.46 eV); D recovered within 3–14 % at each temperature. |

## Limitations

- **Weak diffusion is not resolved by the grid.** ℓ_D = 1 nm is smaller than the pixel size (2 nm), so the weak-diffusion β = 0.37 contains a resolution effect (0.23–0.48 for 1–4 nm pixels). Read it as the effective exponent expected in a 2 nm-pixel AFM measurement, not as a theoretical exponent.
- **λ ≠ 0 with weak diffusion is not solved.** It is numerically unstable (β > 1 or overflow), so the λ scan covers strong diffusion only.
- **Absolute roughness values come from shot noise alone.** Real films have additional noise sources that increase W. Ratios and exponents are more meaningful than absolute values.

## References

- A.-L. Barabási, H. E. Stanley, *Fractal Concepts in Surface Growth*, Cambridge University Press (1995).
- S. M. Cox, P. C. Matthews, "Exponential time differencing for stiff systems", *J. Comput. Phys.* 176, 430 (2002).
- J. W. Cooley, J. W. Tukey, "An algorithm for the machine calculation of complex Fourier series", *Math. Comput.* 19, 297 (1965).
- C. Dasgupta, S. Das Sarma, J. M. Kim, "Controlled instability and multiscaling in models of epitaxial growth", *Phys. Rev. E* 54, R4552 (1996).
- G. Marsaglia, "Xorshift RNGs", *J. Stat. Softw.* 8(14), 1 (2003).
- J. S. Villarrubia, "Algorithms for scanned probe microscope image simulation, surface reconstruction, and tip estimation", *J. Res. NIST* 102, 425 (1997).
- T. D. B. Jacobs, T. Junge, L. Pastewka, "Quantitative characterization of surface topography using spectral analysis", *Surf. Topogr.: Metrol. Prop.* 5, 013001 (2017).

## Citation and license

Harun Dişli, Gebze Technical University, 2026. MIT License (see `LICENSE`).
