# CLASS Integration Guide: Holographic Entropy Memory Kernel
## Replacing CDM in the Boltzmann Solver

**Target:** CLASS v3.2.x (`github.com/lesgourg/class_public`)
**Reference:** Smith, C. G. Jr., *A Holographic Grammar of the Expanding Universe*, v9.7, 2026

---

## Overview

This guide walks through the modifications needed to replace Cold Dark Matter
in CLASS with the holographic entropy memory kernel. The patch touches five
areas of the CLASS codebase:

| File | Modification |
|------|-------------|
| `Makefile` | Add `holographic_kernel.c` to the build |
| `include/perturbations.h` | Add kernel parameters and history buffer to workspace |
| `source/input.c` | Parse new `.ini` parameters |
| `source/perturbations.c` | Replace δ_cdm with δ_holo in Poisson equation |
| `source/background.c` | Handle Ω_cdm = 0 with holographic density |

New files (copy from this package):
- `include/holographic_kernel.h`
- `source/holographic_kernel.c`

---

## Step 0: Prerequisites

```bash
git clone https://github.com/lesgourg/class_public.git
cd class_public
git checkout v3.2.2   # or latest v3.2.x
make clean
```

Copy the new files into the CLASS tree:
```bash
cp /path/to/patch/include/holographic_kernel.h include/
cp /path/to/patch/source/holographic_kernel.c source/
cp /path/to/patch/params/holographic.ini .
```

---

## Step 1: Makefile

Add `holographic_kernel.o` to the object list and compilation rules.

Find the line listing source objects (search for `SOURCE` or the `.o` list):

```makefile
# In the TOOLS or SOURCE variable, add:
SOURCE = ... perturbations.o ... holographic_kernel.o
```

Add the compilation rule (near other source rules):

```makefile
source/holographic_kernel.o: source/holographic_kernel.c include/holographic_kernel.h
	$(CC) $(OPTFLAG) $(CCFLAG) $(INCLUDES) -c $< -o $@
```

Ensure the `INCLUDES` variable contains `-Iinclude`.

---

## Step 2: include/perturbations.h

### 2a. Include the header

Near the top, after other `#include` directives:

```c
#include "holographic_kernel.h"
```

### 2b. Add parameters to `struct perturbs`

Inside `struct perturbs` (the main perturbation structure), add a block
for holographic kernel parameters. Search for members like `has_source_delta_cdm`
and add nearby:

```c
  /* ---- Holographic entropy kernel (replaces CDM) ---- */
  struct holographic_parameters holo_params;
```

### 2c. Add history buffer to the perturbation workspace

Inside `struct perturb_workspace` (usually `ppw`), add:

```c
  /* Holographic kernel: baryon history buffer for this k-mode */
  struct holo_history_buffer holo_history;
```

### 2d. Add a flag to track activation

Inside `struct perturbs`, add:

```c
  short use_holographic_kernel;
```

---

## Step 3: source/input.c

### 3a. Parse holographic parameters

In the function `input_read_parameters()` (or `input_read_parameters_general()`
depending on CLASS version), add a block after the CDM parameters are parsed:

```c
/* ---- Holographic kernel parameters ---- */
class_call(parser_read_string(pfc, "use_holographic_kernel", &string1, &flag1, errmsg),
           errmsg, errmsg);
if ((flag1 == _TRUE_) && (strstr(string1, "y") != NULL)) {
  ppt->use_holographic_kernel = _TRUE_;
} else {
  ppt->use_holographic_kernel = _FALSE_;
}

if (ppt->use_holographic_kernel == _TRUE_) {

  /* tau_gamma in Gyr, convert to seconds */
  class_read_double("tau_gamma", ppt->holo_params.tau_gamma);
  ppt->holo_params.tau_gamma *= 3.15576e16;  /* Gyr -> seconds */

  /* If not specified, will be derived from H_0 in holographic_kernel_init */
  if (ppt->holo_params.tau_gamma <= 0.0) {
    ppt->holo_params.tau_gamma = -1.0;  /* flag for auto-derivation */
  }

  class_read_double("kernel_alpha", ppt->holo_params.kernel_alpha);
  if (ppt->holo_params.kernel_alpha <= 0.0) {
    ppt->holo_params.kernel_alpha = 1.2;  /* default */
  }

  class_call(parser_read_string(pfc, "vonNeumann_weighting", &string1, &flag1, errmsg),
             errmsg, errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1, "y") != NULL)) {
    ppt->holo_params.use_vonNeumann_weighting = _TRUE_;
  } else {
    ppt->holo_params.use_vonNeumann_weighting = _FALSE_;
  }

  class_call(parser_read_string(pfc, "include_phase_transitions", &string1, &flag1, errmsg),
             errmsg, errmsg);
  if ((flag1 == _TRUE_) && (strstr(string1, "y") != NULL)) {
    ppt->holo_params.include_phase_transitions = _TRUE_;
  } else {
    ppt->holo_params.include_phase_transitions = _FALSE_;
  }

  ppt->holo_params.use_holographic_kernel = _TRUE_;
  ppt->holo_params.T_QCD = -1.0;   /* auto-set in init */
  ppt->holo_params.T_EW = -1.0;    /* auto-set in init */
  ppt->holo_params.spike_amplitude = -1.0;  /* auto-set in init */

  /* Override: if holographic kernel is on, CDM must be zero */
  if (pba->Omega0_cdm > 1e-10) {
    printf("WARNING: use_holographic_kernel=yes but omega_cdm > 0.\n");
    printf("         Setting omega_cdm = 0 (kernel replaces CDM).\n");
    pba->Omega0_cdm = 0.0;
  }

  /* Store the CDM fraction the kernel needs to reproduce */
  ppt->holo_params.Omega_holo = 0.264;  /* Planck 2018 Omega_cdm */
}
```

---

## Step 4: source/perturbations.c

This is the core modification. Three functions need changes.

### 4a. perturb_init() — Initialise the kernel

In `perturb_init()`, after the background module has been initialised,
add the kernel initialisation call:

```c
if (ppt->use_holographic_kernel == _TRUE_) {
  class_call(holographic_kernel_init(ppr, pba, &(ppt->holo_params)),
             ppt->error_message,
             ppt->error_message);
  if (ppt->perturbations_verbose > 0)
    printf(" -> holographic entropy kernel activated (replacing CDM)\n");
}
```

### 4b. perturb_workspace_init() — Allocate history buffer

In the function that initialises the per-k workspace (often called
at the beginning of each k-mode integration), add:

```c
if (ppt->use_holographic_kernel == _TRUE_) {
  class_call(holo_history_alloc(&(ppw->holo_history)),
             ppt->error_message,
             ppt->error_message);
}
```

### 4c. perturb_derivs() — Record baryon history at each timestep

Inside `perturb_derivs()` (the ODE right-hand-side function), after
δ_b is available from the current integration state, add a call to
push the current state into the history buffer:

```c
if (ppt->use_holographic_kernel == _TRUE_) {

  /* Current baryon perturbation */
  double delta_b_current = y[ppw->pv->index_pt_delta_b];

  /* Scale factor */
  double a_current = pvecback[pba->index_bg_a];

  /* Proper time: integrate dt = a·dτ
   * CLASS tracks conformal time τ; proper time t = ∫ a dτ
   * We approximate: t ≈ previous_t + a · Δτ
   * A more precise approach uses the background module's
   * proper time interpolation. */
  double tau_current = pvecback[pba->index_bg_conf_distance];
  /* Actually, use the tau passed to perturb_derivs */

  /* Plasma temperature: T(a) = T_cmb / a  in Kelvin, convert to eV */
  double T_eV = (pba->T_cmb / a_current) * 8.617333e-5;

  /* Proper time from background (if available) or approximate */
  double proper_t_current;
  if (pba->index_bg_proper_time >= 0) {
    proper_t_current = pvecback[pba->index_bg_proper_time];
  } else {
    /* Approximate: t ~ 1/(2H) in radiation era, 2/(3H) in matter era */
    double H = pvecback[pba->index_bg_H];
    proper_t_current = (a_current < 3e-4) ? 0.5/H : 0.6667/H;
  }

  holo_history_push(&(ppw->holo_history),
                    tau,  /* conformal time (function argument) */
                    proper_t_current,
                    delta_b_current,
                    a_current,
                    T_eV);
}
```

### 4d. perturb_einstein() — Replace δ_cdm with δ_holo in Poisson equation

This is the critical substitution. In `perturb_einstein()`, CLASS assembles
the source terms for the Einstein equations. The Poisson equation connects
metric perturbations to density perturbations.

**Find the CDM contribution to the Poisson equation.** In Newtonian gauge,
look for code assembling:

```c
/* Standard CLASS: something like */
ppw->delta_rho += pba->rho_cdm * delta_cdm;
/* or in synchronous gauge: */
ppw->pvecmetric[...] = ... + rho_cdm * delta_cdm ...;
```

**Replace with the holographic kernel output:**

```c
if (ppt->use_holographic_kernel == _TRUE_) {

  /* Compute δ_holo from the memory kernel convolution */
  double delta_holo = 0.0;
  double H_now = pvecback[pba->index_bg_H];
  double a_now = pvecback[pba->index_bg_a];
  double proper_t_now;
  if (pba->index_bg_proper_time >= 0) {
    proper_t_now = pvecback[pba->index_bg_proper_time];
  } else {
    proper_t_now = (a_now < 3e-4) ? 0.5/H_now : 0.6667/H_now;
  }

  holo_compute_delta_holo(
    &(ppw->holo_history),
    tau,
    proper_t_now,
    H_now,
    k,
    &(ppt->holo_params),
    &delta_holo
  );

  /* Effective holographic density: use Omega_holo × ρ_crit as the density */
  double rho_holo = ppt->holo_params.Omega_holo * pvecback[pba->index_bg_rho_crit];

  /* Add to the Poisson equation source (where ρ_cdm·δ_cdm would go) */
  ppw->delta_rho += rho_holo * delta_holo;

  /* The velocity perturbation θ_holo is negligible for the kernel
   * (the scaffold is non-oscillating by construction), so:
   * θ_holo ≈ 0  — no contribution to the momentum constraint. */

} else {

  /* Standard CDM contribution (unchanged) */
  ppw->delta_rho += pvecback[pba->index_bg_rho_cdm] * delta_cdm;
  /* ... original CDM velocity terms ... */

}
```

**IMPORTANT:** In synchronous gauge (CLASS default), the modifications are
analogous but involve the synchronous metric variables `h` and `η` rather
than `ψ`. The key point is the same: replace the CDM density perturbation
source with the kernel output. Search for `rho_cdm` or `delta_cdm` in
`perturb_einstein()` and substitute.

### 4e. perturb_free() — Clean up

In `perturb_free()`, add:

```c
if (ppt->use_holographic_kernel == _TRUE_) {
  holographic_kernel_free(&(ppt->holo_params));
}
```

---

## Step 5: source/background.c

### 5a. Handle Ω_cdm = 0 gracefully

CLASS may divide by ρ_cdm or check Ω_cdm > 0 in various places. With the
holographic kernel, Ω_cdm = 0 but the total matter budget must still close.

In `background_solve()`, after computing the density fractions, add:

```c
/* If holographic kernel is active, add Omega_holo to the total
 * matter budget for the Friedmann equation. The kernel provides
 * gravitational scaffolding equivalent to CDM at the background level. */
if (ppt->use_holographic_kernel == _TRUE_) {
  /* The background evolution uses an effective CDM-like component
   * that tracks the kernel's time-averaged density.
   * For the Friedmann equation: H^2 = (8πG/3)(ρ_b + ρ_holo + ρ_γ + ρ_ν + ρ_Λ)
   *
   * At the background level, ρ_holo scales as a^{-3} (like CDM),
   * because the kernel's time-averaged effect on the expansion rate
   * is indistinguishable from pressureless matter.
   *
   * The PERTURBATION-level difference is where the physics lives:
   * δ_holo is non-oscillating (from the kernel), unlike a fluid δ_cdm. */
  pba->Omega0_cdm_eff = ppt->holo_params.Omega_holo;
}
```

**NOTE:** The cleanest approach may be to keep `omega_cdm` at its standard
value (0.12) for background evolution but disable the CDM *perturbation*
equations. This avoids modifying the Friedmann equation while still
replacing CDM's perturbative role. In this case:

- Background: standard ΛCDM with Ω_cdm = 0.264 (provides correct H(z))
- Perturbations: δ_cdm and θ_cdm are NOT evolved; δ_holo replaces them

To implement this cleaner approach, in `input.c`:
```c
/* Keep background CDM for Friedmann equation */
/* pba->Omega0_cdm stays at standard value */

/* But flag perturbation module to skip CDM fluid equations
 * and use holographic kernel instead */
ppt->use_holographic_kernel = _TRUE_;
ppt->skip_cdm_perturbations = _TRUE_;
```

And in `perturb_derivs()`, wrap the CDM evolution equations:
```c
if (ppt->skip_cdm_perturbations == _FALSE_) {
  /* Standard CDM perturbation evolution: δ'_cdm, θ'_cdm */
  dy[ppw->pv->index_pt_delta_cdm] = ...;
  dy[ppw->pv->index_pt_theta_cdm] = ...;
} else {
  /* Holographic kernel: CDM perturbations are not evolved */
  dy[ppw->pv->index_pt_delta_cdm] = 0.0;
  dy[ppw->pv->index_pt_theta_cdm] = 0.0;
}
```

---

## Step 6: Build and Test

```bash
make clean
make

# Run with holographic kernel
./class holographic.ini

# Run standard ΛCDM comparison
./class explanatory.ini

# Compare C_ℓ outputs
python scripts/compare_cl.py output/holographic_cl.dat output/lcdm_cl.dat
```

---

## Step 7: Validation Protocol (Appendix A)

Run the comparison script (provided in `scripts/validate_holo.py`) to check:

| Test | Target | Pass criterion |
|------|--------|---------------|
| Peak 1 height | Match Planck | Within 5% |
| Peak 1/Peak 2 ratio | ~2.5–3× | Within 20% of ΛCDM |
| Peak 1/Peak 3 ratio | ~2.0–2.5× | Within 20% of ΛCDM |
| Peak positions (ℓ) | Match Planck | Within 1% |
| Damping tail (ℓ > 1500) | Silk damping preserved | No over-deepening |
| Low-ℓ (ℓ < 10) | Slight power deficit | Qualitative check |

**If the peak ratios are qualitatively correct but quantitatively off:**
- Adjust `kernel_alpha` in the range [1.0, 1.5]
- Check `tau_gamma`: should be ~1/H_0. Try 10–20 Gyr range.
- Verify von Neumann weighting is active (essential for Silk damping).

**If the mechanism completely fails:**
- Check that δ_holo is actually non-zero (add debug prints)
- Verify history buffer is being populated (check `holo_history.count`)
- Ensure the kernel integral is not dominated by very early times
  (the von Neumann suppression should prevent this)

---

## Architecture Summary

```
                         ┌─────────────────────┐
                         │   CLASS Background   │
                         │   (standard ΛCDM     │
                         │   with Ω_cdm for H(z))│
                         └──────────┬──────────┘
                                    │ H(z), a(τ), T(τ)
                                    ▼
┌────────────────────────────────────────────────────────┐
│                CLASS Perturbations                      │
│                                                         │
│  For each k-mode:                                       │
│                                                         │
│  1. Evolve δ_b, θ_b, δ_γ, θ_γ, δ_ν, θ_ν  (standard)  │
│                                                         │
│  2. At each timestep:                                   │
│     holo_history_push(τ, t, δ_b, a, T)                 │
│                                                         │
│  3. In Einstein equations:                              │
│     ┌─────────────────────────────────────────┐         │
│     │  δ_holo = ∫ δ_b(τ') · G_ret(t-t')      │         │
│     │            · f_vN(T') · spike(T') dτ'   │         │
│     │                                          │         │
│     │  Replace: ρ_cdm·δ_cdm → ρ_holo·δ_holo  │         │
│     └─────────────────────────────────────────┘         │
│                                                         │
│  4. Assemble C_ℓ (standard transfer → spectra)          │
└────────────────────────────────────────────────────────┘
```

---

## Troubleshooting

**Q: CLASS crashes with segfault**
A: Check that `holo_history_alloc()` is called before `holo_history_push()`.
   The buffer must be initialised for each k-mode.

**Q: C_ℓ output is zero or NaN**
A: Add debug prints in `holo_compute_delta_holo()` to verify the integral
   is finite. Common cause: `proper_t` not correctly computed.

**Q: Peak asymmetry is too weak**
A: Increase `kernel_alpha` (try 1.3–1.5) or check that the normalisation
   `N_amplitude = Omega_cdm / Omega_b` is correctly applied.

**Q: Peak asymmetry is too strong / Silk damping destroyed**
A: Verify `vonNeumann_weighting = yes`. The thermal suppression is essential
   to prevent over-deepening of gravitational wells pre-recombination.

**Q: Build errors about undefined `_SUCCESS_`, `_TRUE_`**
A: Ensure `#include "common.h"` is present in `holographic_kernel.h` and
   that the include path is set correctly in the Makefile.

---

## Performance Notes

The memory kernel convolution is O(N_history) per k-mode per timestep.
With default settings (_HOLO_MAX_HISTORY_ = 8192), this adds ~30-50%
to CLASS runtime compared to standard ΛCDM. To reduce:
- Decrease `_HOLO_MAX_HISTORY_` to 4096 (loses some early-universe memory)
- Subsample the history (push every 2nd or 4th timestep)
- Pre-compute the kernel table and interpolate

For production runs, the kernel values `G_ret(Δt)` can be tabulated at
initialisation and looked up via binary search, reducing the cost to
O(N_history) multiplications + one lookup per sample.

---

## Citation

If this patch produces results used in a publication, please cite:

```bibtex
@article{Smith2026,
  author  = {Smith, C. G., Jr.},
  title   = {A Holographic Grammar of the Expanding Universe},
  year    = {2026},
  version = {9.7},
  note    = {Section 5, Appendix A}
}
```
