/**
 * holographic_kernel.h
 *
 * Holographic Entropy Memory Kernel — replacement for Cold Dark Matter
 * in the CLASS Boltzmann solver.
 *
 * Implements the power-law memory kernel from:
 *   Smith, C. G. Jr., "A Holographic Grammar of the Expanding Universe", v9.7, 2026.
 *
 * Physical mechanism:
 *   The kernel integrates baryon density perturbations δ_b(k,τ) over conformal
 *   time using a causal, power-law retarded Green's function. The result ρ_holo
 *   is a non-oscillating gravitational scaffold that replaces CDM in the Poisson
 *   equation. Von Neumann entropy weighting suppresses thermalized plasma
 *   contributions, preserving Silk damping.
 *
 * Usage in CLASS:
 *   Set Omega_cdm = 0, use_holographic_kernel = yes in your .ini file.
 *   The kernel replaces ρ_cdm · δ_cdm in the Einstein/Poisson equations.
 *
 * Targets: CLASS v3.2.x (github.com/lesgourg/class_public)
 */

#ifndef __HOLOGRAPHIC_KERNEL__
#define __HOLOGRAPHIC_KERNEL__

#include "common.h"

/* ============================================================
 * Compile-time constants
 * ============================================================ */

/** Maximum number of stored history samples per k-mode */
#define _HOLO_MAX_HISTORY_ 8192

/** Minimum number of history points before kernel integration activates */
#define _HOLO_MIN_HISTORY_ 4

/** Speed of light in m/s (for unit conversions) */
#define _HOLO_C_SI_ 2.99792458e8

/** Boltzmann constant in J/K */
#define _HOLO_KB_SI_ 1.380649e-23

/** Reduced Planck constant in J·s */
#define _HOLO_HBAR_SI_ 1.054571817e-34

/* ============================================================
 * Parameter structure
 * ============================================================ */

/**
 * Structure holding all holographic kernel parameters.
 * Stored in the perturbations workspace and initialised from input.
 */
struct holographic_parameters {

  /** Master switch */
  short use_holographic_kernel;   /**< _TRUE_ or _FALSE_ */

  /** Power-law memory kernel timescale τ_γ [seconds].
   *  Default: ~14 Gyr = 4.415e17 s. Should match 1/H_0.
   *  Derived from S_dS / ṁ (horizon entropy / entropy production rate). */
  double tau_gamma;

  /** Power-law index α of the memory kernel.
   *  Range: 1.0–1.5. From 2D CFT fluctuation-dissipation theorem.
   *  α = 1.0 gives scale-free (1/t) tail; α = 1.5 gives faster decay. */
  double kernel_alpha;

  /** Von Neumann entropy weighting switch */
  short use_vonNeumann_weighting;

  /** Phase transition decoherence spike switch */
  short include_phase_transitions;

  /** Holographic acceleration scale a_u = cH_0/(2π) [m/s²].
   *  Computed from H_0 at initialisation — not a free parameter. */
  double a_u_0;

  /** Effective holographic density parameter Ω_holo.
   *  Dynamically computed to match the CDM density fraction that
   *  the kernel replaces. Used for background evolution. */
  double Omega_holo;

  /** QCD transition temperature [eV] */
  double T_QCD;

  /** Electroweak transition temperature [eV] */
  double T_EW;

  /** Amplitude of phase-transition decoherence spikes.
   *  Dimensionless multiplicative factor on ṁ at transition. */
  double spike_amplitude;

};

/* ============================================================
 * History buffer for one k-mode
 * ============================================================ */

/**
 * Circular buffer storing the conformal-time history of δ_b
 * for a single Fourier mode k. Used to evaluate the memory
 * kernel convolution integral at each timestep.
 */
struct holo_history_buffer {

  int count;                              /**< Number of stored samples */
  int head;                               /**< Index of most recent entry */
  double tau[_HOLO_MAX_HISTORY_];         /**< Conformal time samples */
  double proper_t[_HOLO_MAX_HISTORY_];    /**< Proper time samples [s] */
  double delta_b[_HOLO_MAX_HISTORY_];     /**< Baryon density perturbation */
  double a_scale[_HOLO_MAX_HISTORY_];     /**< Scale factor at each sample */
  double T_plasma[_HOLO_MAX_HISTORY_];    /**< Plasma temperature [eV] */

};

/* ============================================================
 * Public API
 * ============================================================ */

/**
 * Initialise holographic kernel parameters from the input structure.
 *
 * @param ppr   pointer to precision structure
 * @param pba   pointer to background structure (for H_0)
 * @param pholo pointer to holographic_parameters (output)
 * @return _SUCCESS_ or _FAILURE_
 */
int holographic_kernel_init(
  struct precision * ppr,
  struct background * pba,
  struct holographic_parameters * pholo
);

/**
 * Allocate and zero a history buffer for one k-mode.
 *
 * @param buf  pointer to holo_history_buffer
 * @return _SUCCESS_ or _FAILURE_
 */
int holo_history_alloc(
  struct holo_history_buffer * buf
);

/**
 * Record a new (τ, δ_b, a, T) sample into the history buffer.
 *
 * @param buf      history buffer
 * @param tau      conformal time
 * @param proper_t proper time [seconds]
 * @param delta_b  baryon density perturbation at this (k, τ)
 * @param a        scale factor
 * @param T_eV     plasma temperature in eV
 */
void holo_history_push(
  struct holo_history_buffer * buf,
  double tau,
  double proper_t,
  double delta_b,
  double a,
  double T_eV
);

/**
 * Evaluate the retarded power-law memory kernel.
 *
 *   G_ret(Δt) = 1 / (1 + Δt/τ_γ)^α
 *
 * @param delta_t  proper time difference t - t' [seconds]
 * @param pholo    holographic parameters
 * @return kernel value in [0, 1]
 */
double holo_kernel_Gret(
  double delta_t,
  struct holographic_parameters * pholo
);

/**
 * Evaluate the von Neumann entropy weighting factor.
 *
 * Thermalized plasma (high T) is suppressed; coherent structures
 * (low T, long coherence length) contribute fully.
 *
 *   f_vN = min(1, (ξ_thermal / R_H)^p)  with smooth interpolation
 *
 * In practice, implemented as a smooth suppression function of T/T_rec
 * that transitions from ~0 (hot plasma) to ~1 (post-recombination matter).
 *
 * @param T_eV         plasma temperature [eV]
 * @param a            scale factor
 * @param H            Hubble parameter [1/s]
 * @param k            wavenumber [1/Mpc]
 * @param pholo        holographic parameters
 * @return weighting factor in [0, 1]
 */
double holo_vonNeumann_weight(
  double T_eV,
  double a,
  double H,
  double k,
  struct holographic_parameters * pholo
);

/**
 * Compute the phase-transition decoherence spike factor.
 *
 * Returns a multiplicative enhancement to the decoherence rate ṁ
 * near the QCD and electroweak transition temperatures.
 *
 * @param T_eV   plasma temperature [eV]
 * @param pholo  holographic parameters
 * @return spike factor >= 1.0
 */
double holo_phase_transition_spike(
  double T_eV,
  struct holographic_parameters * pholo
);

/**
 * Compute the holographic entropy density perturbation δ_holo(k, τ)
 * by integrating the memory kernel over the stored baryon history.
 *
 * This is the main computational routine, called once per k-mode
 * per timestep in place of δ_cdm.
 *
 *   δ_holo(k, τ) = ∫₀^τ δ_b(k, τ') · G_ret(t(τ) - t(τ')) · f_vN(τ') dτ'
 *                  × [1 + spike(T(τ'))]
 *
 * Normalised so that the time-averaged result matches the required
 * Ω_cdm equivalent gravitational density.
 *
 * @param buf          history buffer for this k-mode
 * @param tau_now      current conformal time
 * @param proper_t_now current proper time [s]
 * @param H_now        current Hubble parameter [1/s]
 * @param k            wavenumber [1/Mpc]
 * @param pholo        holographic parameters
 * @param delta_holo   output: effective holographic density perturbation
 * @return _SUCCESS_ or _FAILURE_
 */
int holo_compute_delta_holo(
  struct holo_history_buffer * buf,
  double tau_now,
  double proper_t_now,
  double H_now,
  double k,
  struct holographic_parameters * pholo,
  double * delta_holo
);

/**
 * Compute the H(z)-dependent holographic acceleration scale.
 *
 *   a_u(z) = c·H(z) / (2π)
 *
 * At high redshift, this is larger, providing stronger gravitational
 * binding — consistent with JWST early-galaxy observations.
 *
 * @param H  Hubble parameter at current redshift [1/s or 1/Mpc]
 * @return a_u in [m/s²]
 */
double holo_acceleration_scale(double H);

/**
 * Free any allocated resources.
 *
 * @param pholo  holographic parameters
 * @return _SUCCESS_ or _FAILURE_
 */
int holographic_kernel_free(
  struct holographic_parameters * pholo
);

#endif /* __HOLOGRAPHIC_KERNEL__ */
