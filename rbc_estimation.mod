// ============================================================================
// RBC ESTIMATION EXTENSION: EURO AREA DATA
// ============================================================================
// Standalone extension of the existing RBC project.
//
// Data: euro_area_rbc_dynare.csv
// Sample: 2000Q1-2025Q4 (104 quarterly observations)
// Observables:
//   y_obs = 100 * HP-filtered log real GDP
//   n_obs = 100 * HP-filtered log total hours worked
//
// Estimated objects:
//   rho_z, rho_phi, stderr(eps_z), stderr(eps_phi)
// Structural preference/technology parameters remain calibrated.
// ============================================================================

var
    c           (long_name = 'Consumption')
    n           (long_name = 'Employment / hours')
    k           (long_name = 'End-of-period capital')
    y           (long_name = 'Output')
    z           (long_name = 'Total factor productivity')
    phi_t       (long_name = 'Labour-disutility shock')
    kn          (long_name = 'Capital per worker')
    cn          (long_name = 'Consumption per worker')
    rk          (long_name = 'Rental rate of capital')
    R           (long_name = 'Gross real interest rate')
    w           (long_name = 'Real wage')
    invest      (long_name = 'Investment')
    y_obs       (long_name = 'Observed cyclical real GDP, percent')
    n_obs       (long_name = 'Observed cyclical hours, percent')
;

varexo
    eps_z       (long_name = 'Productivity innovation')
    eps_phi     (long_name = 'Labour-disutility innovation')
;

parameters
    beta        (long_name = 'Discount factor')
    alph        (long_name = 'Capital share')
    delta       (long_name = 'Depreciation rate')
    sig         (long_name = 'Inverse intertemporal elasticity')
    frisch      (long_name = 'Inverse Frisch elasticity')
    chi         (long_name = 'Steady-state labour-disutility scale')
    n_target    (long_name = 'Target steady-state employment')
    rho_z       (long_name = 'TFP persistence')
    rho_phi     (long_name = 'Labour-disutility persistence')
;

// --------------------------------------------------------------------------
// Calibrated structural parameters: same baseline as rbc_core.mod
// --------------------------------------------------------------------------
beta     = 0.99;
alph     = 0.33;
delta    = 0.025;
sig      = 2;
frisch   = 1;
n_target = 1/3;

// Initial values for the parameters that will be estimated.
rho_z    = 0.90;
rho_phi  = 0.50;
chi      = 1;

model;
    [name = 'Gross return on capital']
    R = rk + 1 - delta;

    [name = 'Euler equation']
    c^(-sig) = beta*c(+1)^(-sig)*R(+1);

    [name = 'Intratemporal labour supply']
    w*c^(-sig) = chi*phi_t*n^frisch;

    [name = 'Consumption per worker']
    cn = c/n;

    [name = 'Production function']
    y = z*k(-1)^alph*n^(1-alph);

    [name = 'Real wage']
    w = (1-alph)*z*kn^alph;

    [name = 'Rental rate of capital']
    rk = alph*z*kn^(alph-1);

    [name = 'Capital per worker']
    kn = k(-1)/n;

    [name = 'Resource constraint']
    c + k = y + (1-delta)*k(-1);

    [name = 'Investment definition']
    invest = k - (1-delta)*k(-1);

    [name = 'TFP process']
    log(z) = rho_z*log(z(-1)) + eps_z;

    [name = 'Labour-disutility process']
    log(phi_t) = rho_phi*log(phi_t(-1)) + eps_phi;

    // Measurement equations. The empirical series are measured as percentage
    // log deviations from trend, so the model observables use the same units.
    [name = 'GDP measurement equation']
    y_obs = 100*(log(y)-log(steady_state(y)));

    [name = 'Hours measurement equation']
    n_obs = 100*(log(n)-log(steady_state(n)));
end;

steady_state_model;
    z       = 1;
    phi_t   = 1;
    R       = 1/beta;
    rk      = R - (1-delta);
    kn      = (rk/alph)^(1/(alph-1));
    w       = (1-alph)*kn^alph;
    cn      = kn^alph - delta*kn;
    chi     = w*cn^(-sig)*n_target^(-(sig+frisch));
    n       = n_target;
    k       = kn*n;
    c       = cn*n;
    y       = kn^alph*n;
    invest  = delta*k;
    y_obs   = 0;
    n_obs   = 0;
end;

// Starting shock variances. These are subsequently estimated.
shocks;
    var eps_z   = (0.01)^2;
    var eps_phi = (0.01)^2;
end;

steady;
check;
model_diagnostics;

// Exactly two observables and two structural shocks: this avoids stochastic
// singularity in the baseline estimation specification.
varobs y_obs n_obs;

// --------------------------------------------------------------------------
// Bayesian priors
// --------------------------------------------------------------------------
// Persistence parameters have Beta priors on (0,1).
// Shock standard deviations have positive inverse-Gamma priors.
estimated_params;
    rho_z,          beta_pdf,      0.90, 0.05;
    rho_phi,        beta_pdf,      0.50, 0.15;
    stderr eps_z,   inv_gamma_pdf, 0.01, 0.005;
    stderr eps_phi, inv_gamma_pdf, 0.01, 0.005;
end;

// Start the optimizer from the calibrated values above.
estimated_params_init(use_calibration);
end;

// --------------------------------------------------------------------------
// Bayesian estimation
// --------------------------------------------------------------------------
// mode_compute=6 uses Dynare's Monte-Carlo-based mode finder and tunes the
// Metropolis-Hastings proposal scale. 20,000 draws are suitable for a first
// project run; increase substantially for final posterior inference.
estimation(
    datafile='euro_area_rbc_dynare.csv',
    first_obs=1,
    nobs=104,
    mode_compute=6,
    mh_replic=20000,
    mh_nblocks=2,
    mh_drop=0.5,
    smoother,
    bayesian_irf,
    irf=40
) y_obs n_obs y n z phi_t;

// ============================================================================
// END
// ============================================================================
