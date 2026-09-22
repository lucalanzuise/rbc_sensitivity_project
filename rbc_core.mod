// ============================================================================
// RBC CORE MODEL: PRODUCTIVITY AND LABOUR-DISUTILITY SHOCKS
// ============================================================================
// Clean, parameterised RBC model used for the sensitivity analysis.
//
// Key design choice for comparative statics:
//   chi is recalibrated whenever sig or frisch changes so that steady-state
//   employment remains fixed at n_target. This isolates curvature/elasticity
//   effects from changes in the deterministic steady-state allocation.
// ============================================================================

var
    c           (long_name = 'Consumption')
    n           (long_name = 'Employment')
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

beta     = 0.99;
alph     = 0.33;
delta    = 0.025;
sig      = 2;
frisch   = 1;
n_target = 1/3;
rho_z    = 0.90;
rho_phi  = 0.50;

// Initial value. The steady_state_model block recalibrates chi analytically
// whenever sig or frisch changes, keeping steady-state employment fixed.
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
end;

shocks;
    var eps_z   = (0.01)^2;
    var eps_phi = (0.01)^2;
end;

steady;
check;
model_diagnostics;

// --------------------------------------------------------------------------
// 1. Baseline: IRFs, HP-filtered theoretical moments, correlations and FEVD
// --------------------------------------------------------------------------
stoch_simul(order=1, irf=40, loglinear, hp_filter=1600,
            contemporaneous_correlation, nograph, noprint)
            y c n k invest w z phi_t;

core = struct();
core.description = 'RBC model with TFP and labour-disutility shocks';
core.variables = {'y','c','n','k','invest','w','z','phi_t'};
core.shocks = {'eps_z','eps_phi'};
core.shock_std = [0.01 0.01];
core.baseline.parameters = struct('beta',beta,'alph',alph,'delta',delta, ...
                                  'sig',sig,'frisch',frisch, ...
                                  'chi',get_param_by_name('chi'), ...
                                  'n_target',n_target, ...
                                  'rho_z',rho_z,'rho_phi',rho_phi);
core.baseline.irfs = oo_.irfs;
core.baseline.steady_state = oo_.steady_state;
all_endogenous_names = cellstr(M_.endo_names);
% oo_.var, oo_.autocorr and oo_.contemporaneous_correlation contain only
% the variables requested after stoch_simul, in that reported-variable order.
% Keep the full endogenous positions separately for objects such as
% oo_.steady_state, which is indexed over all endogenous variables.
core.baseline.endogenous_indices = zeros(1,length(core.variables));
for variable_index = 1:length(core.variables)
    core.baseline.endogenous_indices(variable_index) = ...
        find(strcmp(all_endogenous_names,core.variables{variable_index}));
end
core.baseline.variable_indices = 1:length(core.variables);
core.baseline.covariance = ...
    oo_.var(core.baseline.variable_indices,core.baseline.variable_indices);
for lag_index = 1:length(oo_.autocorr)
    core.baseline.autocorrelation{lag_index} = ...
        oo_.autocorr{lag_index}(core.baseline.variable_indices, ...
                                core.baseline.variable_indices);
end
core.baseline.correlation = ...
    oo_.contemporaneous_correlation(core.baseline.variable_indices, ...
                                    core.baseline.variable_indices);
core.baseline.variance_decomposition = oo_.variance_decomposition;

// --------------------------------------------------------------------------
// 2. Preference-curvature comparison: sigma in {1, 2, 5, 10}
//    chi is recalibrated so n_ss = n_target for every sigma.
// --------------------------------------------------------------------------
core.sigma.values = [1 2 5 10];
core.sigma.chi = zeros(size(core.sigma.values));
core.sigma.steady_state_employment = zeros(size(core.sigma.values));
n_index = find(strcmp(all_endogenous_names,'n'));
for experiment_index = 1:length(core.sigma.values)
    current_sig = core.sigma.values(experiment_index);
    current_frisch = 1;
    set_param_value('sig',current_sig);
    set_param_value('frisch',current_frisch);
    steady(noprint);
    stoch_simul(order=1, irf=40, loglinear, nograph, noprint)
                y c n k invest w;
    if info(1) ~= 0
        error('Dynare failed during the sigma sensitivity exercise.');
    end
    core.sigma.chi(experiment_index) = get_param_by_name('chi');
    core.sigma.steady_state_employment(experiment_index) = ...
        oo_.steady_state(n_index);
    core.sigma.irfs{experiment_index} = oo_.irfs;
end

// --------------------------------------------------------------------------
// 3. Labour-supply comparison: inverse Frisch elasticity in {0.5, 1, 2, 5}
//    chi is recalibrated so n_ss = n_target for every frisch.
// --------------------------------------------------------------------------
core.frisch.values = [0.5 1 2 5];
core.frisch.chi = zeros(size(core.frisch.values));
core.frisch.steady_state_employment = zeros(size(core.frisch.values));
for experiment_index = 1:length(core.frisch.values)
    current_sig = 2;
    current_frisch = core.frisch.values(experiment_index);
    set_param_value('sig',current_sig);
    set_param_value('frisch',current_frisch);
    steady(noprint);
    stoch_simul(order=1, irf=40, loglinear, nograph, noprint)
                y c n k invest w;
    if info(1) ~= 0
        error('Dynare failed during the labour-supply sensitivity exercise.');
    end
    core.frisch.chi(experiment_index) = get_param_by_name('chi');
    core.frisch.steady_state_employment(experiment_index) = ...
        oo_.steady_state(n_index);
    core.frisch.irfs{experiment_index} = oo_.irfs;
end

// --------------------------------------------------------------------------
// 4. Symmetric shock-persistence comparisons
// --------------------------------------------------------------------------
set_param_value('sig',2);
set_param_value('frisch',1);

core.rho_z.values = [0.50 0.70 0.90 0.98];
for experiment_index = 1:length(core.rho_z.values)
    set_param_value('rho_z',core.rho_z.values(experiment_index));
    steady(noprint);
    stoch_simul(order=1, irf=40, loglinear, nograph, noprint)
                y c n k invest w z;
    if info(1) ~= 0
        error('Dynare failed during the TFP-persistence exercise.');
    end
    core.rho_z.irfs{experiment_index} = oo_.irfs;
end

set_param_value('rho_z',0.90);
core.rho_phi.values = [0.50 0.70 0.90 0.98];
for experiment_index = 1:length(core.rho_phi.values)
    set_param_value('rho_phi',core.rho_phi.values(experiment_index));
    steady(noprint);
    stoch_simul(order=1, irf=40, loglinear, nograph, noprint)
                y c n k invest w phi_t;
    if info(1) ~= 0
        error('Dynare failed during the labour-shock-persistence exercise.');
    end
    core.rho_phi.irfs{experiment_index} = oo_.irfs;
end

// --------------------------------------------------------------------------
// 5. Two-dimensional sigma-frisch sensitivity grid with fixed n_ss
// --------------------------------------------------------------------------
set_param_value('rho_phi',0.50);
core.grid.sigma_values = [1 2 5 10];
core.grid.frisch_values = [0.5 1 2 5];
core.grid.peak_output = zeros(length(core.grid.frisch_values), ...
                              length(core.grid.sigma_values));
core.grid.peak_employment = zeros(length(core.grid.frisch_values), ...
                                  length(core.grid.sigma_values));
core.grid.chi = zeros(length(core.grid.frisch_values), ...
                      length(core.grid.sigma_values));
core.grid.steady_state_employment = zeros(length(core.grid.frisch_values), ...
                                          length(core.grid.sigma_values));

for row_index = 1:length(core.grid.frisch_values)
    current_frisch = core.grid.frisch_values(row_index);
    for column_index = 1:length(core.grid.sigma_values)
        current_sig = core.grid.sigma_values(column_index);
        set_param_value('frisch',current_frisch);
        set_param_value('sig',current_sig);
        steady(noprint);
        stoch_simul(order=1, irf=40, loglinear, nograph, noprint) y n;
        if info(1) ~= 0
            error('Dynare failed during the sigma-frisch grid exercise.');
        end
        core.grid.chi(row_index,column_index) = get_param_by_name('chi');
        core.grid.steady_state_employment(row_index,column_index) = ...
            oo_.steady_state(n_index);
        core.grid.peak_output(row_index,column_index) = ...
            max(abs(oo_.irfs.y_eps_z));
        core.grid.peak_employment(row_index,column_index) = ...
            max(abs(oo_.irfs.n_eps_z));
    end
end

// --------------------------------------------------------------------------
// 6. Long stochastic simulation under the baseline calibration
// --------------------------------------------------------------------------
set_param_value('sig',2);
set_param_value('frisch',1);
set_param_value('rho_z',0.90);
set_param_value('rho_phi',0.50);
steady(noprint);
set_dynare_seed(20260922);
stoch_simul(order=1, periods=2200, drop=200, irf=0, loglinear,
            nograph, noprint) y c n k invest w z phi_t;

all_endogenous_names = cellstr(M_.endo_names);
simulation_indices = zeros(1,length(core.variables));
for variable_index = 1:length(core.variables)
    simulation_indices(variable_index) = ...
        find(strcmp(all_endogenous_names,core.variables{variable_index}));
end
core.simulation.variable_names = core.variables;
core.simulation.log_series = oo_.endo_simul(simulation_indices,:);
core.simulation.steady_state = oo_.steady_state(simulation_indices);
core.simulation.drop = 200;
core.simulation.seed = 20260922;

if ~exist('results','dir')
    mkdir('results');
end
save(fullfile('results','core_results.mat'),'core');
