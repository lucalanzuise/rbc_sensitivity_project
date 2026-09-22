// ============================================================================
// RBC FISCAL EXTENSION: EXOGENOUS GOVERNMENT-PURCHASE SHOCK
// ============================================================================
// Government purchases are financed with lump-sum taxes under a balanced
// government budget. Because taxes are non-distortionary, the aggregate
// resource constraint summarizes the private and public budget constraints.
// The tax variable is included explicitly for accounting transparency.
// ============================================================================

var
    c n k y z phi_t g tax kn cn rk R w invest
;

varexo
    eps_z eps_phi eps_g
;

parameters
    beta alph delta sig frisch chi n_target rho_z rho_phi rho_g g_share gbar
;

beta     = 0.99;
alph     = 0.33;
delta    = 0.025;
sig      = 2;
frisch   = 1;
n_target = 1/3;
rho_z    = 0.90;
rho_phi  = 0.50;
rho_g    = 0.90;
g_share  = 0.20;
gbar     = 0.10; // Updated consistently inside steady_state_model.

// Initial value. The steady_state_model block recalibrates chi analytically
// so the fiscal model shares the target employment normalization.
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

    [name = 'Aggregate resource constraint']
    c + invest + g = y;

    [name = 'Capital accumulation']
    k = (1-delta)*k(-1) + invest;

    [name = 'Balanced-budget lump-sum taxation']
    tax = g;

    [name = 'TFP process']
    log(z) = rho_z*log(z(-1)) + eps_z;

    [name = 'Labour-disutility process']
    log(phi_t) = rho_phi*log(phi_t(-1)) + eps_phi;

    [name = 'Government-purchases process']
    log(g/gbar) = rho_g*log(g(-1)/gbar) + eps_g;
end;

steady_state_model;
    z       = 1;
    phi_t   = 1;
    R       = 1/beta;
    rk      = R - (1-delta);
    kn      = (rk/alph)^(1/(alph-1));
    w       = (1-alph)*kn^alph;
    cn      = (1-g_share)*kn^alph - delta*kn;
    chi     = w*cn^(-sig)*n_target^(-(sig+frisch));
    n       = n_target;
    k       = kn*n;
    c       = cn*n;
    y       = kn^alph*n;
    invest  = delta*k;
    gbar    = g_share*y;
    g       = gbar;
    tax     = g;
end;

shocks;
    var eps_z   = (0.01)^2;
    var eps_phi = (0.01)^2;
    var eps_g   = (0.01)^2;
end;

steady;
check;
model_diagnostics;

stoch_simul(order=1, irf=40, loglinear, hp_filter=1600,
            contemporaneous_correlation, nograph, noprint)
            y c n k invest w g tax z phi_t;

fiscal = struct();
fiscal.description = 'RBC extension with exogenous government purchases';
fiscal.variables = {'y','c','n','k','invest','w','g','tax','z','phi_t'};
fiscal.shocks = {'eps_z','eps_phi','eps_g'};
fiscal.shock_std = [0.01 0.01 0.01];
fiscal.baseline.parameters = struct('beta',beta,'alph',alph,'delta',delta, ...
                                    'sig',sig,'frisch',frisch, ...
                                    'chi',get_param_by_name('chi'), ...
                                    'n_target',n_target, ...
                                    'rho_z',rho_z,'rho_phi',rho_phi, ...
                                    'rho_g',rho_g,'g_share',g_share);
fiscal.baseline.irfs = oo_.irfs;
fiscal.baseline.steady_state = oo_.steady_state;
all_endogenous_names = cellstr(M_.endo_names);
% oo_.var, oo_.autocorr and oo_.contemporaneous_correlation contain only
% the variables requested after stoch_simul, in that reported-variable order.
% Keep the full endogenous positions separately for objects such as
% oo_.steady_state, which is indexed over all endogenous variables.
fiscal.baseline.endogenous_indices = zeros(1,length(fiscal.variables));
for variable_index = 1:length(fiscal.variables)
    fiscal.baseline.endogenous_indices(variable_index) = ...
        find(strcmp(all_endogenous_names,fiscal.variables{variable_index}));
end
fiscal.baseline.variable_indices = 1:length(fiscal.variables);
fiscal.baseline.covariance = ...
    oo_.var(fiscal.baseline.variable_indices,fiscal.baseline.variable_indices);
for lag_index = 1:length(oo_.autocorr)
    fiscal.baseline.autocorrelation{lag_index} = ...
        oo_.autocorr{lag_index}(fiscal.baseline.variable_indices, ...
                                fiscal.baseline.variable_indices);
end
fiscal.baseline.correlation = ...
    oo_.contemporaneous_correlation(fiscal.baseline.variable_indices, ...
                                    fiscal.baseline.variable_indices);
fiscal.baseline.variance_decomposition = oo_.variance_decomposition;

// Persistence of government spending and the implied fiscal multipliers.
fiscal.rho_g.values = [0.00 0.50 0.90 0.98];
for experiment_index = 1:length(fiscal.rho_g.values)
    set_param_value('rho_g',fiscal.rho_g.values(experiment_index));
    steady(noprint);
    stoch_simul(order=1, irf=40, loglinear, nograph, noprint)
                y c n k invest w g tax;
    if info(1) ~= 0
        error('Dynare failed during the government-persistence exercise.');
    end
    fiscal.rho_g.irfs{experiment_index} = oo_.irfs;
    endo_names = cellstr(M_.endo_names);
    y_index = find(strcmp(endo_names,'y'));
    g_index = find(strcmp(endo_names,'g'));
    output_to_spending_ratio = oo_.steady_state(y_index)/oo_.steady_state(g_index);
    fiscal.rho_g.impact_multiplier(experiment_index) = ...
        output_to_spending_ratio*oo_.irfs.y_eps_g(1)/oo_.irfs.g_eps_g(1);
    fiscal.rho_g.cumulative_20q_multiplier(experiment_index) = ...
        output_to_spending_ratio*sum(oo_.irfs.y_eps_g(1:20))/ ...
        sum(oo_.irfs.g_eps_g(1:20));
end

set_param_value('rho_g',0.90);
steady(noprint);

if ~exist('results','dir')
    mkdir('results');
end
save(fullfile('results','fiscal_results.mat'),'fiscal');
