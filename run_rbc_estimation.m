%% Euro-area RBC Bayesian estimation extension
% Runs the standalone Dynare estimation and exports compact diagnostics.

clear;
clc;
close all;

project_directory = fileparts(mfilename('fullpath'));
figures_directory = fullfile(project_directory,'figures');
tables_directory  = fullfile(project_directory,'tables');

if ~exist(figures_directory,'dir'); mkdir(figures_directory); end
if ~exist(tables_directory,'dir');  mkdir(tables_directory);  end

cd(project_directory);

if exist('dynare','file') ~= 2
    error(['Dynare is not on the MATLAB path. Start MATLAB through Dynare ', ...
           'or add the Dynare matlab folder before running this file.']);
end

required_files = {'rbc_estimation.mod','euro_area_rbc_dynare.csv'};
for file_index = 1:length(required_files)
    if exist(fullfile(project_directory,required_files{file_index}),'file') ~= 2
        error('Missing required file: %s',required_files{file_index});
    end
end

fprintf('\nRunning Bayesian estimation of the RBC model on euro-area data...\n');
dynare rbc_estimation.mod noclearall
cd(project_directory);

%% 1. Posterior parameter summary
parameter_name = {'rho_z';'rho_phi';'stderr_eps_z';'stderr_eps_phi'};
parameter_label = {'TFP persistence';'Labour-shock persistence'; ...
                   'TFP shock std. dev.';'Labour-shock std. dev.'};
prior_mean = [0.90;0.50;0.01;0.01];

posterior_mean = [ ...
    posterior_stat(oo_,'parameters','rho_z','mean'); ...
    posterior_stat(oo_,'parameters','rho_phi','mean'); ...
    posterior_stat(oo_,'shocks_std','eps_z','mean'); ...
    posterior_stat(oo_,'shocks_std','eps_phi','mean')];

hpd_lower = [ ...
    posterior_stat(oo_,'parameters','rho_z','lower'); ...
    posterior_stat(oo_,'parameters','rho_phi','lower'); ...
    posterior_stat(oo_,'shocks_std','eps_z','lower'); ...
    posterior_stat(oo_,'shocks_std','eps_phi','lower')];

hpd_upper = [ ...
    posterior_stat(oo_,'parameters','rho_z','upper'); ...
    posterior_stat(oo_,'parameters','rho_phi','upper'); ...
    posterior_stat(oo_,'shocks_std','eps_z','upper'); ...
    posterior_stat(oo_,'shocks_std','eps_phi','upper')];

posterior_table = table(parameter_name,prior_mean,posterior_mean,hpd_lower,hpd_upper, ...
    'VariableNames',{'Parameter','PriorMean','PosteriorMean','HPDLower','HPDUpper'});
writetable(posterior_table,fullfile(tables_directory,'bayesian_posterior_summary.csv'));

figure('Color','w','Position',[100 100 1050 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for i = 1:4
    nexttile;
    hold on;
    plot([hpd_lower(i) hpd_upper(i)],[1 1],'-','LineWidth',2);
    plot(posterior_mean(i),1,'o','MarkerSize',8,'LineWidth',1.5);
    plot(prior_mean(i),1,'x','MarkerSize',9,'LineWidth',1.5);
    title(parameter_label{i});
    yticks([]);
    grid on;
    legend({'Posterior HPD interval','Posterior mean','Prior mean'}, ...
           'Location','best');
end
sgtitle('Bayesian estimation: prior means and posterior estimates');
save_figure(fullfile(figures_directory,'14_bayesian_parameter_estimates.png'));

%% 2. Smoothed structural shocks
try
    tfp_shock = smoothed_shock_mean(oo_,'eps_z');
    labour_shock = smoothed_shock_mean(oo_,'eps_phi');
    T = min(length(tfp_shock),length(labour_shock));
    quarter_labels = make_quarter_labels(2000,1,T);

    figure('Color','w','Position',[100 100 1100 600]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(1:T,tfp_shock(1:T),'LineWidth',1.2);
    yline(0,'-');
    title('Smoothed productivity innovation');
    ylabel('Innovation');
    grid on;
    apply_quarter_ticks(quarter_labels);

    nexttile;
    plot(1:T,labour_shock(1:T),'LineWidth',1.2);
    yline(0,'-');
    title('Smoothed labour-disutility innovation');
    ylabel('Innovation');
    grid on;
    apply_quarter_ticks(quarter_labels);

    sgtitle('Estimated structural shocks, 2000Q1-2025Q4');
    save_figure(fullfile(figures_directory,'15_smoothed_structural_shocks.png'));
catch ME
    warning('Could not export smoothed-shock figure: %s',ME.message);
end

%% 3. Posterior-mean Bayesian IRFs
try
    irf_variables = {'y','n','c','invest'};
    irf_titles = {'Output','Hours','Consumption','Investment'};
    shock_names = {'eps_z','eps_phi'};
    shock_titles = {'Productivity shock','Labour-disutility shock'};

    for s = 1:length(shock_names)
        figure('Color','w','Position',[100 100 1050 650]);
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        for v = 1:length(irf_variables)
            response = posterior_irf_mean(oo_,irf_variables{v},shock_names{s});
            nexttile;
            plot(0:length(response)-1,100*response,'LineWidth',1.5);
            yline(0,'-');
            title(irf_titles{v});
            xlabel('Quarters after shock');
            ylabel('% deviation');
            grid on;
        end
        sgtitle(['Posterior-mean IRFs: ' shock_titles{s}]);
        save_figure(fullfile(figures_directory, ...
            sprintf('%02d_posterior_irfs_%s.png',15+s,shock_names{s})));
    end
catch ME
    warning('Could not export posterior-IRF figures: %s',ME.message);
end

fprintf('\nEstimation extension completed successfully.\n');
fprintf('Posterior table: %s\n',fullfile(tables_directory,'bayesian_posterior_summary.csv'));
fprintf('Figures: %s\n',figures_directory);

%% Local functions
function value = posterior_stat(oo,object_type,name,stat)
    switch lower(stat)
        case 'mean'
            S = oo.posterior_mean;
        case 'lower'
            S = oo.posterior_hpdinf;
        case 'upper'
            S = oo.posterior_hpdsup;
        otherwise
            error('Unknown posterior statistic: %s',stat);
    end

    candidate_objects = {object_type};
    if strcmp(object_type,'shocks_std')
        candidate_objects = {'shocks_std','shock_std','stderr'};
    end

    value = NaN;
    for j = 1:length(candidate_objects)
        obj = candidate_objects{j};
        if isfield(S,obj) && isfield(S.(obj),name)
            value = S.(obj).(name);
            return;
        end
    end
    if isfield(S,name)
        value = S.(name);
        return;
    end
    error('Posterior statistic not found for %s.',name);
end

function x = smoothed_shock_mean(oo,name)
    S = oo.SmoothedShocks;
    if isfield(S,'Mean') && isfield(S.Mean,name)
        x = S.Mean.(name);
    elseif isfield(S,'mean') && isfield(S.mean,name)
        x = S.mean.(name);
    elseif isfield(S,name)
        x = S.(name);
    else
        error('Smoothed shock %s not found.',name);
    end
    x = x(:);
end

function x = posterior_irf_mean(oo,var_name,shock_name)
    field_name = [var_name '_' shock_name];
    S = oo.PosteriorIRF.dsge;
    if isfield(S,'Mean') && isfield(S.Mean,field_name)
        x = S.Mean.(field_name);
    elseif isfield(S,'mean') && isfield(S.mean,field_name)
        x = S.mean.(field_name);
    else
        error('Posterior IRF %s not found.',field_name);
    end
    x = x(:);
end

function labels = make_quarter_labels(start_year,start_quarter,T)
    labels = cell(T,1);
    year = start_year;
    quarter = start_quarter;
    for t = 1:T
        labels{t} = sprintf('%dQ%d',year,quarter);
        quarter = quarter + 1;
        if quarter == 5
            quarter = 1;
            year = year + 1;
        end
    end
end

function apply_quarter_ticks(labels)
    T = length(labels);
    ticks = unique([1:20:T T]);
    xticks(ticks);
    xticklabels(labels(ticks));
    xtickangle(0);
end

function save_figure(output_file)
    output_directory = fileparts(output_file);
    if ~exist(output_directory,'dir'); mkdir(output_directory); end
    try
        exportgraphics(gcf,output_file,'Resolution',220);
    catch
        print(gcf,output_file,'-dpng','-r220');
    end
    close(gcf);
end
