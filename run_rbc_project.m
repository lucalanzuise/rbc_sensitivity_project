%% RBC DSGE extension project
% Run this file from MATLAB after Dynare has been added to the MATLAB path.
% It solves both models and creates all figures and tables used in the project.

clear;
clc;
close all;

project_directory = fileparts(mfilename('fullpath'));
cd(project_directory);

if exist('dynare','file') ~= 2
    error(['Dynare is not on the MATLAB path. Start MATLAB through Dynare ', ...
           'or add the Dynare matlab folder before running this file.']);
end

results_directory = fullfile(project_directory,'results');
figures_directory = fullfile(project_directory,'figures');
tables_directory = fullfile(project_directory,'tables');

output_directories = {results_directory,figures_directory,tables_directory};
for directory_index = 1:length(output_directories)
    if ~exist(output_directories{directory_index},'dir')
        mkdir(output_directories{directory_index});
    end
end

fprintf('\n1/3 Solving the core RBC model and running sensitivity exercises...\n');
dynare rbc_core.mod noclearall nolog
cd(project_directory);

fprintf('\n2/3 Solving the fiscal extension...\n');
dynare rbc_fiscal.mod noclearall nolog
cd(project_directory);

fprintf('\n3/3 Creating figures and tables...\n');
load(fullfile(results_directory,'core_results.mat'),'core');
load(fullfile(results_directory,'fiscal_results.mat'),'fiscal');

plot_variables = {'y','c','n','k','invest','w'};
plot_titles = {'Output','Consumption','Employment','Capital','Investment','Real wage'};

%% 1. Baseline impulse responses
plot_irf_panel(core.baseline.irfs,'eps_z',plot_variables,plot_titles, ...
    'Baseline response to a one-standard-deviation productivity shock', ...
    fullfile(figures_directory,'01_baseline_tfp_irfs.png'));

plot_irf_panel(core.baseline.irfs,'eps_phi',plot_variables,plot_titles, ...
    'Baseline response to a one-standard-deviation labour-disutility shock', ...
    fullfile(figures_directory,'02_baseline_labour_irfs.png'));

%% 2. Parameter sensitivity
plot_sensitivity_panel(core.sigma.irfs,core.sigma.values,'eps_z', ...
    plot_variables,plot_titles,'\sigma', ...
    'TFP responses under alternative inverse intertemporal elasticities', ...
    fullfile(figures_directory,'03_sigma_sensitivity.png'));

plot_sensitivity_panel(core.frisch.irfs,core.frisch.values,'eps_z', ...
    plot_variables,plot_titles,'\varphi', ...
    'TFP responses under alternative inverse Frisch elasticities', ...
    fullfile(figures_directory,'04_labour_elasticity_sensitivity.png'));

plot_sensitivity_panel(core.rho_z.irfs,core.rho_z.values,'eps_z', ...
    plot_variables,plot_titles,'\rho_z', ...
    'TFP responses under alternative productivity persistence', ...
    fullfile(figures_directory,'05_tfp_persistence_sensitivity.png'));

plot_sensitivity_panel(core.rho_phi.irfs,core.rho_phi.values,'eps_phi', ...
    plot_variables,plot_titles,'\rho_\phi', ...
    'Labour-shock responses under alternative shock persistence', ...
    fullfile(figures_directory,'06_labour_persistence_sensitivity.png'));

%% 3. Direct shock comparison
plot_shock_comparison_panel(core.baseline.irfs, ...
    'eps_z','Productivity shock','eps_phi','Labour-disutility shock', ...
    {'y','c','n','invest','w','k'}, ...
    {'Output','Consumption','Employment','Investment','Real wage','Capital'}, ...
    'Baseline transmission: productivity versus labour-disutility shocks', ...
    fullfile(figures_directory,'07_tfp_vs_labour_shock.png'));

%% 4. Quantitative response metrics
baseline_metrics = build_metrics_table(core.baseline.irfs, ...
    plot_variables,{'eps_z','eps_phi'},40);
writetable(baseline_metrics,fullfile(tables_directory,'baseline_irf_metrics.csv'));

sigma_metrics = build_sensitivity_metrics_table(core.sigma.irfs, ...
    core.sigma.values,'sigma','eps_z',plot_variables,40);
writetable(sigma_metrics,fullfile(tables_directory,'sigma_sensitivity_metrics.csv'));

frisch_metrics = build_sensitivity_metrics_table(core.frisch.irfs, ...
    core.frisch.values,'varphi','eps_z',plot_variables,40);
writetable(frisch_metrics,fullfile(tables_directory,'frisch_sensitivity_metrics.csv'));

rho_z_metrics = build_sensitivity_metrics_table(core.rho_z.irfs, ...
    core.rho_z.values,'rho_z','eps_z',plot_variables,40);
writetable(rho_z_metrics,fullfile(tables_directory,'rho_z_sensitivity_metrics.csv'));

rho_phi_metrics = build_sensitivity_metrics_table(core.rho_phi.irfs, ...
    core.rho_phi.values,'rho_phi','eps_phi',plot_variables,40);
writetable(rho_phi_metrics,fullfile(tables_directory,'rho_phi_sensitivity_metrics.csv'));

%% 5. Steady-state control check for preference experiments
number_sigma = length(core.sigma.values);
number_frisch = length(core.frisch.values);
experiment_type = [repmat({'sigma'},number_sigma,1); ...
                   repmat({'varphi'},number_frisch,1)];
parameter_value = [core.sigma.values(:); core.frisch.values(:)];
steady_state_employment = [core.sigma.steady_state_employment(:); ...
                           core.frisch.steady_state_employment(:)];
labour_disutility_scale = [core.sigma.chi(:); core.frisch.chi(:)];
target_employment = repmat(core.baseline.parameters.n_target, ...
                           number_sigma+number_frisch,1);
absolute_error = abs(steady_state_employment-target_employment);
if max(absolute_error) > 1e-10
    error('Steady-state employment was not held fixed in a preference experiment.');
end
steady_state_control_table = table(experiment_type,parameter_value, ...
    labour_disutility_scale,steady_state_employment,target_employment, ...
    absolute_error, ...
    'VariableNames',{'Experiment','ParameterValue','LabourDisutilityScaleChi', ...
                     'SteadyStateEmployment','TargetEmployment','AbsoluteError'});
writetable(steady_state_control_table, ...
    fullfile(tables_directory,'steady_state_employment_control.csv'));

%% 6. Sigma-varphi sensitivity heatmaps
figure('Color','w','Position',[100 100 1100 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(core.grid.sigma_values,core.grid.frisch_values, ...
        100*core.grid.peak_output);
set(gca,'YDir','normal');
colorbar;
xlabel('Inverse intertemporal elasticity, \sigma');
ylabel('Inverse Frisch elasticity, \varphi');
title('Peak absolute output response (%)');

nexttile;
imagesc(core.grid.sigma_values,core.grid.frisch_values, ...
        100*core.grid.peak_employment);
set(gca,'YDir','normal');
colorbar;
xlabel('Inverse intertemporal elasticity, \sigma');
ylabel('Inverse Frisch elasticity, \varphi');
title('Peak absolute employment response (%)');

sgtitle('Sensitivity to household preference parameters (fixed steady-state employment)');
save_figure(fullfile(figures_directory,'08_sigma_varphi_heatmaps.png'));

%% 7. HP-filtered theoretical moments, correlations and autocorrelations
standard_deviation = 100*sqrt(diag(core.baseline.covariance));
hp_moments_table = table(core.variables(:),standard_deviation, ...
    'VariableNames',{'Variable','HPFilteredStdDevPercent'});

maximum_autocorrelation_lag = min(5,length(core.baseline.autocorrelation));
autocorrelation_table = table(core.variables(:),'VariableNames',{'Variable'});
for lag_index = 1:maximum_autocorrelation_lag
    lag_values = diag(core.baseline.autocorrelation{lag_index});
    variable_name = sprintf('Autocorrelation%d',lag_index);
    hp_moments_table.(variable_name) = lag_values;
    autocorrelation_table.(variable_name) = lag_values;
end
writetable(hp_moments_table, ...
    fullfile(tables_directory,'baseline_hp_filtered_moments.csv'));
writetable(autocorrelation_table, ...
    fullfile(tables_directory,'baseline_autocorrelations.csv'));

correlation_table = array2table(core.baseline.correlation, ...
    'VariableNames',core.variables,'RowNames',core.variables);
writetable(correlation_table,fullfile(tables_directory,'baseline_correlations.csv'), ...
    'WriteRowNames',true);

variance_decomposition_table = array2table( ...
    core.baseline.variance_decomposition, ...
    'VariableNames',{'ProductivityShockPercent','LabourShockPercent'}, ...
    'RowNames',core.variables);
writetable(variance_decomposition_table, ...
    fullfile(tables_directory,'baseline_variance_decomposition.csv'), ...
    'WriteRowNames',true);

plot_variance_decomposition(core.baseline.variance_decomposition, ...
    core.variables,{'y','c','n','invest'}, ...
    {'Output','Consumption','Employment','Investment'}, ...
    {'Productivity shock','Labour-disutility shock'}, ...
    'Baseline unconditional variance decomposition', ...
    fullfile(figures_directory,'09_baseline_variance_decomposition.png'));

%% 8. A reproducible stochastic sample
burn_in = core.simulation.drop;
simulation = core.simulation.log_series;
if size(simulation,2) > burn_in
    simulation = simulation(:,burn_in+1:end);
end
simulation_deviations = 100*(simulation-log(core.simulation.steady_state(:)));
sample_length = min(200,size(simulation_deviations,2));

figure('Color','w','Position',[100 100 1100 680]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for variable_index = 1:length(plot_variables)
    stored_index = find(strcmp(core.simulation.variable_names, ...
                               plot_variables{variable_index}));
    nexttile;
    plot(1:sample_length, ...
         simulation_deviations(stored_index,1:sample_length), ...
         'LineWidth',1.1);
    yline(0,'Color',[0.45 0.45 0.45]);
    title(plot_titles{variable_index});
    xlabel('Quarter');
    ylabel('% deviation');
    grid on;
end
sgtitle('Illustrative stochastic simulation under the baseline calibration');
save_figure(fullfile(figures_directory,'10_stochastic_simulation.png'));

%% 9. Government-purchases extension
plot_irf_panel(fiscal.baseline.irfs,'eps_g', ...
    {'y','c','n','invest','w','g'}, ...
    {'Output','Consumption','Employment','Investment','Real wage', ...
     'Government purchases'}, ...
    'Response to a one-standard-deviation government-purchases shock', ...
    fullfile(figures_directory,'11_government_spending_irfs.png'));

plot_sensitivity_panel(fiscal.rho_g.irfs,fiscal.rho_g.values,'eps_g', ...
    {'y','c','n','invest','w','g'}, ...
    {'Output','Consumption','Employment','Investment','Real wage', ...
     'Government purchases'}, ...
    '\rho_g','Fiscal responses under alternative spending persistence', ...
    fullfile(figures_directory,'12_fiscal_persistence_sensitivity.png'));

fiscal_rho_g_metrics = build_sensitivity_metrics_table(fiscal.rho_g.irfs, ...
    fiscal.rho_g.values,'rho_g','eps_g', ...
    {'y','c','n','invest','w','g'},40);
writetable(fiscal_rho_g_metrics, ...
    fullfile(tables_directory,'fiscal_rho_g_sensitivity_metrics.csv'));

fiscal_multiplier_table = table(fiscal.rho_g.values(:), ...
    fiscal.rho_g.impact_multiplier(:), ...
    fiscal.rho_g.cumulative_20q_multiplier(:), ...
    'VariableNames',{'GovernmentPersistence','ImpactMultiplier', ...
                     'CumulativeMultiplier20Q'});
writetable(fiscal_multiplier_table, ...
    fullfile(tables_directory,'fiscal_multipliers.csv'));

fiscal_variance_decomposition_table = array2table( ...
    fiscal.baseline.variance_decomposition, ...
    'VariableNames',{'ProductivityShockPercent','LabourShockPercent', ...
                     'GovernmentShockPercent'}, ...
    'RowNames',fiscal.variables);
writetable(fiscal_variance_decomposition_table, ...
    fullfile(tables_directory,'fiscal_variance_decomposition.csv'), ...
    'WriteRowNames',true);

plot_variance_decomposition(fiscal.baseline.variance_decomposition, ...
    fiscal.variables,{'y','c','n','invest'}, ...
    {'Output','Consumption','Employment','Investment'}, ...
    {'Productivity shock','Labour-disutility shock','Government shock'}, ...
    'Fiscal model unconditional variance decomposition', ...
    fullfile(figures_directory,'13_fiscal_variance_decomposition.png'));

%% 10. Research-question map
research_question = (1:5)';
question = { ...
    'How does the intertemporal elasticity of substitution affect consumption smoothing and capital accumulation?'; ...
    'How does the Frisch elasticity affect the labour-market response to shocks?'; ...
    'How does shock persistence change the magnitude and duration of business-cycle fluctuations?'; ...
    'How do productivity and labour-preference shocks generate different macroeconomic dynamics?'; ...
    'Which shocks account for fluctuations in output, consumption, employment and investment?'};
main_statistic = { ...
    'Impact/peak consumption plus cumulative investment and capital responses'; ...
    'Impact/peak employment, output and wage responses'; ...
    'Impact, signed peak, time to peak, half-life, cumulative 20Q and cumulative 40Q responses'; ...
    'Side-by-side IRFs to equally sized one-standard-deviation innovations'; ...
    'Unconditional variance shares under the calibrated shock variances'};
primary_figure = { ...
    '03_sigma_sensitivity.png'; ...
    '04_labour_elasticity_sensitivity.png'; ...
    '05_tfp_persistence_sensitivity.png and 06_labour_persistence_sensitivity.png'; ...
    '07_tfp_vs_labour_shock.png'; ...
    '09_baseline_variance_decomposition.png'};
primary_table = { ...
    'sigma_sensitivity_metrics.csv'; ...
    'frisch_sensitivity_metrics.csv'; ...
    'rho_z_sensitivity_metrics.csv and rho_phi_sensitivity_metrics.csv'; ...
    'baseline_irf_metrics.csv'; ...
    'baseline_variance_decomposition.csv'};
research_question_map = table(research_question,question,main_statistic, ...
    primary_figure,primary_table, ...
    'VariableNames',{'ResearchQuestion','Question','MainStatistic', ...
                     'PrimaryFigure','PrimaryTable'});
writetable(research_question_map, ...
    fullfile(tables_directory,'research_question_map.csv'));

fprintf('\nProject completed successfully.\n');
fprintf('Figures: %s\n',figures_directory);
fprintf('Tables:  %s\n',tables_directory);
fprintf(['Note: variance decompositions describe the calibrated model under ', ...
         'the assumed shock variances; they are not empirical estimates.\n']);

%% Local functions
function plot_irf_panel(irfs,shock_name,variable_names,variable_titles, ...
                        overall_title,output_file)
    figure('Color','w','Position',[100 100 1100 680]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    for variable_index = 1:length(variable_names)
        field_name = [variable_names{variable_index} '_' shock_name];
        response = 100*irfs.(field_name);
        quarters = 0:length(response)-1;
        nexttile;
        plot(quarters,response,'LineWidth',1.6);
        yline(0,'Color',[0.45 0.45 0.45]);
        title(variable_titles{variable_index});
        xlabel('Quarters after shock');
        ylabel('% deviation');
        grid on;
    end
    sgtitle(overall_title);
    save_figure(output_file);
end

function plot_sensitivity_panel(irf_collection,parameter_values,shock_name, ...
                                variable_names,variable_titles,parameter_label, ...
                                overall_title,output_file)
    colours = lines(length(parameter_values));
    legend_labels = arrayfun(@(value) sprintf('%s = %g',parameter_label,value), ...
                             parameter_values,'UniformOutput',false);
    figure('Color','w','Position',[100 100 1100 680]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    for variable_index = 1:length(variable_names)
        nexttile;
        hold on;
        for experiment_index = 1:length(parameter_values)
            field_name = [variable_names{variable_index} '_' shock_name];
            response = 100*irf_collection{experiment_index}.(field_name);
            quarters = 0:length(response)-1;
            plot(quarters,response,'LineWidth',1.45, ...
                 'Color',colours(experiment_index,:));
        end
        yline(0,'Color',[0.45 0.45 0.45]);
        title(variable_titles{variable_index});
        xlabel('Quarters after shock');
        ylabel('% deviation');
        grid on;
        if variable_index == 1
            legend(legend_labels,'Location','best');
        end
    end
    sgtitle(overall_title);
    save_figure(output_file);
end

function plot_shock_comparison_panel(irfs,shock_1,shock_1_label, ...
                                     shock_2,shock_2_label, ...
                                     variable_names,variable_titles, ...
                                     overall_title,output_file)
    figure('Color','w','Position',[100 100 1100 680]);
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    for variable_index = 1:length(variable_names)
        field_1 = [variable_names{variable_index} '_' shock_1];
        field_2 = [variable_names{variable_index} '_' shock_2];
        response_1 = 100*irfs.(field_1);
        response_2 = 100*irfs.(field_2);
        quarters_1 = 0:length(response_1)-1;
        quarters_2 = 0:length(response_2)-1;
        nexttile;
        hold on;
        plot(quarters_1,response_1,'LineWidth',1.6);
        plot(quarters_2,response_2,'--','LineWidth',1.6);
        yline(0,'Color',[0.45 0.45 0.45]);
        title(variable_titles{variable_index});
        xlabel('Quarters after shock');
        ylabel('% deviation');
        grid on;
        if variable_index == 1
            legend({shock_1_label,shock_2_label},'Location','best');
        end
    end
    sgtitle(overall_title);
    save_figure(output_file);
end

function plot_variance_decomposition(decomposition,all_variable_names, ...
                                     selected_variables,selected_titles, ...
                                     shock_titles,overall_title,output_file)
    row_indices = zeros(1,length(selected_variables));
    for variable_index = 1:length(selected_variables)
        row_indices(variable_index) = ...
            find(strcmp(all_variable_names,selected_variables{variable_index}));
    end
    selected_decomposition = decomposition(row_indices,:);
    if max(sum(selected_decomposition,2)) <= 1.01
        selected_decomposition = 100*selected_decomposition;
    end

    figure('Color','w','Position',[100 100 900 520]);
    bar(selected_decomposition,'stacked');
    xticks(1:length(selected_titles));
    xticklabels(selected_titles);
    ylabel('% of unconditional variance');
    ylim([0 100]);
    legend(shock_titles,'Location','eastoutside');
    title(overall_title);
    grid on;
    save_figure(output_file);
end

function metrics_table = build_metrics_table(irfs,variable_names,shock_names,horizon)
    number_of_rows = length(variable_names)*length(shock_names);
    variable_column = cell(number_of_rows,1);
    shock_column = cell(number_of_rows,1);
    impact_column = zeros(number_of_rows,1);
    peak_column = zeros(number_of_rows,1);
    peak_quarter_column = zeros(number_of_rows,1);
    half_life_column = nan(number_of_rows,1);
    cumulative_20q_column = zeros(number_of_rows,1);
    cumulative_40q_column = zeros(number_of_rows,1);
    row_index = 0;

    for shock_index = 1:length(shock_names)
        for variable_index = 1:length(variable_names)
            row_index = row_index + 1;
            field_name = [variable_names{variable_index} '_' ...
                          shock_names{shock_index}];
            response = 100*irfs.(field_name);
            response = response(1:min(horizon,length(response)));
            [impact_response,peak_response,peak_quarter,half_life, ...
             cumulative_20q,cumulative_40q] = response_metrics(response);

            variable_column{row_index} = variable_names{variable_index};
            shock_column{row_index} = shock_names{shock_index};
            impact_column(row_index) = impact_response;
            peak_column(row_index) = peak_response;
            peak_quarter_column(row_index) = peak_quarter;
            half_life_column(row_index) = half_life;
            cumulative_20q_column(row_index) = cumulative_20q;
            cumulative_40q_column(row_index) = cumulative_40q;
        end
    end

    metrics_table = table(variable_column,shock_column,impact_column, ...
        peak_column,peak_quarter_column,half_life_column, ...
        cumulative_20q_column,cumulative_40q_column, ...
        'VariableNames',{'Variable','Shock','ImpactPercent','PeakPercent', ...
                         'PeakQuarter','HalfLifeAfterPeak', ...
                         'CumulativeResponse20Q','CumulativeResponse40Q'});
end

function metrics_table = build_sensitivity_metrics_table(irf_collection, ...
    parameter_values,parameter_name,shock_name,variable_names,horizon)

    number_of_rows = length(parameter_values)*length(variable_names);
    parameter_column = repmat({parameter_name},number_of_rows,1);
    value_column = zeros(number_of_rows,1);
    variable_column = cell(number_of_rows,1);
    shock_column = repmat({shock_name},number_of_rows,1);
    impact_column = zeros(number_of_rows,1);
    peak_column = zeros(number_of_rows,1);
    peak_quarter_column = zeros(number_of_rows,1);
    half_life_column = nan(number_of_rows,1);
    cumulative_20q_column = zeros(number_of_rows,1);
    cumulative_40q_column = zeros(number_of_rows,1);
    row_index = 0;

    for experiment_index = 1:length(parameter_values)
        for variable_index = 1:length(variable_names)
            row_index = row_index + 1;
            field_name = [variable_names{variable_index} '_' shock_name];
            response = 100*irf_collection{experiment_index}.(field_name);
            response = response(1:min(horizon,length(response)));
            [impact_response,peak_response,peak_quarter,half_life, ...
             cumulative_20q,cumulative_40q] = response_metrics(response);

            value_column(row_index) = parameter_values(experiment_index);
            variable_column{row_index} = variable_names{variable_index};
            impact_column(row_index) = impact_response;
            peak_column(row_index) = peak_response;
            peak_quarter_column(row_index) = peak_quarter;
            half_life_column(row_index) = half_life;
            cumulative_20q_column(row_index) = cumulative_20q;
            cumulative_40q_column(row_index) = cumulative_40q;
        end
    end

    metrics_table = table(parameter_column,value_column,variable_column, ...
        shock_column,impact_column,peak_column,peak_quarter_column, ...
        half_life_column,cumulative_20q_column,cumulative_40q_column, ...
        'VariableNames',{'Parameter','ParameterValue','Variable','Shock', ...
                         'ImpactPercent','PeakPercent','PeakQuarter', ...
                         'HalfLifeAfterPeak','CumulativeResponse20Q', ...
                         'CumulativeResponse40Q'});
end

function [impact_response,peak_response,peak_quarter,half_life, ...
          cumulative_20q,cumulative_40q] = response_metrics(response)
    impact_response = response(1);
    [~,peak_index] = max(abs(response));
    peak_response = response(peak_index);
    peak_quarter = peak_index-1;

    if abs(peak_response) < eps
        half_life = 0;
    else
        half_threshold = 0.5*abs(peak_response);
        half_index = find(abs(response(peak_index:end)) <= half_threshold,1);
        if isempty(half_index)
            half_life = NaN;
        else
            half_life = half_index-1;
        end
    end

    cumulative_20q = sum(response(1:min(20,length(response))));
    cumulative_40q = sum(response(1:min(40,length(response))));
end

function save_figure(output_file)
    output_directory = fileparts(output_file);
    if ~isempty(output_directory) && ~exist(output_directory,'dir')
        mkdir(output_directory);
    end
    try
        exportgraphics(gcf,output_file,'Resolution',220);
    catch
        print(gcf,output_file,'-dpng','-r220');
    end
    close(gcf);
end
