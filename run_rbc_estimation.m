%% Euro-area RBC Bayesian estimation extension
% Run this file from MATLAB after Dynare has been added to the MATLAB path.

clear;
clc;
close all;

project_directory = fileparts(mfilename('fullpath'));
cd(project_directory);

if exist('dynare','file') ~= 2
    error(['Dynare is not on the MATLAB path. Start MATLAB through Dynare ', ...
           'or add the Dynare matlab folder before running this file.']);
end

required_files = {'rbc_estimation.mod','euro_area_rbc_dynare.csv'};
for file_index = 1:length(required_files)
    if exist(required_files{file_index},'file') ~= 2
        error('Missing required file: %s',required_files{file_index});
    end
end

fprintf('\nRunning Bayesian estimation of the RBC model on euro-area data...\n');
dynare rbc_estimation.mod noclearall

fprintf('\nEstimation run completed.\n');
