%% Use this demo.m file to explore different options available when using the functions in this package.

% **The expected time to complete this demo is approximately five minutes**


%% Create synthetic data sets
% This will produce a spatial, multivariate, and temporal data set using a
% lung ventilation model that undergoes spatial tipping, on a 10 x 10 grid
% with 400 time steps.
clear, clc, close all

t1 = tic;

time_steps = 400;

grid_length = 10;

[spatial_data, multivariate_data, temporal_data] = synthetic_data('MODEL1', grid_length, 0.01, 0.01, time_steps, 1.2, 0);


%% Compute and plot EWS using the run_ews() function
% Using the run_ews function without any additional inputs other than a data
% set will compute EWS using the default settings.

% Compute spatial EWS
run_ews(spatial_data); 

% Compute multivariate EWS
run_ews(multivariate_data);

% Compute temporal EWS
run_ews(temporal_data);


%% Perform a sensitivity analysis using the run_ews() function
% To perform a sensitivity analysis requires using the name-value pair input
%   ('sens', 'yes').
%
% Will perform a sensitivity analysis of the EWS computed for the inputted
% data against varying the bandwidth and window size. 

[output] = run_ews(spatial_data, 'sens', 'yes');
sens_fig(output);

% The same can also be achieved using the following with 
%   ('sens_plot', 'yes'), as an name-value pair input, for example:

% run_ews(spatial_data, 'sens', 'yes', 'sens_plot', 'yes');

% by default the run_ews function does not produce sensitivity figures.


%% Create sensitivity plots using test_data.mat
% Displays sensitivity plots computed for a 'large' spatial data set (20x20x2000) for 400
% combinations of window size and bandwidth ('res', 20).

load('test_data.mat');
[fig] = sens_fig(test_data);

%% Changing the default options
% To see the default options for the run_ews function, use the run_ews function with no input arguments

default_options = run_ews();
disp(default_options)

% For an explanation of the different inputs available, see run_ews.m

%% EXAMPLES
%% Change the 'resolution' of the sensitivity images
% The 'res' argument changes the number of bandwidth-windowsize pairs used
% to compute the EWS. The total number of pairs used is res^2, hence
% increasing the resolution increases the computational time. In this
% example we also demonstrate changing many of the other available options.

% Generate data
time_steps = 400;
grid_length = 10;
[spatial_data, ~] = synthetic_data('MODEL1', grid_length, 0.01, 0.01, time_steps, 1.2, 0);

[output] = run_ews(spatial_data, 'res', 3, 'sens', 'yes', 'spat_detrending', 'gaussian', 'temp_detrending', 'gaussian', 'bif_par_bounds', [0, 1.2], 'tau_bounds', [0.6, 1]);
sens_fig(output);

%% Use parallel computing to speed things up    **(Expected run time: approximately 1-2 minutes)**
% It is also possible to adjust the number of cores used using the
% ('num_cores', x) name-value pair input, by default the maximum number of
% logical cores will be used.

% Parallel computing will only be used if the Parallel Computing Toolbox is
% installed, otherwise run_ews.m will use a normal for loop instead of a
% parfor loop. For this demo, the parallel computing option will only be
% demonstrated if the Parallel Computing Toolbox is installed.
t2 = tic;
addons = matlab.addons.installedAddons;
parallel_installed = any(strcmp(addons.Name, 'Parallel Computing Toolbox'));
if parallel_installed

    run_ews(spatial_data, 'use_parallel', 'yes', 'sens_plot', 'yes', 'sens', 'yes', 'res', 6, 'spat_detrending', 'gaussian', 'temp_detrending', 'gaussian', 'bif_par_bounds', [0, 1.2], 'tau_bounds', [0.6, 1]);

end
t2_end = toc(t2);

%% Create data sets displaying spatial tipping for SLDS and SPDE model.
% Use the synthetic_data function to create figures of the data sets

synthetic_data('MODEL1', 100, 0.05, 0.01, 1000, 1.4, 1, 1);

synthetic_data('MODEL7', 100, 0.05, 0.01, 1000, 1.4, 1, 1);


%% Create sensitivity figures for all detrending methods    **(Expected run time: approximately 90s)**
% Compare sensitivity figures for all available detrending methods for a
% single data set

t3 = tic;
% Create a 'small' synthetic spatial data set
[spatial_data, ~] = synthetic_data('MODEL1', 10, 0.01, 0.01, 400, 1.2, 0);

detrending_methods = {'gaussian', 'linear', 'remove_mean'};
for i = 1:numel(detrending_methods)

    spat_detrending = detrending_methods{i};
    temp_detrending = detrending_methods{i};

    [output] = run_ews(spatial_data, ...
    'sens', 'yes', ...
    'bif_par_bounds', [0, 1.2], ...
    'tau_bounds', [0.6, 1.03], ...
    'window_percent', [3, 50], ...
    'bandwidth_percent', [3, 50], ...
    'res', 4, ...
    'sens_plot', 'no', ...
    'use_parallel', 'no', ...
    'spat_detrending', spat_detrending, ...
    'temp_detrending', temp_detrending);

    % Make folder to save sensitivity data and figures in
    folder_name = 'sensitivity_data';
    if not(isfolder(folder_name))
        mkdir(folder_name);
    end

    % Save sensitivity data in folder
    data_file_name = fullfile(folder_name, sprintf('SLDS_%s_sensitivity_data.mat', detrending_methods{i}));
    save(data_file_name, 'output', '-mat');

    % Create figures
    [fig] = sens_fig(output);

    % Save figures in folder
    addons = matlab.addons.installedAddons;
    export_fig_installed = any(strcmp(addons.Name, 'export_fig'));
    if ~export_fig_installed
        fprintf('export_fig not installed. To save sensitivity figures as .pdf use export_fig.\n');
        
        fig_file_name = fullfile(folder_name, sprintf('SLDS_%s_sensitivity_fig.fig', detrending_methods{i}));
        savefig(fig_file_name);

    else
        fig_file_name = fullfile(folder_name, sprintf('SLDS_%s_sensitivity_fig.pdf', detrending_methods{i}));
        export_fig(fig, fig_file_name, '-pdf');
    end
end
t3_end = toc(t3);
t1_end = toc(t1);


%% Demo run times

fprintf('Total run time to use parallel computing is %f seconds.\n', t2_end);
fprintf('Total run time to create all sensitivity figures is %f seconds.\n', t3_end);
fprintf('Total demo run time is %f minutes.\n', t1_end/60);
%% (ADDITIONAL) Use the run_all.m and create_all_figs.m functions to reproduce the data sets and some figures from our paper

% WARNING:  **Using the following will take a while when using few cores. The
%           following was designed to be left running on a computational server 
%           with at least 20 cores.**

% The run_all.m function can be used to create a folder containing the
% sensitivity data for each detrending method for the SLDS and SPDE models
% discussed in our paper. 

% The create_all_figs.m function can be used to make a folder
% containing the sensitivity figures for each sensitivity data created
% using the run_all.m function. 

% The run_all.m function is designed for use on a compuational server that 
% does not have the ability to produce graphics, hence the function creates a folder that
% can be transfered to a seperate computer where the create_all_figs.m file
% can be used to run the sens_fig.m function on each of the sensitivity
% data. These figures are then saved in a new folder in the current
% directory.

% For example, using the following:
%
% run_all(4, 'all_sens_data');
% create_all_figs('all_sens_data');
%
% would produce the sensitivity data for each model and for each detrending
% method using 4 cores and store the output in a folder named
% 'all_sens_data'. Then a folder named 'sensitivity_figs' would be created
% containing the sensitivity plots for each data set in the 'all_sens_data' 
% folder.
