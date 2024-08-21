function run_all(num_cores, varargin)

%%% USER GUIDE ------------------------------------------------------------

% Inputs:
% num_cores (optional): Number of logical cores used when creating a
%                       parallel pool of workers used to compute
%                       sensitivity data. If empty, run_all.m will use the
%                       maximum number of available cores to compute the
%                       sensitivity data.
%
% folder_name (optional): name of folder where sensitivity data is stored.
%                       The default folder name is 'sensitivity_data'.

Models = {'MODEL1', 'MODEL7'};

if nargin == 0
    num_cores = []; % Uses maximum number of available cores
end

for m = 1:numel(Models)

    % Create a synthetic spatial data set
    [spatial_data, ~] = synthetic_data(Models{m}, 20, 0.01, 0.01, 2000, 1.2, 0);

    % Create tau_bounds based on model
    if m == 1
        tau_bounds = [0.6, 1.03];
    elseif m == 2
        tau_bounds = [0.52, 0.95];
    end

    detrending_methods = {'gaussian', 'linear', 'remove_mean'};
    for i = 1:numel(detrending_methods)
    
        spat_detrending = detrending_methods{i};
        temp_detrending = detrending_methods{i};
    
        [output] = run_ews(spatial_data, ...
        'sens', 'yes', ...
        'bif_par_bounds', [0, 1.2], ...
        'tau_bounds', tau_bounds, ...
        'window_percent', [3, 50], ...
        'bandwidth_percent', [3, 50], ...
        'res', 20, ...
        'sens_plot', 'no', ...
        'use_parallel', 'yes', ...
        'num_cores', num_cores, ...
        'spat_detrending', spat_detrending, ...
        'temp_detrending', temp_detrending);
    
        % Make folder to save sensitivity data
        folder_name = 'sensitivity_data';
        if nargin == 2
            folder_name = varargin{1};
        end
        if not(isfolder(folder_name))
            mkdir(folder_name);
        end
    
        % Save sensitivity data in folder
        data_file_name = fullfile(folder_name, sprintf('%s_%s_sensitivity_data.mat', Models{m}, detrending_methods{i}));
        save(data_file_name, 'output', '-mat');

    end
end
end