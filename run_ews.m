function [output] = run_ews(data, varargin)

% USER GUIDE --------------------------------------------------------------
%
% Computes early warning signals (EWS) and plots them with corresponding Kendall's
% tau values. Alternatively, performs a sensitivity analysis of the
% robustness of the EWS to the detrending parameters.
%
% This function uses name-value pair inputs for optional inputs, for example:
%
%   run_ews(data, 'temp_detrending', 'gaussian'),
%
% will compute and plot the early warning signals for the inputted data,
% whilst also changing the default setting for temporal detrending from
% linear to gaussian detrending.
%
%   [output] = run_ews(data, 'sens', 'yes', 'use_parallel', 'yes', 'num_cores', 6)
% 
% will perform a sensitivity analysis on the inputted data for the default detrending options using MATLAB's parallel computing
% with 6 cores and store the sensitivity and significance data in the output. To see the default options use:
%
%   run_ews()
%
% with no input arguments.

% INPUTS
% 
% bandwidth_percent:    Percentage of total length of the data used as a 
%                       smoothing bandwidth when detrending.
%                       Can either be a scalar or a vector containing the
%                       lower and upper values when performing a
%                       sensitivity analysis to different values.
%
% window_percent:       Percentage of total length of the data used as a
%                       rolling window when estimating EWS. Either a scalar
%                       or vector value.
%
% detrending:           Either 'gaussian', 'remove_mean', 'linear'.  
% -temp_detrending:     Type of detrending used to detrend the temporal
%                       data sets,
% -spat_detrending:     type of detrending used to detrend the spatial or
%                       multivariate data sets. 
% 
% bif_par_bounds:       If the data set used is produced by linearly
%                       varying a bifurcation parameter for example \mu, then the upper 
%                       and lower bounds can be specified as a vector,
%                       [\mu_{min}, \mu_{max}]. EWS will be plotted against
%                       the underlying bifurcation parameter.
%
% tau_bounds:           Compute Kendall's tau on a subset of the data. If
%                       an underlying bifurcation parameter is not used
%                       specify tau bounds as integer values, 
%                       [\tau_{min}, \tau_{max}]. For example, if the total
%                       number of time steps is 600, specifying tau bounds 
%                       as [100, 450] will compute Kendall's tau on this 
%                       interval of the data. If a bifurcation parameter is
%                       used then the tau bounds can be real
%                       valued. For example, if bif_par_bounds = [2, 4],
%                       then tau_bounds = [2.3, 3.1] will work.
%                           
% sens:                 Perform sensitivity analysis. Either 'yes', or
%                       'no'.
%
% res:                  Resolution of the sensitivity analysis. Computes
%                       the values of Kendall's tau for res^2 combinations of the
%                       bandwidth_percentage and window_percentage between the upper and lower
%                       bounds.
%
% plotting:             Plots EWS. Either 'yes' or 'no'.
%
% sens_plot:            Plots sensitivity figures. Either 'yes' or 'no'.
%
% use_parallel:         Use MATLAB's parfor when performing sensitivity
%                       analysis. Either 'yes' or 'no'.
%
% num_cores:            Number of logical cores used when creating a
%                       parallel pool of workers.

% OUTPUTS
%
% default_options:      run_ews(), with no input arguments.
% 
% Early warning signals and Kendall's tau values for a specific choice
% of parameters / and figure containing EWS and corresponding Kendall's tau values.
%
% Sensitivity analyis data / and figures.
% -------------------------------------------------------------------------

% Default options
default_options = struct('res', 8, ...
    'sens', 'no', ...
    'bandwidth_percent', [5, 50], ...
    'window_percent', [5, 50], ...
    'use_parallel', 'no', ...
    'num_cores', [], ...
    'sens_plot', 'no');

% Get default options from other functions
default_ews_func_options = ews_calculator();
default_Kendall_func_options = ews_fig();

% Update default options using default options from other functions
fields_ews = fieldnames(default_ews_func_options);
fields_Kendall = fieldnames(default_Kendall_func_options);
for i = 1:numel(fields_ews)
    if ~isfield(default_options, fields_ews{i})
        default_options.(fields_ews{i}) = default_ews_func_options.(fields_ews{i});
    end
end
for i = 1:numel(fields_Kendall)
    if ~isfield(default_options, fields_Kendall{i})
        default_options.(fields_Kendall{i}) = default_Kendall_func_options.(fields_Kendall{i});
    end
end

output = struct();
if nargin == 0
    output = default_options;

else

    % Update default options from user input
    if nargin > 1
        for i = 1:2:length(varargin) 
            default_options.(varargin{i}) = varargin{i+1}; 
        end
    end
    
    % Create cell array of arguments to pass into other functions
    args = create_arguments_cell(default_options);
    
    % Create variables from default options structure
    field_names = fieldnames(default_options);
    for i = 1:numel(field_names)
        variable_name = field_names{i};
        eval([variable_name ' = default_options.' variable_name ';']);
    end
    
    % Runs early warning signal calculator and computes statistics
    if strcmp(sens, 'no')
        [EWS] = ews_calculator(data, args{:});
        [statistics] = ews_fig(EWS, args{:});
        output.EWS = EWS;
        output.statistics = statistics;
    end
    
    % Performs sensitivity analysis for a choice of detrending methods
    if strcmp(sens, 'yes')
        
        % Prevents error (plotting the EWS res^2 times from the ews_fig function)
        if strcmp(plotting, 'yes')
            args_indx_plotting = find(strcmp(args, 'plotting'));
            args{args_indx_plotting + 1} = 'no';
        end
    
        % Create ranges of bandwidth and window size to vary over
        if numel(bandwidth_percent) ~=2 || numel(window_percent) ~=2
            error('Too few or too many inputs for bandwidth or window size. Each argument must be a vector containing the upper and lower bound when computing sensitivity data.'); 
        end
        bandwidth_range = linspace(bandwidth_percent(1), bandwidth_percent(2), res);
        window_size_range = linspace(window_percent(1), window_percent(2), res);
    
        ews = ews_calculator(data);
        data_type = ews.params.data_type;
        num_EWS = num_EWS_func(data_type);

        % Initialise arrays
        sens_data = zeros(res, res, num_EWS);
        sig_data = zeros(res, res, num_EWS);

        for jj = 1:res
            
            args_indx_ws = find(strcmp(args, 'window_percent'));
            args{args_indx_ws + 1} = window_size_range(jj);
    
            % If bandwidth sensitivity exists, loop over bandwidth range
            if ~strcmp(temp_detrending, 'linear') || strcmp(spat_detrending, 'gaussian')    
                
                if strcmp(use_parallel, 'yes')

                    % Check to see if Parallel computing toolbox is
                    % installed
                    addons = matlab.addons.installedAddons;
                    parallel_installed = any(strcmp(addons.Name, 'Parallel Computing Toolbox'));
                    if ~parallel_installed
                        warning(['Parallel Computing Toolbox is not installed. To use parallel computing' ...
                            'install Parallel Computing Toolbox']);
                    end
    
                    % Remove bandwidth from arguments so parfor can work
                    args_indx_bw = find(strcmp(args, 'bandwidth_percent'));
                    if ~isempty(args_indx_bw)
                        args{args_indx_bw} = 'not_bandwidth';
                        args{args_indx_bw + 1} = [];
                    end
                    
                    % If a parallel pool doesnt exist, then make one
                    if isempty(gcp("nocreate"))
                        if ~isempty(num_cores)
                            numcores = num_cores;
                        else
                            warning('Maximum number of logical cores will be used.');
                            numcores = feature('numcores');
                        end
                        parpool(numcores)
                    end

                    parfor ii = 1:res
                        bw = bandwidth_range(ii);
                        [EWS] = ews_calculator(data, 'bandwidth_percent', bw, args{:});
                        [statistics] = ews_fig(EWS, 'bandwidth_percent', bw, args{:});
            
                        sens_data(ii, jj, :) = statistics.taus;
                        sig_data(ii, jj, :) = statistics.p_vals;

                    end
    
                elseif strcmp(use_parallel, 'no')
    
                    for ii = 1:res

                        args_indx_bw = find(strcmp(args, 'bandwidth_percent'));
                        args{args_indx_bw + 1} = bandwidth_range(ii);
              
                        [EWS] = ews_calculator(data, args{:});
                        [statistics] = ews_fig(EWS, args{:});
            
                        sens_data(ii, jj, :) = statistics.taus;
                        sig_data(ii, jj, :) = statistics.p_vals;
                    end

                end

            elseif strcmp(temp_detrending, 'linear') && ~strcmp(spat_detrending, 'gaussian')
            
                [EWS] = ews_calculator(data, args{:});
                [statistics] = ews_fig(EWS, args{:});
                
                for ii = 1:res
                    sens_data(ii, jj, :) = statistics.taus;
                    sig_data(ii, jj, :) = statistics.p_vals;
                end
            end
        end

        delete(gcp("nocreate"))
    
        % Create outputs
        output.sens_data = sens_data;
        output.sig_data = sig_data;
        output.window_size_range = window_size_range;
        output.bandwidth_range = bandwidth_range;
        output.data_type = data_type;
    
        % Plotting
        if ~strcmp(temp_detrending, 'linear') || strcmp(spat_detrending, 'gaussian')
            contour = 'yes';    
        elseif strcmp(temp_detrending, 'linear') && ~strcmp(spat_detrending, 'gaussian')
            contour = 'no';
        end
        output.contour = contour;

        if strcmp(sens_plot, 'yes')
            sens_fig(output);
        end

    end
end
end


% functions
function [args] = create_arguments_cell(default_options)
    args = {};
    field_names = fieldnames(default_options);
    for i = 1:numel(field_names)  
        val = default_options.(field_names{i});
        if isempty(val)
            args = [args, field_names{i}, {[]}];
        else
            args = [args, field_names{i}, val];
        end
    end
end

function [num_EWS] = num_EWS_func(data_type)
    if strcmp(data_type, 'spatial')
        num_EWS = 9;
    elseif strcmp(data_type, 'multivariate')
        num_EWS = 8;
    else
        num_EWS = 4;
    end
end