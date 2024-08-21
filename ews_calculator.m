function [EWS] = ews_calculator(data, varargin)

default_options = struct('bandwidth_percent', 5, ...
    'window_percent', 5, ...
    'temp_detrending', 'linear', ...
    'spat_detrending', 'remove_mean');

if nargin == 0
    output = default_options;
    EWS = output;
else

    % Update default options structure from user input.
    if nargin > 1
        for i = 1:2:length(varargin) 
            default_options.(varargin{i}) = varargin{i+1}; 
        end
    end

    % Check dimensions and determine 'data_type'
    num_dim = numel(size(data));
    if num_dim == 3
        default_options.data_type = 'spatial';
    elseif num_dim == 2
        % compute # elements in each dimension, if either equals one then
        % data_type = 'temporal'
        dim_lengths = size(data);
        if dim_lengths(1) == 1 || dim_lengths(2) == 1
            default_options.data_type = 'temporal';
        else
            default_options.data_type = 'multivariate';
        end
    end
    % fprintf('%s data type detected. Computing %s EWS.\n', default_options.data_type, default_options.data_type)

    data_type = default_options.data_type;
    % Adjust rows and columns (if needed)
    dim_lengths = size(data);
    % Assumes # time steps > spatial dimensions
    if strcmp(data_type, 'spatial')
        dim_lengths = size(data);
        if dim_lengths(3) < dim_lengths(1)
            data = permute(data, [3, 2, 1]);
        elseif dim_lengths(3) < dim_lengths(2)
            data = permute(data, [1, 3, 2]);
        end
    elseif strcmp(data_type, 'multivariate')
        if dim_lengths(2) < dim_lengths(1)
            data = data';
        elseif dim_lengths(1) == dim_lengths(2)
            warning(['Multivariate data type detected. Number of variables equals number of time steps.' ...
                ' Please check that input data contains variables along the rows, and time steps along the columns.']);
        end
    elseif strcmp(data_type, 'temporal')
        if dim_lengths(1) < dim_lengths(2)
            data = transpose(data);
        end
    end
    
    % Create variables using fields in default options
    field_names = fieldnames(default_options);
    for i = 1:numel(field_names)
        variable_name = field_names{i};
        eval([variable_name ' = default_options.' variable_name ';']);
    end
    
    % Prepare data and assign variables
    if strcmp(data_type, 'spatial')
        spatial_data = data;
        numSteps = length(spatial_data(1,1,:));
        N = length(spatial_data(:,1,1));
        temporal_data = transpose(mean(reshape(spatial_data, [], numSteps))); 
    elseif strcmp(data_type, 'multivariate')
        multivariate_data = data;
        temporal_data = transpose(mean(multivariate_data));
        N = length(multivariate_data(:, 1));
        % For purposes of dealing with different cases, let
        spatial_data = multivariate_data;
        % however distinguish use using 'data_type' variable.
    elseif strcmp(data_type, 'temporal')
        temporal_data = data;
    end
    
    % Compute some parameters
    total_time = numel(temporal_data);

    if numel(bandwidth_percent) > 1
        bandwidth_percent = bandwidth_percent(1);
    end
    if numel(window_percent) > 1
        window_percent = window_percent(1);
    end
    bandwidth = round(bandwidth_percent/100 * total_time);
    window_size = round(window_percent/100 * total_time);
    
    % Detrending data
    [temp_residuals] = temporal_detrending(temporal_data, temp_detrending, bandwidth);
     
    if strcmp(data_type, 'spatial') || strcmp(data_type, 'multivariate')
        [spatial_residuals] = spatial_detrending(data_type, spatial_data, temporal_data, spat_detrending, bandwidth, N);
    end
    
    % Compute EWS
    [temp_EWS] = temp_EWS_calc(window_size, temp_detrending, total_time, temp_residuals);
    early_warning_signals = temp_EWS;

    if strcmp(data_type, 'spatial') || strcmp(data_type, 'multivariate')
        [spat_EWS] = spat_EWS_calc(N, data_type, spat_detrending, total_time, spatial_residuals, temporal_data, window_size);
        early_warning_signals = [temp_EWS, spat_EWS];
    end
        
    % outputs
    EWS = struct();
    EWS.early_warning_signals = early_warning_signals;
    EWS.params = default_options;

end

end

% functions
% Detrend temporal data
function [temp_residuals] = temporal_detrending(temporal_data, detrending, bandwidth)

    total_time = length(temporal_data);
    if strcmp(detrending, 'linear')
        detrended_temp_data = zeros(total_time, 1);
    elseif strcmp(detrending, 'gaussian')
        detrended_temp_data = smoothdata(temporal_data, 'gaussian', bandwidth);
    elseif strcmp(detrending, 'remove_mean')
        detrended_temp_data = smoothdata(temporal_data, 'movmean', bandwidth);
    end
    temp_residuals = temporal_data - detrended_temp_data;
end

% Detrending spatial/multivariate data
function [spatial_residuals] = spatial_detrending(data_type, spatial_data, temporal_data, detrending, bandwidth, N)
    
    total_time = length(temporal_data);
    if strcmp(data_type, 'spatial') || strcmp(data_type, 'multivariate')
        if strcmp(data_type, 'spatial')
            reshaped_spatial_data = reshape(spatial_data, [], total_time);
        else
            reshaped_spatial_data = spatial_data;
        end
        if strcmp(detrending, 'linear')
            if strcmp(data_type, 'spatial')
                detrended_spat_data = zeros(N^2, total_time);
            else
                detrended_spat_data = zeros(N, total_time);
            end
        elseif strcmp(detrending, 'gaussian')
            detrended_spat_data = smoothdata(temporal_data, 'gaussian', bandwidth)';
        elseif strcmp(detrending, 'remove_mean')
            detrended_spat_data = mean(reshaped_spatial_data);
        end
        spatial_residuals = reshaped_spatial_data - detrended_spat_data;
    end
end

% Linear regression
function [linear_fit] = linear_regression(ith_window)
    window_size = length(ith_window);
    X = [ones(window_size, 1), (1:window_size)'];
    coeffs = (X' * X) \ (X' * ith_window);
    linear_fit = coeffs(1) + coeffs(2).*(1:window_size);
end


% compute temporal early warning signals
function [temp_EWS] = temp_EWS_calc(window_size, detrending, total_time, temp_residuals)
    num_temp_EWS = 4;
    temp_EWS = zeros(total_time, num_temp_EWS);
    for i = window_size:total_time
        ith_window = temp_residuals(i+1-window_size:i);

        if strcmp(detrending, 'linear')
            linear_fit = linear_regression(ith_window);
            ith_window = ith_window - linear_fit';
        end
        n = window_size;
        % std
        temp_EWS(i,1) = std(ith_window);
        % skewness
        temp_EWS(i,2) = (sqrt((n - 1) ./ n) .* n ./ (n - 2)) .* mean((ith_window - mean(ith_window)).^3)./mean((ith_window - mean(ith_window)).^2).^1.5;
        % 'R' type autocorrelation (Taken from generic_ews code)
        s = var(ith_window);
        mu = mean(ith_window);
        lag = 1;
        Xt = ith_window(1:end - lag);
        Xtk = ith_window(1 + lag:end);
        temp_EWS(i,3) = 1 / (n - 1) * ( sum((Xt - mu) .* (Xtk - mu)) / s);
        % Autoregression
        data = ith_window(2:end);
        data_lag1 = ith_window(1:end-1);
        Matrix = [ones(length(data_lag1),1), data_lag1];
        coeffs = (Matrix' * Matrix) \ (Matrix' * data);
        temp_EWS(i,4) = coeffs(2);
    end
end

% compute spatial early warning signals
function [spat_EWS] = spat_EWS_calc(N, data_type, detrending, total_time, spatial_residuals, temporal_data, window_size)
            
    if strcmp(data_type, 'multivariate')
        num_spatial_EWS = 4;
    elseif strcmp(data_type, 'spatial')
        num_spatial_EWS = 5;
    end
    spat_EWS = zeros(total_time, num_spatial_EWS);

    % Compute eigenvalues of covariance matrix
    for i = window_size : total_time
        ith_spat_window = (spatial_residuals(:, i+1-window_size:i))';
        ith_temp_window = temporal_data(i+1-window_size:i);

        if strcmp(detrending, 'linear')
            % In this case, the temporal data is the spatial mean. 
            linear_fit = linear_regression(ith_temp_window);
            % Uniformly remove the linearly fit to the spatial mean from
            % each lattice site across a window
            ith_spat_window = ith_spat_window - linear_fit';
        end
        covariance_matrix = cov(ith_spat_window);
        [~, D] = eig(covariance_matrix);
        largest_eig = max(diag(D));
        spat_EWS(i,1) = largest_eig;
        spat_EWS(i,2) = max(diag(D))/ sqrt(sum( diag(D).^2 ));
    end
    % Spatial Variance and Skewness
    for i = 1:total_time
        if strcmp(data_type, 'spatial')
            spatial_mean = sum(spatial_residuals(:,i))/(N^2);
            spatial_variance = (1/N^2) * sum( (spatial_residuals(:,i) - spatial_mean).^2 );
            spatial_std = sqrt(spatial_variance);
            spatial_skewness = (1/N^2) * (1/spatial_std^3) * sum((spatial_residuals(:,i) - spatial_mean).^3);
        elseif strcmp(data_type, 'multivariate')
            spatial_mean = sum(spatial_residuals(:,i))/(N);
            spatial_variance = (1/N) * sum( (spatial_residuals(:,i) - spatial_mean).^2 );
            spatial_std = sqrt(spatial_variance);
            spatial_skewness = (1/N) * (1/spatial_std^3) * sum((spatial_residuals(:,i) - spatial_mean).^3);
        end
        spat_EWS(i,3) = spatial_variance;
        spat_EWS(i,4) = spatial_skewness;
    end
    % Spatial Correlation
    if strcmp(data_type, 'spatial')
        for ii = 1:total_time
            X = spatial_residuals(:,ii);
            spatial_mean = sum(X)/(N^2);
            product = (X - spatial_mean) * (X - spatial_mean)';
            neighbour_matrix = zeros(N^2,N^2);
            [i,j] = meshgrid(1:N^2, 1:N^2);
            % For periodic boundary conditions
            neighbours = (mod(i+1, N^2) == mod(j, N^2)) | (mod(i-1, N^2) == mod(j, N^2)) | (mod(i+N, N^2) == mod(j, N^2)) | (mod(i-N, N^2) == mod(j, N^2));
            neighbour_matrix(neighbours) = 1;
            unsummed_numerator = neighbour_matrix .* product;
            unsummed_denomenator = (X - spatial_mean).^2;
            unormalized_spatial_correlation = sum(sum(unsummed_numerator)) / sum(unsummed_denomenator);     
            normalization = (N^2) / sum(sum(neighbour_matrix));
            spatial_correlation = normalization * unormalized_spatial_correlation;
            spat_EWS(ii, 5) = spatial_correlation;
        end
    end
end