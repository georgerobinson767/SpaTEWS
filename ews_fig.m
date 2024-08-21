function [statistics] = ews_fig(EWS_structure, varargin)

default_options = struct('bif_par_bounds', [], ...
    'tau_bounds', [], ...
    'plotting', 'yes');

if nargin == 0
    output = default_options;
    statistics = output;
else

    if nargin > 1
        for i = 1:2:length(varargin) 
            default_options.(varargin{i}) = varargin{i+1}; 
        end
    end
    
    % EWS variables
    ews = EWS_structure.early_warning_signals;
    total_time = length(ews(:,1));
    data_type = EWS_structure.params.data_type;

    num_EWS = num_EWS_func(data_type);
    
    window_percent = EWS_structure.params.window_percent;
    
    % This functions variables
    field_names = fieldnames(default_options);
    for i = 1:numel(field_names)
        variable_name = field_names{i};
        expression = [variable_name '= default_options.' variable_name ';'];
        eval(expression);
    end

    % % Alternatively to eval()
    % bif_par_bounds = default_options.bif_par;
    % tau_bounds = default_options.tau_bounds;
    % plotting = default_options.plotting;

    % Errors
    if ~isempty(tau_bounds)
        if isempty(bif_par_bounds) && ~isinteger(tau_bounds)
            error(['Bifurcation parameters have not been specified and the bounds for computing Kendalls tau are not integers.' ...
                ' Either specify the bifurcation parameter bounds if using them, or use integer values for the tau bounds.']);
        end
    end

    if ~isempty(tau_bounds) && ~isempty(bif_par_bounds)
        if (tau_bounds(1) < bif_par_bounds(1)) || (bif_par_bounds(2) < tau_bounds(2))
            error('Bounds for computing Kendalls tau exceed bifurcation parameter bounds.')
        end
    end

    if numel(window_percent) > 1
        window_percent = window_percent(1);
    end
    window_size = round((window_percent/100)*total_time);
    
    % Compute Kendall's tau values for EWS ------------------------------------
    [indx_ws, indx] = create_index(bif_par_bounds, tau_bounds, window_size, total_time);
    
    % Initialise array
    taus = zeros(num_EWS, 1);
    H_vals = zeros(num_EWS, 1);
    p_vals = zeros(num_EWS, 1);
    
    for i = 1:num_EWS
        if i <= 6
            ord_vec = indx_ws(1):indx_ws(2);
        else
            ord_vec = indx(1):indx(2);
        end
    
        % Kendall's tau
        taus(i) = corr(ord_vec', ews(ord_vec, i), 'Type','Kendall', 'Rows','complete');
    
        % Modified Mann-Kendall test
        alpha = 0.05;
        alpha_ac = alpha;
        [p, H] = Modified_MannKendall_test(ews(ord_vec, i), alpha, alpha_ac);
        p_vals(i) = p;
        H_vals(i) = H;
    end
    
    taus = round(taus .* 10^3) / 10^3;
    p_vals = round(p_vals .* 10^3) / 10^3;
    
    % Outputs
    statistics = struct();
    statistics.taus = taus;
    statistics.p_vals = p_vals;
    statistics.H_vals = H_vals;
    statistics.num_EWS = num_EWS;
    
    if strcmp(plotting, 'yes')
        
        % Determine number of subplots based on data_type
        figure('WindowState', 'maximized');
        % figure;
        if strcmp(data_type, 'spatial')
            rows = 3;
            columns = 3;
        elseif strcmp(data_type, 'multivariate')
            rows = 2;
            columns = 4;
        elseif strcmp(data_type, 'temporal')
            rows = 2;
            columns = 2;
        end
        %ha = tight_subplot(rows, columns, [0.1, 0.05], 0.05, 0.05);
        
        EWSignals = {'Standard Deviation', 'Skewness', 'acf', 'AR(1)',...
            '$\sigma_1$', '$\sigma_1 / \sqrt{\sigma_1^2 + \cdots + \sigma_n^2}$', ...
            'Spatial Variance', 'Spatial Skewness', 'Spatial Correlation'};   
    
        for i = 1:num_EWS
            %axes(ha(i));
            
            if i <= 6
                ord_vec = indx_ws(1):indx_ws(2);
                ord_vec_total = window_size:total_time;
            else
                ord_vec = indx(1):indx(2);
                ord_vec_total = 1:total_time;
            end
            
            % Create time vectors from ordinal vectors where fowards: time (\in mathbb{z}) -> bifurcation parameter (\in \mathbb{R}) (if nessecary)
            time_vec = ord_vec;
            time_vec_total = ord_vec_total;
            time_full = 1:total_time;
            if ~isempty(bif_par_bounds)
                time_vec = fowards(time_vec, bif_par_bounds, total_time);
                time_vec_total = fowards(time_vec_total, bif_par_bounds, total_time);
                time_full = fowards(time_full, bif_par_bounds, total_time);
            end
    
            % Plot
            subplot(rows, columns, i);
            if i == 1 || i == 5 || i == 7   
                semilogy(time_vec_total', ews(ord_vec_total, i), 'b', 'LineWidth', 1);
                hold on
                semilogy(time_vec', ews(ord_vec, i), 'r', 'LineWidth', 2);
                hold off
            else
                plot(time_vec_total', ews(ord_vec_total, i), 'b', 'LineWidth', 1);
                hold on
                plot(time_vec', ews(ord_vec, i), 'r', 'LineWidth', 2);
                hold off
            end
    
            x_min = time_full(1);
            x_max = time_full(end);
            xlim([x_min, x_max]);
            ylabel(EWSignals{i}, 'Interpreter','latex', 'fontsize', 16);
            title(['\tau = ', num2str(taus(i))], 'FontSize',24);
            pbaspect([2,1,1]);
        end
    end
end
end


% Functions
% Invert tau indicies in 'bifurcation paramater' space into time indexes
% starting at 1 
function [indices] = invert(total_time, bounds, Y)
    bif_min = bounds(1);
    bif_max = bounds(2);
    slope = (bif_max - bif_min)/(total_time - 1);
    indices = (Y - bif_min)./slope + 1;
    indices = round(indices);
end

% convert integer indexing into the 'bifurcation parameter' space
function [indexed_bif_par] = fowards(X, bounds, total_time)
    bif_min = bounds(1);
    bif_max = bounds(2);
    slope = (bif_max - bif_min)/(total_time - 1);
    indexed_bif_par = slope.*(X-1) + bif_min;
end

% How many EWS are there?
function [num_EWS] = num_EWS_func(data_type)
    if strcmp(data_type, 'spatial')
        num_EWS = 9;
    elseif strcmp(data_type, 'multivariate')
        num_EWS = 8;
    else
        num_EWS = 4;
    end
end

% Compute indices based on several conditions
function [indx_ws, indx] = create_index(bif_par_bounds, tau_bounds, window_size, total_time)
    if ~isempty(tau_bounds)
        if ~isempty(bif_par_bounds)
            indices = invert(total_time, bif_par_bounds, tau_bounds); % need to find what the tau_bounds (in bifurcation parameter space) are w.r.t the indices [1, ..., t_f] for the EWS. 
        else
            indices = tau_bounds;
        end
        if indices(1) < window_size
            indx_ws = [window_size, indices(2)];
        else
            indx_ws = indices;
        end
        indx = indices;
    else
        indx_ws = [window_size, total_time];
        indx = [1, total_time];
    end
end