function [spatial_data, reshaped_spatial_data, temporal_data] = synthetic_data(model_name, N, epsilon, sigma, required_size, kappa_max, plotting_value, varargin)

%%% USER GUIDE
%
% - model_name: select the model out of 'MODEL1', ... , 'MODEL9', 'harvest'
%
% - N: length of lattice (N by N)
%
% - epsilon: time scale seperation between fast and slow subsytem
%
% - sigma: standard deviation of noise
%
% - required_size: specifies the size of the subset of data produced from
%       the Euler-Maruyama method used as the final number of points within the
%       output. For example: required_size = 2000.
%
% - kappa_max: maximum value of underlying bifurcation parameter.
%
% - plotting_value: 0 or 1, plots snapshots of solution and the spatial
%       mean in a single figure.
%
% - varargin: Use varargin to specify kappa_min, the minimum value of the
%       underlying bifurcation parameter. 

% Example:
%
%   [spatial_data, multivariate_data, temporal_data] = synthetic_data('MODEL1', 20, 0.01, 0.01, 2000, 1.2, 1);
%
% Produces a 20 x 20 lattice of spatial data with 2000 time steps, for the
% first model with a normally distibuted white noise with standard
% deviation 0.01, for a fast-slow system with time scale seperation epsilon
% = 0.01, for kappa upto 1.2. 
% 
%--------------------------------------------------------------------------

    rng(910);
    % Parameters
    P_b0 = 7.25;
    k = 14.065;
    A = 0.6077;
    P_I = 0.9643;
    % sigma = 0.2;
    dR_s = 0.25;
    R_s = 0.5;
    varepsilon = 0.092;
    P_1 = (varepsilon - R_s)/dR_s + P_I;
    P_2 = (1 - R_s)/dR_s + P_I;

    kappa_min = 0;
    if nargin == 8
        kappa_min = varargin{1};
    end

    % Time
    tf = (kappa_max-kappa_min)/epsilon;
    dt = 0.01;
    time = 0: dt: tf;
    numSteps = length(time);
    kappa_vals = kappa_min : epsilon * dt : kappa_max;   
    % Initialize array
    r = zeros(N,N, numSteps);
    
    % MODELS
    % Chris's model
    if strcmp(model_name, 'MODEL1')
        heq_vals = Homogeneous_Equilibrium(model_name, 'sigmoidal', kappa_vals, numSteps);
        r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps-1
            kappa = kappa_vals(i);
            r_t = r(:,:,i);
            P_b = P_b0*N^2/sum(sum(r_t.^4));
            local_coupling = r_t.^4 + circshift(r_t, [0 1]).^4 + circshift(r_t, [0 -1]).^4 + circshift(r_t, [1 0]).^4 + circshift(r_t, [-1 0]).^4;
            P = P_b - kappa * k./r_t + P_b*A*(1-r_t + 1.5*(1-r_t).^2).*local_coupling;
            R = (1 + exp(-(P-P_I))).^(-1);
            drdt = R - r_t;
            r(:,:,i+1) = r(:,:,i) + drdt*dt + sigma*sqrt(dt)*randn(N,N);
        end   
    end

    % Replacing the luminal radius-transmural pressure relationship (everything else is the same) with piecewise linear relationship
    if strcmp(model_name, 'MODEL2')
        heq_vals = Homogeneous_Equilibrium(model_name, 'linear', kappa_vals, numSteps);
        r(:,:,1) = heq_vals(1); %+ sigma*randn(N,N);  
        % Euler-Maruyama
        for i = 1:numSteps-1
            kappa = kappa_vals(i);
            r_t = r(:,:,i);
            P_b = P_b0*N^2/sum(sum(r_t.^4));
            local_coupling = r_t.^4 + circshift(r_t, [0 1]).^4 + circshift(r_t, [0 -1]).^4 + circshift(r_t, [1 0]).^4 + circshift(r_t, [-1 0]).^4;
            P = P_b - kappa * k./r_t + P_b*A*(1-r_t + 1.5*(1-r_t).^2).*local_coupling;
            R = zeros(size(P));
            condition1 = (P < P_1);
            condition2 = (P >= P_1) & (P <= P_2);
            condition3 = (P > P_2);
            R(condition1) = varepsilon;
            R(condition2) = dR_s * (P(condition2)-P_I) + R_s;
            R(condition3) = 1;
            drdt = R - r_t;
            r(:,:,i+1) = r(:,:,i) + drdt*dt + sigma*sqrt(dt)*randn(N,N);   
        end
    end
    
    % Replacing the local coupling with linearized local coupling (WITHOUT collecting terms into "reaction-diffusion" type equation) (MODEL 3)
    if strcmp(model_name, 'MODEL3')
        heq_vals = Homogeneous_Equilibrium(model_name, 'linear', kappa_vals, numSteps);
        r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps-1
            kappa = kappa_vals(i);
            r_t = r(:,:,i);
            heq = heq_vals(i);
            P_b = P_b0*N^2/sum(sum(r_t.^4));
            local_coupling = r_t.^4 + 4*heq^3 *(circshift(r_t, [0 1]) + circshift(r_t, [0 -1]) + circshift(r_t, [1 0]) + circshift(r_t, [-1 0])) -12 * heq^4;
            P = P_b - kappa * k./r_t + P_b*A*(1-r_t + 1.5*(1-r_t).^2).*local_coupling;
            R = zeros(size(P));
            condition1 = (P < P_1);
            condition2 = (P >= P_1) & (P <= P_2);
            condition3 = (P > P_2);
            R(condition1) = varepsilon;
            R(condition2) = dR_s * (P(condition2)-P_I) + R_s;
            R(condition3) = 1;
            drdt = R - r_t;
            r(:,:,i+1) = r(:,:,i) + drdt*dt + sigma*sqrt(dt)*randn(N,N);
        end
    end

    % Linearized local coupling and rewriting as a "reaction-diffusion" type system (MODEL 4)
    if strcmp(model_name, 'MODEL4')
        heq_vals = Homogeneous_Equilibrium(model_name, 'linear', kappa_vals, numSteps);
        r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps-1
            r_t = r(:,:,i);
            heq = heq_vals(i);
            kappa = kappa_vals(i);
            P_b = (P_b0*N^2)/(sum(sum(r_t.^4)));
            Discrete_Laplacian = (circshift(r(:,:,i), [0 1]) + circshift(r(:,:,i), [0 -1]) + circshift(r(:,:,i), [1 0]) + circshift(r(:,:,i), [-1 0]) - 4 .* r(:,:,i));
            Local_coupling = r_t.^4 + 4.*heq.^3 .* (Discrete_Laplacian + 4.*r_t - 3.*heq);
            Quadratic_term = (1-r_t + 1.5.*(1-r_t).^2);
            P = P_b - kappa *(k./r_t) + P_b.*A.*Local_coupling.*Quadratic_term;
            condition1 = (P < P_1);
            condition2 = (P >= P_1) & (P <= P_2);
            condition3 = (P > P_2);
            Diffusion = zeros(size(P));
            Diffusion(condition1) = 0;
            Diffusion(condition2) = dR_s * P_b * A * Quadratic_term(condition2) * 4 * heq^3;
            Diffusion(condition3) = 0;
            Reaction = zeros(size(P));
            Reaction(condition1) = varepsilon - r_t(condition1);
            Reaction(condition2) = dR_s*(P_b - kappa*k./r_t(condition2) - P_I + P_b*A*(1-r_t(condition2)+1.5*(1-r_t(condition2)).^2).*(r_t(condition2).^4 + 4*heq^3 * (4*r_t(condition2) - 3*heq))) + R_s - r_t(condition2);
            Reaction(condition3) = 1 - r_t(condition3);
            dr = Diffusion .* Discrete_Laplacian + Reaction;
            r(:,:,i+1) = r(:,:,i) + dr*dt + sigma*sqrt(dt)*randn(N,N);
        end
    end

    % Diffusion for all x and t (non-constant) (MODEL5)
        if strcmp(model_name, 'MODEL5')
            heq_vals = Homogeneous_Equilibrium(model_name, 'linear', kappa_vals, numSteps); 
            r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
            % Euler-Maruyama
            for i = 1:numSteps-1
                r_t = r(:,:,i);
                heq = heq_vals(i);
                kappa = kappa_vals(i);
                P_b = (P_b0*N^2)/(sum(sum(r_t.^4)));
                Discrete_Laplacian = (circshift(r(:,:,i), [0 1]) + circshift(r(:,:,i), [0 -1]) + circshift(r(:,:,i), [1 0]) + circshift(r(:,:,i), [-1 0]) - 4 .* r(:,:,i));
                Local_coupling = r_t.^4 + 4.*heq.^3 .* (Discrete_Laplacian + 4.*r_t - 3.*heq);
                Quadratic_term = (1-r_t + 1.5.*(1-r_t).^2);
                P = P_b - kappa *(k./r_t) + P_b.*A.*Local_coupling.*Quadratic_term;
                condition1 = (P < P_1);
                condition2 = (P >= P_1) & (P <= P_2);
                condition3 = (P > P_2);
                Diffusion = dR_s * P_b * A * Quadratic_term * 4 * heq^3;        
                Reaction = zeros(size(P));
                Reaction(condition1) = varepsilon - r_t(condition1);
                Reaction(condition2) = dR_s*(P_b - kappa*k./r_t(condition2) - P_I + P_b*A*(1-r_t(condition2)+1.5*(1-r_t(condition2)).^2).*(r_t(condition2).^4 + 4*heq^3 * (4*r_t(condition2) - 3*heq))) + R_s - r_t(condition2);
                Reaction(condition3) = 1 - r_t(condition3);
                dr = Diffusion .* Discrete_Laplacian + Reaction;
                r(:,:,i+1) = r(:,:,i) + dr*dt + sigma*sqrt(dt)*randn(N,N);
            end
        end

    % Linearized local coupling with constant diffusion (MODEL 6)
    if strcmp(model_name, 'MODEL6')        
        heq_vals = Homogeneous_Equilibrium(model_name, 'linear', kappa_vals, numSteps);
        r(:,:,1) = heq_vals(1); % + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps-1
            r_t = r(:,:,i);    
            heq = heq_vals(i);
            kappa = kappa_vals(i);
            P_b = P_b0*N^2/sum(sum(r_t.^4));
            local_coupling = r_t.^4 + 4*heq^3 *(circshift(r_t, [0 1]) + circshift(r_t, [0 -1]) + circshift(r_t, [1 0]) + circshift(r_t, [-1 0])) -12 * heq^4;
            Quadratic_term = (1-r_t + 1.5.*(1-r_t).^2);
            P = P_b - kappa * k./r_t + P_b*A*(1-r_t + 1.5*(1-r_t).^2).*local_coupling;
            condition1 = (P < P_1);
            condition2 = (P >= P_1) & (P <= P_2);
            condition3 = (P > P_2);    
            Reaction = zeros(size(P));
            Reaction(condition1) = varepsilon - r_t(condition1);
            Reaction(condition2) = dR_s * (P_b - (kappa *k)./r_t(condition2) + P_b *  A * Quadratic_term(condition2) .* (r_t(condition2).^4 + 4 * heq^3 * (4*r_t(condition2) - 3*heq)) - P_I) + R_s - r_t(condition2);
            Reaction(condition3) = 1 - r_t(condition3);    
            Discrete_Lapalacian = circshift(r_t, [1 0]) + circshift(r_t, [-1 0]) + circshift(r_t, [0 1]) + circshift(r_t, [0 -1]) - 4*r_t;    
            dr = 1 * Discrete_Lapalacian + Reaction;
            r(:,:,i+1) = r(:,:,i) + dr*dt + sigma*sqrt(dt)*randn(N,N);
        end
    end

    % No piecewise function-constant diffusion (MODEL 7)
    if strcmp(model_name, 'MODEL7')
        heq_vals = Homogeneous_Equilibrium(model_name, 'linear', kappa_vals, numSteps);
        r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps - 1
            r_t = r(:,:,i);    
            kappa = kappa_vals(i);
            P_b = (P_b0 * (N)^2)/(sum(sum(r_t.^4)));       
            heq = heq_vals(i);   
            Discrete_Laplacian = (circshift(r(:,:,i), [0 1]) + circshift(r(:,:,i), [0 -1]) + circshift(r(:,:,i), [1 0]) + circshift(r(:,:,i), [-1 0]) - 4 .* r(:,:,i));     
            Quadratic_term = (1-r_t + 1.5.*(1-r_t).^2);      
            Diffusion_matrix = 1;
            Reaction_matrix = dR_s * (P_b - (kappa *k)./r_t + P_b *  A * Quadratic_term .* (r_t.^4 + 4 * heq^3 * (4*r_t - 3*heq)) - P_I) + R_s - r_t;    
            dr = Diffusion_matrix .* Discrete_Laplacian + Reaction_matrix;   
            r(:,:,i+1) = r(:,:,i) + dt .* dr + sigma .* sqrt(dt) .* randn(N, N);
            r(:,:,i+1) = max(r(:,:,i+1), varepsilon);
        end
    end

    % Chris's Model with constant base pressure
    if strcmp(model_name, 'MODEL8')
        heq_vals = Homogeneous_Equilibrium(model_name, 'constant', kappa_max, numSteps);
        r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps-1
            r_t = r(:,:,i);
            P_b = P_b0;
            local_coupling = r_t.^4 + circshift(r_t, [0 1]).^4 + circshift(r_t, [0 -1]).^4 + circshift(r_t, [1 0]).^4 + circshift(r_t, [-1 0]).^4;
            P = P_b - epsilon*time(i) * k./r_t + P_b*A*(1-r_t + 1.5*(1-r_t).^2).*local_coupling;
            R = (1 + exp(-(P-P_I))).^(-1);
            drdt = R - r_t;
            r(:,:,i+1) = r(:,:,i) + drdt*dt + sigma*sqrt(dt)*randn(N,N);
        end
    end

    % No piecewise function-constant diffusion (MODEL 9) with constant base pressure
    if strcmp(model_name, 'MODEL9')
        heq_vals = Homogeneous_Equilibrium(model_name, 'constant', kappa_max, numSteps);
        r(:,:,1) = heq_vals(1) + sigma*randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps - 1
            r_t = r(:,:,i);    
            P_b = P_b0;    
            heq = heq_vals(i);
            Discrete_Laplacian = (circshift(r(:,:,i), [0 1]) + circshift(r(:,:,i), [0 -1]) + circshift(r(:,:,i), [1 0]) + circshift(r(:,:,i), [-1 0]) - 4 .* r(:,:,i));     
            Quadratic_term = (1-r_t + 1.5.*(1-r_t).^2);   
            Diffusion_matrix = 1;
            Reaction_matrix = dR_s * (P_b - (epsilon* time(i) *k)./r_t + P_b *  A * Quadratic_term .* (r_t.^4 + 4 * heq^3 * (4*r_t - 3*heq)) - P_I) + R_s - r_t;   
            dr = Diffusion_matrix .* Discrete_Laplacian + Reaction_matrix;    
            if r(:,:,i) > 0
                r(:,:,i+1) = r(:,:,i) + dt .* dr + sigma .* sqrt(dt) .* randn(N, N); %+ sigma .* sqrt(dt) .* randn(N+1, N+1)
            else
                r(:,:,i+1) = 0;
            end
        end
    end

    % Harvesting model
    if strcmp(model_name, 'harvest')
        % sigma = 0.1;
        R = 0.2;
        K = 10;
        a = 1;    
        kappa_min = 2;
        Eqn = @(x) a*x.*(1-x./K) - kappa_min.*(x.^2 ./(x.^2 + 1));
        heq_val = fzero(Eqn, 7.3);      
        r(:,:,1) = heq_val; % + sigma * randn(N,N);
        % Euler-Maruyama
        for i = 1:numSteps-1
            X_t = r(:,:,i);
            c = kappa_vals(i);
            Discrete_Laplacian = (circshift(r(:,:,i), [0 1]) + circshift(r(:,:,i), [0 -1]) + circshift(r(:,:,i), [1 0]) + circshift(r(:,:,i), [-1 0]) - 4 .* r(:,:,i));
            dXdt = a*X_t .* (1 - X_t ./ K) - c .* (X_t.^2 ./ (X_t.^2 + 1)) + R .* Discrete_Laplacian;
            r(:,:,i+1) = r(:,:,i) + dXdt * dt + sigma * sqrt(dt) * randn(N,N);
        end
    end
    
    % Store outputs
    spatial_data = r;
    reshaped_spatial_data = reshape(spatial_data, [], numSteps);
    temporal_data = transpose(mean(reshape(r, [], numSteps)));

    % Check if tight_subplot() is installed
    addons = matlab.addons.installedAddons;
    tight_subplot_installed = any(strcmp(addons.Name, 'tight_subplot(Nh, Nw, gap, marg_h, marg_w)'));

    % Plotting
    if plotting_value == 1
        % Use subplot() for testing figures, use tight_subplot() (i.e ha =
        ... and axes(ha(i)) for producing nice figures, Note: tight_subplot() is
        % slower than suplot() )
        figure('WindowState', 'maximized');
        number_subplots = 5;
        columns = 5;
        rows = 1;

        % Create/initialize axes
        if tight_subplot_installed
            ha = tight_subplot(rows, columns+1, [0, 0.03], [0.05, 0.05], [0.05, 0.05]);
        else
            ax = zeros(columns+1,1);
        end
        
        for i = 1:number_subplots

            % subplot axes
            if tight_subplot_installed
                axes(ha(i));
            else
                ax(i) = subplot(rows, columns+1, i);
            end

            imagesc(r(:,:, 1 + (i-1) * floor(numSteps/(number_subplots-1))));
            colormap("turbo")
            clim([0 1]);
            if strcmp(model_name, 'MODEL1')
                title(['$\kappa = $', num2str(kappa_vals(1 + (i-1) * floor(numSteps/(number_subplots-1))))], 'Interpreter','latex', 'FontSize', 28);
            end
            axis xy
            pbaspect([1,1,1]);
        end

        % Plot the spatial mean
        if tight_subplot_installed
            axes(ha(columns+1));
        else
            ax(columns+1) = subplot(rows, columns+1, columns+1);
        end
        plot(kappa_vals, temporal_data, 'LineWidth', 2);

        xlabel('$\kappa$', 'Interpreter','latex', 'FontSize', 20);
        if strcmp(model_name, 'MODEL1')
            title('Spatial Mean', 'Interpreter','latex', 'FontSize', 28);
        end
        xlim([kappa_min, kappa_max]);
        ylim([0, 1]);
        pbaspect([1,1,1]);    

        % filename = sprintf('%s_numerical_sln.pdf', model_name);
        % export_fig(filename, '-pdf');

    elseif plotting_value == 0
        % Dont plot
    end

    if ~required_size == 0
        subset_indices = round(linspace(1, numSteps, required_size));
        temporal_data = temporal_data(subset_indices);
        spatial_data = spatial_data(:,:, subset_indices);
        heq_vals = heq_vals(subset_indices);
        kappa_vals = kappa_vals(subset_indices);
        reshaped_spatial_data = reshaped_spatial_data(:, subset_indices);
    elseif required_size == 0
        % dont take a subset of the data (this is for when epsilon is
        % small)
    end

    % kappa_upper = find_k_upper(kappa_vals, heq_vals, temporal_data, 4*sigma/N);

end

% Determine kappa_upper
% function [kappa_upper] = find_k_upper(kappa_vals, heq_vals, temporal_data, threshold)
%     diff_vec = abs(temporal_data - heq_vals);
%     exceed_threshold_indices = find(diff_vec > threshold);
%     % disp(exceed_threshold_indices);
%     upper_index = exceed_threshold_indices(1);
%     kappa_upper = kappa_vals(upper_index);
% end
   
function [heq_vals] = Homogeneous_Equilibrium(model_name, Transmural_Pressure, kappa_vals, numSteps)
    heq_vals = zeros(numSteps,1);
    % Parameters
    P_b0 = 7.25;
    k = 14.065;
    A = 0.6077;
    P_I = 0.9643;
    dR_s = 0.25;
    R_s = 0.5;   
    for t = 1:numSteps
        kappa = kappa_vals(t);
        if strcmp(Transmural_Pressure, 'constant')
            Eqn = @(x) (1 + exp(-P_b0 + k*kappa/x - P_b0*A*5*(1-x+1.5*(1-x)^2)*x^4 + P_I))^(-1) - x;
        elseif strcmp(Transmural_Pressure, 'sigmoidal')
            Eqn = @(x) (1 + exp(-P_b0/x^4 + k*kappa/x - P_b0*A*5*(1-x+1.5*(1-x)^2) + P_I))^(-1) - x;     
        elseif strcmp(Transmural_Pressure, 'linear')
            Eqn = @(x) dR_s*(P_b0/(x^4) - k*kappa/x + P_b0*A*5*(1-x+1.5*(1-x)^2) - P_I) + R_s - x;
        end
        if t == 1
            heq_vals(t) = fzero(Eqn, 0.98);
        else
            heq_vals(t) = fzero(Eqn, heq_vals(t-1));
        end
        if ~strcmp(model_name, 'MODEL7')
            heq_vals(t) = min(heq_vals(t), 1);
        end
    end
end
