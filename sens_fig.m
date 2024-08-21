function [fig] = sens_fig(input)

% Check if tight_subplot() is installed
addons = matlab.addons.installedAddons;
tight_subplot_installed = any(strcmp(addons.Name, 'tight_subplot(Nh, Nw, gap, marg_h, marg_w)'));

% Create variables from input
contour = input.contour;
ban_range = input.bandwidth_range;
win_range = input.window_size_range;
sig_data = input.sig_data;
sens_data = input.sens_data;
data_type = input.data_type;

% Array containing EWS name strings
EWSignals = {'Standard Deviation', 'Skewness', 'acf', 'AR(1)', '$\sigma_1$',...
    '$\sigma_1 / \sqrt{\sigma_1^2 + \dots + \sigma_n^2}$', 'Spatial Variance',...
    'Spatial Skewness', 'Spatial Correlation'};     

% Take absolute value of data to compare trend magnitude
abs_sens_data = abs(sens_data);

% Create figure
fig = figure('windowstate', 'maximized');
rows = 2;
if strcmp(data_type, 'spatial') || strcmp(data_type, 'multivariate')
    columns = 6;
else
    columns = 4;
end

% Create/initialize axes
if tight_subplot_installed
    ha = tight_subplot(rows, columns, [0, 0.05], 0.2, 0.1);
else
    ax = zeros(rows*columns, 1);
end

min_col = 0;
max_col = 1;
for j = 1:rows
    for i = 1:columns
        im_num = (j-1)*columns + i;

        % subplot axes
        if tight_subplot_installed
            axes(ha(im_num));
        else
            ax(im_num) = subplot(rows, columns, im_num);
        end
        
        % Create contour plots if bandwidth sensitivity is included
        if strcmp(contour, 'yes')
            [X,Y] = meshgrid(win_range, ban_range);

            % p-values
            if j == rows
                contourf(X,Y,sig_data(:,:,i), [0.05 1]);
                clim([0, 1]);

            % tau-values
            else
                contourf(X,Y,abs_sens_data(:,:,i));
                clim([min_col, max_col]);
                colormap("hot");
            end
            set(gca, 'FontSize', 20); 

        % Otherwise plot only window size sensitivity
        else
            y_lim_max = max(max(sig_data(1,:,:)));
            if j == rows
                plot(win_range, sig_data(1,:,i), 'r');
                hold on
                alpha = 0.05;
                plot(win_range, alpha .* ones(numel(win_range), 1), 'k', 'LineWidth',2.5);
                hold off
                ylim([0, y_lim_max]);
            else
                plot(win_range, abs_sens_data(1,:,i), 'b');
                ylim([min_col, max_col]);
            end
        end
        set(gca, 'FontSize', 16);
        if strcmp(contour, 'no')
            ylabels = {'$\tau$', '$p$'};
            if i == 1
                ylabel(ylabels{j}, 'Interpreter', 'latex', 'fontsize', 28);
            end
        end
        if j == 1
            title(EWSignals{i}, 'Interpreter', 'latex', 'fontsize', 22);
        end
        xlim([win_range(1), win_range(end)]);
        pbaspect([1,1,1]);
    end    
end

% Adjust colormap for the statistical significance contour plots
indices = columns+1:columns*rows;
if tight_subplot_installed
    set(ha(indices), 'colormap', [0,0,0;1,1,1]);
else
    set(ax(indices), 'colormap', [0,0,0;1,1,1]);
end

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';

if strcmp(contour, 'yes')
    ylabel(han,'Bandwidth (\%)', 'Interpreter', 'latex', 'fontsize', 22);
end
xlabel(han,'Window Size (\%)', 'Interpreter', 'latex', 'fontsize', 22);

if tight_subplot_installed
    han.YLabel.Position = [-0.09, 0.5, 0];
    han.XLabel.Position = [0.5, 0.08, 0];
else
    % subplot's label positions
    han.YLabel.Position = [-0.06, 0.5, 0];
    han.XLabel.Position = [0.5, 0.01, 0];
end

% Add color bar if neseccary
if strcmp(contour, 'yes')
    cbar = colorbar('Position', [0.92, 0.55, 0.015, 0.19]);
    cbar.FontSize = 14;
    clim([min_col, max_col]);

    % Create new set of axis for fig to artificially add a discrete
    % colorbar
    han2=axes(fig,'visible','off');

    black_white = [0,0,0; 1,1,1];
    han2.Colormap = black_white;
    cbar2 = colorbar('Position', [0.92, 0.25, 0.015, 0.19], "TickLabelInterpreter","tex");
    cbar2.FontSize = 14; 
    %cbar2.TickLabelInterpreter = "tex";
    cbar2.Ticks = [0.25, 0.75];
    cbar2.TickLabels = {'p>0.05', 'p<0.05'};

    if ~tight_subplot_installed
        cbar.Position = [0.92, 0.65, 0.015, 0.21];
        cbar2.Position = [0.92, 0.175, 0.015, 0.21];
        annotation('textbox',[0.95, 0.9, 0, 0],'string', '$\tau$','FitBoxToText','on', 'FontSize', 28,'EdgeColor','white', 'Interpreter', 'latex');
        annotation('textbox',[0.95, 0.41, 0, 0],'string', '$p$','FitBoxToText','on', 'FontSize', 28,'EdgeColor','white', 'Interpreter', 'latex');
    else
        annotation('textbox',[0.95, 0.775, 0, 0],'string', '$\tau$','FitBoxToText','on', 'FontSize', 28,'EdgeColor','white', 'Interpreter', 'latex');
        annotation('textbox',[0.95, 0.475, 0, 0],'string', '$p$','FitBoxToText','on', 'FontSize', 28,'EdgeColor','white', 'Interpreter', 'latex');
    end
end
end