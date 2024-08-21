function create_all_figs(varargin)

%%% USER GUIDE ------------------------------------------------------------

% Inputs:
% folder_name (optional): name of folder where sensitivity data is stored.
%                       The default folder name is 'sensitivity_data'.
%                       Folder name must match folder name used in
%                       run_all.m function.


folder_name = 'sensitivity_data';
if nargin == 1
    folder_name = varargin{1};
end
if not(isfolder(folder_name))
    fprintf(['The folder named "', folder_name, '" does not exist in current directory.\n' ...
        ' If a folder name other than "sensitivity_data" was used to create a folder using run_all.m,\n' ...
        ' then make sure the same folder name is used an an input to the create_all_figs() function\n']);
end

% Get the files in the folder
files = dir(fullfile(folder_name, '*.mat'));

for f = 1:numel(files)
    
    % Get the data from the file into the workspace
    input_file_name = files(f).name;

    % Determine if computer is windows, mac, or linux
    if ispc
        slash = '\';
    else
        slash = '/';
    end
    relative_file_path = [folder_name, slash, input_file_name];
    structure = load(relative_file_path);
    
    output = structure.output;

    % Create figures
    [fig] = sens_fig(output);
    
    % Save figures in folder
    % Create new directory
    new_folder = 'sensitivity_figs';
    if not(isfolder(new_folder))
        mkdir(new_folder);
    end

    % Create figure file name using input file name
    [~, name, ~] = fileparts(input_file_name);
    file_name = name;

    % Check if export_fig is installed
    addons = matlab.addons.installedAddons;
    export_fig_installed = any(strcmp(addons.Name, 'export_fig'));
    
    % Save figs in the new folder
    if ~export_fig_installed
        fprintf('export_fig not installed. To save sensitivity figures as .pdf use export_fig.\n');
    
        fig_file_name = fullfile(new_folder, sprintf('%s.fig', file_name));
        savefig(fig_file_name);
    
    else
        fig_file_name = fullfile(new_folder, sprintf('%s.pdf', file_name));
        export_fig(fig, fig_file_name, '-pdf');
    end
end
end
