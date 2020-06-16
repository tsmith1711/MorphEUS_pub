function workspace_list = get_workspaces(directory)
%GET_WORKSPACES Gets a list of workspaces in a directory
%   Returns a horizontal cell array of file names of all workspaces in the
%   directory (relative to CALLER SCRIPT, not this function)

% add the path to the workspaces folder
addpath(directory);

% get all workspaces in the workspace folder
workspace_list = dir(strcat(directory, '/*.mat'));
workspace_list = {workspace_list.name};

end

