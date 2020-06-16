function ver_figsave(handle,fig_name,variable_names,drugs,additional_info)
% saves figure using version information from the git 
% required inputs: handle, fig_name
% handle: hanlde of figure
% fig_name: what to save figure as
% variable_names: optional, variables used to make figure
% drugs: optiona, drugs used to make figure
% additional_info: any additional ionfo to save under extraInfo in the
% UserData
% input a figure handle, save UserData and save figure
% user can later use get_fig_info to easily access this metadata

    % chosen_workspaces is a global variable to make it more easily
    % accessible.  This avoids making it another input.
    global chosen_workspaces

    % use external function to get git info
    git_info = getGitInfo();

    % set up an empty struct to put in figure info
    figure_info = struct;
    
    % can save hash info in the UserData of figure
    try
        figure_info.hash =  git_info.hash; 
    catch
        disp("Unable to get git hash information")
    end 
    
    % save workspacee information
    if ~isempty(chosen_workspaces)
        figure_info.workspace = chosen_workspaces;
    end 
    
    % get variable names if given as input, save in variables 
    if exist('variable_names','var')
        figure_info.variables = variable_names; 
    end
    
    % get drugs if given as input, save in drugs
    if exist('drugs','var')
        figure_info.drugs = drugs; 
    end
    
    % if there is additional description, add it to extraInfo
    if exist('additional_info','var')
        figure_info.extraInfo = additional_info;
    end 
    
    % want to get the time of creation of figure
    figure_info.createdDate = datetime('now');

    % set the figure_info to the UserData
    handle.UserData = figure_info;
    
    % get hex time to add to figure name to get a unique figure name to
    % prevent saving oveer
    hextime = get_hex_time();
    
    % save figure with given figure name 
    savefig(handle,strcat(fig_name,hextime))
end

