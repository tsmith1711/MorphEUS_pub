function fig_info = get_fig_info(fig_handle)
% easy way to look at figure UserInfo obtained 

    % if there was not a figure inputed, default to gcf
    if ~exist('fig_handle','var')
       fig_handle = gcf; 
    end

    % want to get the UserData
    fig_info = fig_handle.UserData;
end

