function selection = checkboxList(promptstring, variables, defaults)
% Variables is horizontal cell array of variablenames
% defaults is a horizontal logical array of defaults for each variable

%% Check for bad input
if length(variables) ~= length(defaults)
    error("Variables and defaults list must be same length.")
end

%% Program
numvars = length(variables);

var_width = 130;
var_height = 15;
var_y_spacing = 2;

ss = get(0,'ScreenSize'); % Bottom left is 0, 0; top right is max x and y
fig_width = 140;
extra_fig_space = 70;
fig_height = ((var_height + var_y_spacing)*numvars) + extra_fig_space;
fig_x = (ss(3)/2) - (fig_width/2); 
fig_y = (ss(4)/2) - (fig_height/2);
fig_pos = [fig_x fig_y fig_width fig_height];

text_y = fig_height - 35;
text_x = 10;
text_height = 30;
text_pos = [text_x text_y fig_width-text_x text_height];

vars_x = 10;
vars_y = text_y - 20;

ok_width = 70;
ok_height = 20;
ok_x = (fig_width/2) - (ok_width/2);
ok_y = 5;
ok_pos = [ok_x ok_y ok_width ok_height];

% Create enclosing figure
h.f = figure('units','pixels','position', fig_pos,...
             'toolbar','none','menu','none');

% Create prompting test
prompt_text = uicontrol('Style','text','String',promptstring,...
        'HorizontalAlignment','center',...
        'Position', text_pos); %#ok        
       
% Create yes/no checkboxes
for i=1:numvars
    x = vars_x;
    y = vars_y - ((var_height + var_y_spacing) * (i - 1));
    
    h.c(i) = uicontrol('style','checkbox','units','pixels',...
                    'position',[x, y, var_width, var_height],'string',variables{i}, 'Value', defaults(i));
end

% Create OK pushbutton   
h.p = uicontrol('style','pushbutton','units','pixels',...
                'position',ok_pos,'string','OK',...
                'callback',@p_call);
            
%% Focus window + update selection with callback through appdata

% make sure we are on screen
movegui(h.f)
set(h.f, 'Visible','on'); drawnow;

try
    % Give default focus to the listbox *after* the figure is made visible
    uicontrol(h.p);
    c = matlab.ui.internal.dialog.DialogUtils.disableAllWindowsSafely();
    uiwait(h.f);
    delete(c);
catch
    if ishghandle(h.f)
        delete(h.f)
    end
end

if isappdata(0,'CheckboxListDialogAppData__')
    ad = getappdata(0,'CheckboxListDialogAppData__');
    selection = cell2mat(ad.selection).';
    rmappdata(0,'CheckboxListDialogAppData__')
else
    % figure was deleted
    selection = [];
end
drawnow; % Update the view to remove the closed figure (g1031998)

%% Callbacks        
    % Pushbutton callback
    function p_call(varargin)
        ad.selection = get(h.c,'Value');
        setappdata(0,'CheckboxListDialogAppData__',ad);
        close(h.f)
    end
end

