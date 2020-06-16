function connectivity_plot(results_table,color_struct,all_categories,map,wbg)
% create connectivity plot from bayes_struct 
% needs color_struct and all_categories (probably created in the knn color
% setup section) in order to create node colors
% needs helper function node_colors
% map -> a colormap, default to a cool red one
% wbg -> boolean on whether to create a white background, default false

%% set default vars
    % by default, do a grey background 
    if ~exist('wbg','var')
        wbg = false;
    end 
    % create a default colormap
    if ~exist('map','var')
        vec = [100; 50; 0;];
        raw = [ 1 1 1; ...
            255/255 104/255 104/255; ...
            47/255 0 0];
        N = 128;
        map = interp1(vec,raw,linspace(100,0,N),'pchip');
    end 
    
%% set up for figure 
% just store anything in S you wanna use in the callback
    % Plot different plots according to slider location.
    % some settings
    fig_w = 1000;
    fig_h = 650;

    % figure handle 
    S.fh = figure('units','pixels',...
                  'position',[300 600 fig_w fig_h],...
                  'name','connectivity plot',...
                  'numbertitle','off');    
              
%     S.ax = axes('unit','pix',...
%                 'position',[20 80 560 510]);

    %% begin connectivity  plot 
    % create digraph 
    drug_names = results_table.DRUG;
    neighbor_names = results_table.NEIGHBOR;
    
    knn_digraph =  digraph(drug_names,neighbor_names, ...
        results_table.DISTANCE);
    % get node colors 
    node_colors = assign_knn_color(knn_digraph,color_struct,all_categories);
    % create knn plot 
    subplot(1,3,1:2) % allow more space for connectivity plot than histogram
    % plot 
    S.LN = plot(knn_digraph,'Layout','force',...
                            'LineWidth',2.5,...
                            'NodeColor',node_colors,...
                            'ShowArrows','off');  

    % if not white background, set a grey background 
    if ~wbg
        set(gca,'Color',[.8 .8 .8]);
    end 
    % save white background variable in callback struct 
    S.wbg = wbg;

    % set up colorbar                     
    colorbar
    colormap(S.fh, map);
    % get min and max distances for colorbar limits
    max_dist = max(results_table.DISTANCE);
    min_dist = min(results_table.DISTANCE);
    % set colorbar limits 
    caxis([min_dist max_dist])
    
    % add colors to nodes 
    S.LN.EdgeCData= knn_digraph.Edges.Weight;
    % make it pretty
    S.LN.MarkerSize = 5;
    S.LN.NodeFontSize = 10;
    % store info in S for later
    S.results_table = results_table;
    S.color_struct = color_struct;
    S.all_categories = all_categories;
    S.map = map;
    
    % set graph title
    total_drugs = length(unique(results_table.DRUG));
    graph_title = strcat("Displaying ",num2str(total_drugs)," of ", num2str(total_drugs)," total drugs");
    
    % get number of connections for subtitle
    total_connections = length(results_table.DISTANCE);
    graph_subtitle = strcat(num2str(total_connections), " total connections");
    % set title 
    title({graph_title,graph_subtitle})
    
    
    %% create gui text                     
    % set text                    
    bgcolor = S.fh.Color;
    S.bl1 = uicontrol('Parent',S.fh,'Style','text','Position',[5,30,23,23],...
                'String','0','BackgroundColor',bgcolor);
    S.bl2 = uicontrol('Parent',S.fh,'Style','text','Position',[fig_w-50,30,23,23],...
                'String',num2str(round(max_dist)),'BackgroundColor',bgcolor);
    S.bl3 = uicontrol('Parent',S.fh,'Style','text','Position',[fig_w/2,10,100,23],...
                'String','Cutoff %','BackgroundColor',bgcolor);
    
    
    %% create histogram 
    % make histogram plot
    subplot(1,3,3)
    S.hist = histogram(results_table.DISTANCE,[0 10 20 30 40 50 60 70 80 90 100]);
    title("histogram")
    ylabel("number of connections")

    %% create gui slider 
    % make slider and assign callback function (this needs to go last)
    S.sl = uicontrol('style','slide',...
                     'unit','pix',...
                     'position',[30 20 fig_w-90 30],...
                     'min',1,'max',max_dist,'val',1,...
                     'sliderstep',[1/20 1/20],...
                     'callback',{@slider_callback,S}); 
                 
    
            
end 
             

function slider_callback(varargin)
% callback function for the connectivity plot 
    % h = slider
    % S = struct with all the good stuff
    [h,S] = varargin{[1,3]};
    
    %% connectivity plot 
    subplot(1,3,1:2)
    % pull out that juicy information we stored in S
    results_table = S.results_table;
    color_struct = S.color_struct;
    all_categories = S.all_categories; 
    map = S.map;

    %% calculate new connectivity plot with cutoff 
    % find new cutoff information based on slider bar
    new_cutoff = h.Value;
    above_inds =  find(results_table.DISTANCE > new_cutoff);
    % keep rows that are above the cutoff 
    successful_rows = results_table(above_inds,:);
    
    all_possible_drugs = unique([unique(successful_rows.DRUG)' unique(successful_rows.NEIGHBOR)']);

    % let's get the number of drugs that are showing up 
    numb_drugs = length(all_possible_drugs);
    total_drugs = length(unique(results_table.DRUG));
    % get min and max distances for colorbar
    min_dist = min(successful_rows.DISTANCE);
    max_dist = max(successful_rows.DISTANCE);
    
    % if the drug list for successful_rows and results_table is not the
    % same
    if ~isequal(all_possible_drugs,unique(results_table.DRUG))
        % this means some got left out
        left_out = setdiff(unique(results_table.DRUG),all_possible_drugs)';    
        % create a new row for that drug
        for q = left_out
            missing_drug = q(1);
            % accomdate different versions of results_table
            try
                new_row = table(q,q,0,q,'VariableNames',{'DRUG','NEIGHBOR','DISTANCE','NEIGHBOR_MOA'});
                successful_rows = vertcat(successful_rows,new_row);
            % accomdate different versions of results_table
            catch
                new_row = table(q,q,0,'VariableNames',{'DRUG','NEIGHBOR','DISTANCE'});
                successful_rows = vertcat(successful_rows,new_row);
            end 
        end 
    end 

     
    % create a new digraph with the successful rows
    drug_names = successful_rows.DRUG;
    neighbor_names = successful_rows.NEIGHBOR;
    
    knn_digraph =  digraph(drug_names,neighbor_names, ...
    successful_rows.DISTANCE);
    % generate node colors 
    node_colors = assign_knn_color(knn_digraph,color_struct,all_categories);
    
    % set string to show cutoff information 
    S.bl3.String = strcat("Cutoff % = ", num2str(new_cutoff));
    % re plot 
    S.LN = plot(knn_digraph,'Layout','force','LineWidth',2.5,'NodeColor',node_colors,...
        'ShowArrows','off');
    % add weights 
    S.LN.EdgeCData= knn_digraph.Edges.Weight;
    
    % make it pretty
    S.LN.MarkerSize = 5;
    S.LN.NodeFontSize = 10;
    % set background color 
    if ~S.wbg
        set(gca,'Color',[.8 .8 .8]);
    end

    % set map color
    colormap(S.fh, map);
    colorbar

    % set limits of colorbar 
    caxis([min_dist max_dist])
    
    %% get title information 
    
    graph_title = strcat("Displaying ",num2str(numb_drugs)," of ", num2str(total_drugs)," total drugs");
    
    % get number of connections for subtitle
    total_connections = length(successful_rows.DISTANCE);
    graph_subtitle = strcat(num2str(total_connections), " total connections");
    title({graph_title,graph_subtitle})
    %% histogram plot 
    % now do the histogram
    subplot(1,3,3)
    S.hist = histogram(successful_rows.DISTANCE,[0 10 20 30 40 50 60 70 80 90 100]);
    title("histogram")
    ylabel("number of connections")
end 