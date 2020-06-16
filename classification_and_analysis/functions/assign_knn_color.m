function node_colors = assign_knn_color(knn_digraph,color_struct,all_categories)
% will use structures color_struct and all_categories to assign colors for
% the nodes of knn_digraph to color the graph in a more readable manner 
% all of the fields in all_categories must be named the same as the
% structures in color_struct for this to work (or things will just get
% colored as other often)

    node_colors = [];
    
    % get all of the different categories in all_categories
    % we want this and not the categories in color_struct because
    % color_struct is the same regardless if we are using large_group or
    % not whereas all_categories will change 
    category_fields = fieldnames(all_categories)';
    
    for row = 1:height(knn_digraph.Nodes)
        % counted variable keeps track and allows things to be put in the
        % "other" category
        counted = false;
        % get current node name 
        name = table2array(knn_digraph.Nodes(row,1));
        
        % loop through all possible moa category fields 
        for j = category_fields
            moa_category = j{1};
            category_list = all_categories.(moa_category);
            if ismember(name,category_list)
                node_colors = [node_colors; color_struct.(moa_category)];
                % set counted to true makes sure we won't label this as
                % other since it was selected at some point 
                counted = true;
            end
        end 
        % if was not in any of the categoires, mark it as other
        if ~counted
            node_colors = [node_colors; color_struct.other];
        end
    end 
end

