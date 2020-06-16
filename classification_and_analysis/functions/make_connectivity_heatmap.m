function fig_handle = make_connectivity_heatmap(heatmap_table,tpose)
    % input heatmap table, output figure handle
    % optional input of transpose
    
    if ~exist('tpose','var')
        tpose = false; 
    end 
    
    fig_handle = figure;
    % if transpose
    if tpose
        imagesc(heatmap_table{:,:}')
    else
        imagesc(heatmap_table{:,:})
    end 
    
    % if transpose, need to switch x and y labels
    if tpose
        x_labels = heatmap_table.Properties.RowNames;
        y_labels = heatmap_table.Properties.VariableNames;        
    else 
        x_labels = heatmap_table.Properties.VariableNames;
        y_labels = heatmap_table.Properties.RowNames;
    end 
    set(gca, 'XTick', 1:length(x_labels))
    set(gca, 'YTick', 1:length(y_labels));
    set(gca, 'XTickLabel', x_labels);
    set(gca, 'YTickLabel', y_labels);
    colorbar
    caxis([0 100])
end

