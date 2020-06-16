function fig_handle = pca_plot(data_table,pc_labels,idx_type,gcs)
%%% Quick code for creating the 3D pca plot
% Automatically applies graph_color_s
% outputs figure handle
% inputs: data table to use (ex. scores_table)
% drugs: (but could also input, for example, images, dates, etc)
% gcs: if to use graph_color_s

    % default to doing gcs
    if ~exist('gcs','var')
       gcs = true; 
    end
    % default to doing drug indicies 
    if ~exist('idx_type','var')
       idx_type = 'DRUG'; 
    end
    % calculate the numeric columns in the table
    numeric_cols = varfun(@isnumeric,data_table,'OutputFormat', 'uniform');
    
    % obtain numeric information
    numeric_data = data_table{:,numeric_cols};
    
    % get drugs
    drugs = unique(data_table.(idx_type))';
    
    % get distingusihable colors (only needed if not plotting by drug)
    c = distinguishable_colors(size(drugs,2),{'w','k'}); 
    
    % get indexes corresponding to idx_type (default DRUG)
    indexes = cell2struct(cell(1,length(drugs)), drugs, 2);

    for i = drugs
        drug = i{1};
        indexes.(drug) = find(strcmp(data_table.(idx_type), drug)).';
    end
    
    fig_handle = figure();
    for i = 1:length(drugs)
        drug = drugs{i};
        scatter3(numeric_data(indexes.(drug),1), numeric_data(indexes.(drug),2), ...
            numeric_data(indexes.(drug),3),40,c(i,:),'o','filled')
        hold on 
    end

    legend(drugs,'Color',[0.8 0.8 0.8],'Interpreter','none');
    xlabel(pc_labels(1));
    ylabel(pc_labels(2));
    zlabel(pc_labels(3));
    grid on;
    ax = gca;
    ax.Color = [0.8 0.8 0.8];
    ax.GridColor = [1 1 1];
    
    % only apply graph_color_s if plotting by drug 
    if strcmp(idx_type, 'DRUG') && gcs
        try
            graph_color_s(gcf)
        catch
            disp("could not apply graph_color_s")
        end 
    end

    
end

