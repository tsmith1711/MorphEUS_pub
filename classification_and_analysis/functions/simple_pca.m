function hpca = simple_pca(data_table,idx_type)
% input a normalized table and we will pca the h*ck out of it

    % default to doing drug indicies 
    if ~exist('idx_type','var')
       idx_type = 'DRUG'; 
    end

    numeric_final_data_cols = varfun(@isnumeric,data_table,'OutputFormat', 'uniform');
    % quick and dirty pca just input a table
    ndata = data_table{:,numeric_final_data_cols}; %numeric data from normalized table 

    % non numeric data
    non_numeric = data_table(:,~numeric_final_data_cols);

    % Columns of the scores variable corresponds to the principal components (from first to last, decreasing order)
    % dataInPrincipalComponentSpace = ndata*coeff
    [~,scores,pcvars] = pca(ndata);
    pcvari = pcvars/sum(pcvars);

    % put data in a table for better organization 
    scores_table = horzcat(non_numeric,array2table(scores));

    % find all principal components and percentages for axes labels 
    pc_labels = cell(1,length(pcvari));
    for i = 1:length(pcvari)
        pc_labels(i) = cellstr(strcat("PC:",num2str(round(i,2)), " ", (num2str(pcvari(round(i,2))*100)), "%"));
    end
    
    hpca = pca_plot(scores_table,pc_labels,idx_type,false);
end

