function new_table = multicompare_table(c,nms)
% input matrix c of the pairwise comparison results from a multiple comparison test 
% and nms -- group names from pairwise comparison results
% converts c into a more readable format using the names in nms
    first_col = c(:,1);
    second_col = c(:,2);
    p_vals = c(:,6);
    
    
    % make a new column that replaces the group numbers with group names
    new_first_col = cell(size(first_col));

    % for each group
    for i = 1:length(nms)
        % replace group number with group name in nms
        new_first_col(first_col==i) = nms(i);
    end
    
    % do the same for the second column    
    new_second_col = cell(size(second_col));

    % for each group
    for i = 1:length(nms)
        % replace group number with group name in nms
        new_second_col(second_col==i) = nms(i);
    end
    
    % now want to go through p vals and bring 
    sig_col = cell(size(p_vals));
    
    sig_col(p_vals < 0.05) = {'SIGNIFICANT'};

    new_table = table(new_first_col,new_second_col,p_vals,sig_col,...
        'VariableNames', {'GROUP1','GROUP2','P_VAL','SIGNIFICANT'});
    
    
end

