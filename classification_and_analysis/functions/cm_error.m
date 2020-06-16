function pct_correct = cm_error(cm)
% gets error from confusion chart object
%   input 1: cm - conusionchart object
%   output 1: pct_correct = percent correct of confusion matrix 

    % get success percent
    % get confusion matrixfrom the confusion chart object
    cmat = cm.NormalizedValues;
    % find the sum of the rows to get true classes
    true_class = sum(cmat,2);
    
    % get class labels in graph
    labels = cm.ClassLabels;
    
    % determine if there is an others category
    contains_other = contains('other',labels);
    
    % calcualte cm error 
    if contains_other
        % find the index of the 'other' category
        other_ind = contains(cm.ClassLabels,'other');
        
        % get all of the correct 'other' category
        %correct_others = cmat(end,end);
        correct_others = cmat(other_ind,other_ind);
        
        % add up along diagonal to get all correct
        all_correct = trace(cmat);
        % good correct does not includee correct from others category
        good_correct = all_correct - correct_others;
        % get all of the drugs not including the 'other' category
        
        % get sum of others
        all_others = true_class(other_ind);
        
        % get sum of all not including others
        all_with_no_others = sum(true_class)-all_others;
        pct_correct = round((good_correct/all_with_no_others)*100,2); 
        
    else
        % no others category, so more straightforward
        % add up along diagonal to get all correct
        all_correct = trace(cmat);

        all_drugs = sum(true_class);
        
        pct_correct = round((all_correct/all_drugs)*100,2);
        
    end

end

