function bayes_struct = make_bayes_struct(drug_list)
% input some drugs and then get a struct to sture bayes information
% input all unique drugs (maybe as knn drugs?)

    % want to get the correct dimension to use of drug_list
    [~,i] = max(size(drug_list));
    
    if i == 1
        drug_list = drug_list';
    end 
    
    % make an empty struct
    bayes_struct = struct;
    for i = drug_list
        drug = i{1};
        % make a sub structure for that drug
        bayes_struct.(drug) = struct;
        % aaaand loop again!!
        for j = drug_list
            drug2 = j{1};
            
            % this will be overwritten as things get added
            bayes_struct.(drug).(drug2) = 0;
        end 
        
    end 
end

