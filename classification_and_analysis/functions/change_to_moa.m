function drug_moas = change_to_moa(drugList,all_categories)
% INPUT list of drugs 
% will read them in and change them to their MoA
% is this like assign MoA? probably
% let me live my life
% all_categories is a struct where ,for example, all_categories.protein has
% an array of all of the different protein drugs

    drug_moas = drugList;
    
    category_fields = fieldnames(all_categories)';
    % go through and find moa of each drug and neighbor in order to create
    % confusion matrix (don't neeed MoA of others and unknowns)
    for k = 1:length(drugList)
        % get drug
        drug1 = drugList(k);
        drug1 = drug1{1};
        % counted variable keeps track and allows things to be put in the
        % "other" category
        counted= false;
        % loop through each moa category in categories
        for j = category_fields
            moa_category = j{1};
            category_list = all_categories.(moa_category);
            % now loop through each drug

            if ismember(drug1,category_list)
                % go through all their drug indicies
                drug_moas(k) = {moa_category};
                % set counted to true makes sure we won't label this as
                % other since it was selected at some point 
                counted = true;
            end
        end 
        % if was not in any of the categorires, mark it as other 
        if ~counted
           drug_moas(k) = {'other'}; 
        end
    end

end

