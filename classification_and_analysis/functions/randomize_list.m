function randomized_drugs = randomize_list(drug_list)
% Input a list of drugs and return that list in a randomized order 
    % copy list of knn_drusg 
    randomized_drugs = drug_list;

    rand_order = randperm(length(drug_list));

    % go through each drug and reorder in the random order defined 
    for p = 1:length(drug_list)
        randomized_drugs(p) = drug_list(rand_order(p));
    end 

end

