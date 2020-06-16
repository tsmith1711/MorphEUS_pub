function drug_info_table = get_drug_info(final_data_table, choice_extension, conditions)
%GET_DRUG_INFO Gets a table of information about the drugs in a final_data_table
%   Gets a unique list of choice_extension, and then from that unique list
%   finds all the corresponding other information for each row:
%       - EXP, ID, DOSE, and combinations of those
%
%   Often, choice_extension ends up getting duplicated by one of these
%   combinations (eg DRUG -> ID_DOSE, and DRUG_EXP -> ID_DOSE_EXP),
%   But it is important that we still respect basing our row organization
%   off the unique values of choice_extension so that we get a full list of
%   all combos.

%% Obtain the table

% We want to base the rows of our table off choice extension, so we get all
% combos
drug_info_table = table;
[uout, uindx, ~] = unique(final_data_table.(choice_extension));
drug_info_table.(choice_extension) = uout;

drug_info_table.ID = final_data_table.ID(uindx);

drug_info_table.DOSE = extractAfter(drug_info_table.(choice_extension), '_');
drug_info_table.DOSE(strcmp(drug_info_table.DOSE, '')) = {'none'};

drug_info_table.EXP = final_data_table.EXP(uindx);

if conditions %We need to fix the dose name using the DRUG col.
    drug_info_table.DRUG = final_data_table.DRUG(uindx);
    
    drug_info_table.DOSE = extractAfter(drug_info_table.DRUG, '_');
    drug_info_table.DOSE(strcmp(drug_info_table.DOSE, '')) = {'none'};
    
    drug_info_table.DRUG = []; % We don't actually want the drug col if it wasnt choice_extension, just use ID or ID_DOSE
end

drug_info_table.ID_EXP = strcat(drug_info_table.ID, '_', drug_info_table.EXP);

drug_info_table.ID_DOSE = strcat(drug_info_table.ID, '_', drug_info_table.DOSE);
drug_info_table.ID_DOSE(contains(drug_info_table.ID_DOSE, 'none')) = extractBefore(drug_info_table.ID_DOSE(contains(drug_info_table.ID_DOSE, 'none')), '_');

drug_info_table.ID_DOSE_EXP = strcat(drug_info_table.ID_DOSE, '_', drug_info_table.EXP);

end

