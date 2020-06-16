function out_table = fix_10drug_ids(data_table)
%FIX_10DRUG_IDS Sets the ID column to drug name without dose
%   Detailed explanation goes here

drug_col = data_table.DRUG;

% Extract drug name before underscore
ids = cellfun(@(d) extractBefore(d, '_'), drug_col, 'UniformOutput', false);

% Fix non-dose drugs
non_dose_drugs = {'INH_control', 'Untreated', 'DMSO'};
for i=non_dose_drugs
    drug = i{1};
    ids(strcmp(drug_col, drug)) = {drug};
end

% if ID col already exists
if any(strcmp(data_table.Properties.VariableNames, 'ID'))
    data_table.ID = ids;
    out_table = data_table;
else
    %if no ID col, insert right after drug
    drug_col_indx = find(strcmp(data_table.Properties.VariableNames, 'DRUG'));
    before_drug_table = data_table(:, 1:drug_col_indx);
    after_drug_table = data_table(:, (drug_col_indx+1):end);
    
    id_col = table();
    id_col.ID = ids;
    
    out_table = [before_drug_table id_col after_drug_table];
end

end

