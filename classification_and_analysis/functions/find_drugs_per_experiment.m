function [drugs_by_experiment, universal_drugs] = find_drugs_per_experiment(data_table)
%FIND_DRUGS_PER_EXPERIMENT Finds the drugs for each experiment in a
%final_data_table
%
%   USED IN conditions pipeline
%
%   Returns: drugs_by_experiment, a struct of experiment names -> drug lists
%            universal_drugs, drugs that are contained in all experiments

if ~any(strcmp(data_table.Properties.VariableNames, 'EXP')) || ~any(strcmp(data_table.Properties.VariableNames, 'DRUG'))
    error("Data tablepassed to find_drugs_by_experiments does not contain either field 'EXP' or 'DRUG'")
end

drugs_by_experiment = struct();

experiments = unique(data_table.EXP).';
for i=experiments
    exp = i{1};
    
    drugs = unique(data_table.DRUG(strcmp(data_table.EXP, exp))).';
    
    drugs_by_experiment.(exp) = drugs;
end

drugs_by_experiment_cell = struct2cell(drugs_by_experiment);

universal_drugs = drugs_by_experiment_cell{1};

if length(experiments) > 1
    indxs = cellfun(@(d) ismember(universal_drugs, d), drugs_by_experiment_cell, 'UniformOutput', false);
    indxs = vertcat(indxs{:});
    indx = sum(indxs) == length(experiments);
    
    universal_drugs = universal_drugs(indx);
end
end

