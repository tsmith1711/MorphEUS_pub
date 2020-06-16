function [doses, doses_by_id, drugs_by_id] = get_dose_list(drugs, ids)
%GET_DOSE_LIST Extracts a list of common doses between a list of drugs
%   INPUT PARAMETERS:
%       - drugs: A horiz cell array of drug-dose treatments, eg 'RIF_025x' or 'DMSO'
%       - ids: A horiz cell array of drug treatments without doses, eg 'RIF' or 'DMSO'
%    
%   OUTPUT PARAMETERS:
%       - doses: A horiz cell array of all the unique dosages contained in drugs
%       - doses_by_id: A struct of ID -> drugs under that ID //calculated with substrings; might be better to calculate with actual table
%       - drugs_by_id: A struct of ID -> doses under that ID

%% Create structs to contain list of drugs and doses for each id

% Initialize empty structs
doses_by_id = cell2struct(cell(1, length(ids)), ids, 2);
drugs_by_id = cell2struct(cell(1, length(ids)), ids, 2);

% Fill structs
for i=1:length(ids)
    id = ids{i};
    
    %find all the drugs for each id
    drugs_by_id.(id) = drugs(contains(drugs, id));
    
    %find the doses for each id
    id_with_underscore = strcat(id, '_');
    
    for j=1:length(drugs_by_id.(id))
        drug = drugs_by_id.(id){j};
        dose = replace(drug, id_with_underscore, '');
        
        % If the dose is just the drug name, then there are no doses for
        % this drug; in this case set the dose list to empty
        if strcmp(dose, drug)
            continue
        else
            doses_by_id.(id){j} = dose;
        end
    end
    
end

%% Find a unique list of doses
doses = struct2cell(doses_by_id);
doses = cellfun(@transpose, doses, 'UniformOutput', false);
doses = unique(vertcat(doses{:})).';
    
