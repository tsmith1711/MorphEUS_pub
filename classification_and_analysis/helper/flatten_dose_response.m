%% Flattens multiple doses into columns
%
%
%

%% Find info about chosen_doses
num_chosen_doses = length(chosen_doses);

%% %% Add ID_EXP and DOSE column

%~~ Find numeric and non-numeric parts of table ~~ %
%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');
% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);
% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols);

%~~ Create new column ~~%
%id_exp
id_exp = cellstr(strcat(final_data_table.ID,"_",final_data_table.EXP));
non_numeric.ID_EXP = id_exp;

%dose
dose_col = cellstr(extractAfter(final_data_table.DRUG, strcat(final_data_table.ID,"_")));
non_numeric.DOSE = dose_col;

% reorganize table with new column
final_data_table = [non_numeric numeric_data];


%% Recalculate numeric columns and data

%find numeric columns of final table
numeric_final_data_cols = varfun(@isnumeric,final_data_table,'OutputFormat', 'uniform');
%find column names
numeric_final_col_names = final_data_table.Properties.VariableNames(numeric_final_data_cols);
% numeric data
numeric_data = final_data_table(:,numeric_final_data_cols);
% Get table of non numeric values
non_numeric = final_data_table(:,~numeric_final_data_cols);

%% Split tables into tree
%Looks somethihg like this:
%
%               rep1   dose1
%          id1- rep2 - dose2
%   exp1 -      rep3   dose3
% -        id2
%   exp2
%
%

exp_list = unique(final_data_table.EXP).';

split_tables = struct();

% Split by EXP
for i=exp_list
    exp = i{1};
    
    %Get subtable for this experiment
    exp_table = final_data_table(strcmp(final_data_table.EXP, exp), :);

    %Branch struct
    split_tables.(exp) = struct();
    
    % Split by ID
    for j=chosen_ids
        id = j{1};

        %Get subtable for this id
        id_table = exp_table(strcmp(exp_table.ID, id), :);

        %Branch struct
        split_tables.(exp).(id) = struct();
        
        % Split by REP
        rep_list = unique(id_table.REP).';
        for k=rep_list
            rep = k{1};
            
            %Get subtable for this rep
            rep_table = id_table(strcmp(id_table.REP, rep), :);

            %Branch struct
            split_tables.(exp).(id).(rep) = struct();
            
            % Split by drug (drug_dose)
            for q=chosen_drugs_by_id.(id)
                drug = q{1};
                
                %Get subtable for this dose
                drug_table = rep_table(strcmp(rep_table.DRUG, drug), :);

                % Populate drug_dose field with table
                split_tables.(exp).(id).(rep).(drug) = drug_table;
            end
        end

    end
end

% %% Ensure rep list is valid
% if sum(strcmp({'rep1', 'rep2', 'rep3'}, rep_list)) ~= length(rep_list)
%     error("Rep list invalid; fix the code you dummy")
% end

%% for Untreated, we can just put all of the untreated into one rep and call it a day
% this is just a fix for that bayes stuff right now le'ts go
if bayesian
    unt_reps = fieldnames(split_tables.drug34_doseResponse.Untreated)';

    % add a new field for all reps

    split_tables.drug34_doseResponse.Untreated.all_reps = struct;

    % want to make this flexible later but now just need to get things going 
    for i = unt_reps
        rep=i{1};
        unt3x = split_tables.drug34_doseResponse.Untreated.(rep).Untreated_3x;
        unt025x = split_tables.drug34_doseResponse.Untreated.(rep).Untreated_025x;

        if strcmp(rep,"rep1")
            split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_3x = unt3x;
            split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_025x = unt025x;
        else
            split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_3x = ...
                vertcat(split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_3x, unt3x);
            split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_025x = ...
                vertcat(split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_025x, unt025x);

        end
        split_tables.drug34_doseResponse.Untreated = ...
            rmfield(split_tables.drug34_doseResponse.Untreated,rep);
    end 

    % go replace the rep column to all have the same thing so it works on next
    % step

    maxlen = length(split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_025x.REP);

    for i = 1:maxlen
        split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_025x.REP(i) = ...
            {'ALL'};
        split_tables.drug34_doseResponse.Untreated.all_reps.Untreated_3x.REP(i) = ...
            {'ALL'};


    end 

end 
%% Fill subtables for drugs that are missing any doses

%loop through exp
for i=exp_list
    exp = i{1};
    
    %loop through id
    for j=chosen_ids
        id = j{1};
        rep_list = fieldnames(split_tables.(exp).(id)).';
        
        %loop through rep
        for k=rep_list
            rep=k{1};
            
            % If there are fewer doses for this ID than for other ones
            if length(chosen_drugs_by_id.(id)) < num_chosen_doses

                % Display to console that we're doing this
                disp(strcat("Filling in horizontal data for ",...
                    id, " at [", exp, ", ", rep, "] (Only have ",...
                    num2str(length(chosen_drugs_by_id.(id))), " dose(s), but need ",...
                    num2str(num_chosen_doses), ")"))

                % Copy the table from the first dose recorded for this drug
                %(99% of the time there will only 1; eg Untreated / DMSO)
                drug_to_copy = chosen_drugs_by_id.(id){1}; 
                table_to_copy = split_tables.(exp).(id).(rep).(drug_to_copy);

                %~~ Insert all fields needed with generated names and copy of same table
                split_tables.(exp).(id).(rep) = struct();
                % Generate names (eg DMSO1, DMSO2)
                names = cellfun(@(x) strcat(id, num2str(x)), num2cell(1:num_chosen_doses), 'UniformOutput', false);
                for w=names
                    name = w{1};
                    split_tables.(exp).(id).(rep).(name) = table_to_copy;
                end
                
                % Fix chosen_doses_by_drug
                chosen_doses_by_drug.(id) = cell2struct(chosen_doses, names, 2);
            end
        end
    end
end



%% Horizontally concatenate doses for each rep and fix table columns
missing_data = "";
combined_tables = struct();
final_combined_table = table();

% set the rng for reproducibility
rng(rng_val);

%loop through exp
for i=exp_list
    exp = i{1};
    
    combined_tables.(exp) = struct();
    
    %loop through id
    for j=chosen_ids
        id = j{1};
        
        combined_tables.(exp).(id) = struct();
        
        rep_list = fieldnames(split_tables.(exp).(id)).';
        %loop through rep
        for k=rep_list
            rep = k{1};
            %disp(strcat("Combining doses for: [", exp, ", ", id, ", ", rep, "]")) 
            
            fields = fieldnames(split_tables.(exp).(id).(rep));
            
            % Subset each dose table to numeric cols and non numeric cols
            non_numeric_tables = structfun(@(t) t(:, ~numeric_final_data_cols), split_tables.(exp).(id).(rep), 'UniformOutput', false);
            numeric_tables = structfun(@(t) t(:, numeric_final_data_cols), split_tables.(exp).(id).(rep), 'UniformOutput', false);

            % Append dose to each numeric column name
            for q=fieldnames(numeric_tables).'
                field = q{1};

                dose = chosen_doses_by_drug.(id).(field);

                numeric_tables.(field).Properties.VariableNames = cellfun(@(c) strcat(c, '_', dose), ...
                    numeric_tables.(field).Properties.VariableNames, 'UniformOutput', false);
            end

            %~~ HorzCat tables ~~%
            
            %find tables with extra rows
            table_rows = structfun(@(t) size(t, 1), split_tables.(exp).(id).(rep), 'UniformOutput', false);
            table_rows_mat = cell2mat(struct2cell(table_rows));
            min_num_rows = min(table_rows_mat);
            
            num_extra_rows = structfun(@(n) n - min_num_rows, table_rows, 'UniformOutput', false);
            indx = structfun(@(n) n > 0, num_extra_rows, 'UniformOutput', true).';
            drugs_with_extra_rows = fields(indx).';
            
            % Remove extra rows
            if ~isempty(drugs_with_extra_rows)
                for q=drugs_with_extra_rows
                    drug = q{1};

                    % Generate random row indexes to remove
                    indexes_to_remove = randperm(table_rows.(drug), num_extra_rows.(drug));

                    % Find img col for display
                    local_imgs = non_numeric_tables.(drug).IMG;

                    % Remove from both numeric and non_numeric tables
                    numeric_tables.(drug)(indexes_to_remove, :) = [];
                    non_numeric_tables.(drug)(indexes_to_remove, :) = [];

                    % Tell user what we have done
                    disp(strcat(drug, " at [", exp, ", ", rep, "] had ", num2str(num_extra_rows.(drug)), ...
                        " too many rows. Removing images (randomly picked): "))
                    removed_local_imgs = local_imgs(indexes_to_remove).';
                    for w = removed_local_imgs
                        disp(strcat(" - ", w{1}))
                    end
                end
            end
            
            % Loop through rows that will be combined, and aggregate
            % non_numeric data from each field in the struct
            combined_non_numeric_header = table();
            for r=1:min_num_rows
                % Get numeric cols - rows=rth row one from each dose
                rth_rows = cellfun(@(t) t(r, :), struct2cell(non_numeric_tables), 'UniformOutput', false);
                original_header = vertcat(rth_rows{:});
                
                
                % Calculate extra cols we want
                
                %imgs (IMG_025x, IMG_3x ...)
                img_cols = table();
                for g=fieldnames(non_numeric_tables).'
                    drug = g{1};
                    img_cols.(strcat('IMG_', chosen_doses_by_drug.(id).(drug))) = convertCharsToStrings(non_numeric_tables.(drug).IMG{r});
                end
                img_cols.IMG = strjoin(img_cols{1,:}, '-@-');
                
                % Date 
                date_col = table();
                if all(strcmp(original_header.DATE{1}, original_header.DATE))
                    date_label = original_header.DATE{1};
                else
                    date_label = strjoin(original_header.DATE, '_&_');
                end
                date_col.DATE = {convertStringsToChars(date_label)};
                
                % dose (DOSE = 025x_&_3x)
                dose_col = table();
                dose_label = strjoin(original_header.DOSE, '_&_');
                if isempty(chosen_doses_by_id.(id))
                    dose_label = 'N/A';
                end
                dose_col.DOSE = {dose_label};
                
                % Drug / ID
                % Set DRUG to ID, so that we don't have to change the
                % pipeline
                original_header.DRUG = original_header.ID;
                
                
                %Remove the cols that we are adding
                cols_to_remove = {'IMG', 'DATE','DOSE', 'BATCH', 'DRUG_DATE', 'DRUG_EXP', 'SET', 'SET_REP'};
                % Only remove if they exist
                for c = cols_to_remove
                    col = c{1};
                    if any(strcmp(col, original_header.Properties.VariableNames))
                        original_header.(col) = [];
                        %disp(strcat("removed ", col));
                    end
                end
                
        %%% DYNAMICALLY REMOVES NON-UNIQUE VARS; DOESN'T WORK WITH DMSO %%%
                % Remove any cols with more than one unique value %%
%                 for h=original_header.Properties.VariableNames
%                     varname = h{1};
%                     if length(unique(original_header.(varname))) > 1
%                         original_header.(varname) = [];
%                         disp(strcat("Removed var ", varname, " because it is not unique across doses."))
%                     end
%                 end

                % Get unique vals of remaining cols
                try
                    original_header = varfun(@unique, original_header);
                catch
                    error("More columns than expected; remove them at around line 255 in flatten_dose_response.m");
                end
                original_header.Properties.VariableNames = extractAfter(original_header.Properties.VariableNames, 'unique_');
                
                % Get summarized row
                new_row = [original_header date_col img_cols dose_col];
                
                
                % Convert all values to string objects
                new_row = varfun(@convertCharsToStrings, new_row);
                new_row.Properties.VariableNames = extractAfter(new_row.Properties.VariableNames, 'convertCharsToStrings_');
                
                % Add to combined_non_numeric_table
                combined_non_numeric_header = [combined_non_numeric_header; new_row];
            end
                
            % Combine numeric tables
            numeric_tables_cell = struct2cell(numeric_tables);
            combined_numeric_table = horzcat(numeric_tables_cell{:});
            
            % Add non_numeric_cols
            combined_table = [combined_non_numeric_header combined_numeric_table];
            
            % Check if table has no rows, and add to combined table if all
            % good.
            if size(combined_table, 1) <= 0
                msg = strcat("Missing data for [", id, ", ", rep, ", ", exp, "]");
                warning(msg);
                missing_data = strcat(missing_data, " | ", msg);
            else
                % Add to aggregate tables
                combined_tables.(exp).(id).(rep) = combined_table;
                final_combined_table = [final_combined_table; combined_table];
            end
        end
    end
end
    
    
    