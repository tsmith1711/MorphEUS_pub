function hist_data = hist_file(drugs,numcols,colnames,data)
% outputs data ready for a histogram
% inputs:
% drugs: cell array of drugs
% numcols: logical array of columns with numeric information
% colnames: names of all variables in data table
% data: data table
% example: syto_hist = hist_file(drugs,numeric_syto_cols, fluorescence_column_names,normalized_syto_groups);


% create empty variable to return
hist_data = []; 

% ignore name columns, just want column with data
numcolnames = colnames(numcols);

    for i = drugs % for each drug 
        drug = i{1};
        hist_data.(drug) = [];
        for l = 1:length(numcolnames) % all numeric columns in data
            variable_title = numcolnames{l};
            hist_data.(drug).(variable_title) = [];
            for j = 1:length(data.(drug)) % all photos for given drug
                tab = data.(drug){j};
                tab = tab{:,numcols};
                onecol = tab(:,l)';
                hist_data.(drug).(variable_title) = [hist_data.(drug).(variable_title) onecol];
            end
        end
    end
end

