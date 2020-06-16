function combined_file = csv_combine(data_path,file_keyname,kept_cols)
% combine multiple csv files
% to be used near the beginning of BCP microbej processing: will go through
% the data folder (indicated by data_path), and find all csv files with the
% given keyword (eg all_bacteria.csv, all_bacteria_2.csv, and
% all_bacteria_3.csv will all be called if the keyword is set to
% 'bacteria'). The script will open those files and then concatenate the
% tables contained 
    % file_keyname: whatever the keyword in the file is, ex
    % 'bacteria','feature1','feature2','transverse'
    
    csv_list = dir(strcat(data_path, '/*.csv'));
    csv_list = {csv_list.name};

    % get the indicies of the csv files with keyword
    file_inds = contains(csv_list,file_keyname);
    % get the names of the  files with the keyword
    file_list = csv_list(file_inds);

    % make an empty struct to store results
    file_struct = struct;
    
    % loop through each file with the keyword 
    for i = file_list
        current_file = i{1};

        % load the csv file 
        single_csv = readtable(strcat(data_path,current_file), 'ReadVariableNames',true,'TreatAsEmpty',{'Infinity','NA'});

        % only keep selected columns (if desired)
        if exist('kept_cols','var')
            single_csv = single_csv(:,kept_cols);
        end 
        
        % get the file name without the csv extensiojn
        current_file_no_csv = extractBefore(current_file,'.csv');

        % save each file as a field in the structure
        file_struct.(current_file_no_csv) = single_csv;
    end 

    % get together all of the csv files
    all_csvs = struct2cell(file_struct);

    % vertically concatenate all of the different csvs
    combined_file = vertcat(all_csvs{:});

    
end

