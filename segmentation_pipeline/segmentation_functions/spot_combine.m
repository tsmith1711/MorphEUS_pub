function filtered_struct = spot_combine(filtered_struct,info_table)
% INPUT: filtered_bacteria, filtered_FM_fluoresence,
% filtered_syto_fluoresence
% info_table: a previous readtable with pad info
% OUTPUT: merge spots given by info_table


    info_array = table2array(info_table);

    % get all unique values
    unq_drugs = unique(info_array)';
    
    
% now want to loop through these unique values 
for i = unq_drugs
    drug = i{1};
    loc_inds = contains(info_array,drug);
    
    total_spots = sum(sum(loc_inds));
    
    % only want to proceed if total_spots > 1 because that means we have
    % two locations to merge 
    %%
    if total_spots > 1
        % find instance of last underscore
        last_und_pos = find(drug == '_', 1, 'last');
        % get the actual name of the drug by getting rid of the number
        % underscore at the end
        drug_name =drug(1:last_und_pos-1);
        
        % get all image names 
        all_imgs = fieldnames(filtered_struct.(drug_name));
        
        % find the location in the array for each
        locations = find(loc_inds)';
        
        
        % create array to store information 
        all_img_info = struct;
    
        % go through each location
        for location = locations
            % convert to subscript indexes
            [row,col] = ind2sub(size(info_array),location);
            % get location in table
            table_drug = info_table(row,col);
            
            % get the spot number (column)
            spot_num = table_drug.Properties.VariableNames{1};
            % get rid of the default x added by MATLAB and convert to
            % number 
            spot_num = extractAfter(spot_num,'x');
            
            % get the spot pad (row) (letter)
            spot_pad = table_drug.Properties.RowNames{1};
            
            % now what helP
            % we have the pad and number info so that means we could get
            % the right image
            
            % get the pad information as it appears in the image
            spot = strcat(spot_pad,"_",spot_num);
            
            % want to see if this spot information appears in one of the
            % images for the drug 
            if any(contains(all_imgs,spot))
                img_loc = contains(all_imgs,spot);
                % get the image corresponding to the spot
                img = all_imgs(img_loc);
                img = img{1};
                % save the information from this image in all_img_info
                all_img_info.(img) = filtered_struct.(drug_name).(img);
                
                % remove the field with this image information from the
                % filtered_struct
                filtered_struct.(drug_name) = rmfield(filtered_struct.(drug_name),img);
            end
            
           
        end % end going through each location 
        
        new_concat_imgs = struct2cell(all_img_info);
        combined_imgs = vertcat(new_concat_imgs{:});
        
        img_merged = strcat(img,'_spot_combined');
        % use last image as placeholder new name I guess
        filtered_struct.(drug_name).(img_merged) = combined_imgs;
        
        
    end  % end if total_spots >1
    
    


end

