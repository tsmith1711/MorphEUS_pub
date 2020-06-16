function merged_set = merge_reps(filtered_set)
% INPUT: filtered_bacteria, filtered_FM_fluoresence,
% filtered_syto_fluoresence
% OUTPUT: merge those images baby 


    % get a list of all the drugs in the structure
    druglist = fieldnames(filtered_set)';

    % new structure to store everything
    merged_set = struct;

    % for each field(drug)
    for i = druglist
        drug = i{1};

        imgs = fieldnames(filtered_set.(drug))';

        % get the set of images for each drug 
        % want to concatenate tables that have the same pad and spot and
        % different img_replicates 
        % images are of format:
        % mm_04252019_ph_techrep1_e_02__9_R3D_BaSiC
        %or
        % mm_04102019_timecourse_set1_rep3_d_02__10_R3D_BaSiC

        % new structure to re-sort
        for j = imgs
            img = j{1};
            [start_indx, end_indx] = regexp(img, '[a-e]_0\d__(1[0-2]|[0-9])');
            
            if ~isempty(start_indx)
                % gets d_02__10
                classifier = img(start_indx:end_indx);

                % gets d
                pad = classifier(1);
                % gets 10
                spot = classifier(7:end);

                % get suffix
                suffix = img(end_indx+2:end);

                % gets 02
                %img_replicate = extractBefore(extractAfter(classifier,"_"),"_");

                % use start_indx to help get set-rep -- will work for both techrep
                % gets mm_04102019_timecourse_set1_rep3
                img_before_classifier = img(1:start_indx-2);

                % combine pad and spot
                pad_spot = strcat(pad,"_",spot);

                % get final titl
                header_pad_spot = strcat(img_before_classifier,"_",pad_spot,"_",suffix);

                % remove row names to prevent conflict;
                filtered_set.(drug).(img).Properties.RowNames = {};

                %now we need to check if this already exists -- if so, we want to
                %concateenate this table with the existing one 
                if isfield(merged_set,drug)
                    % if the field exists already, concatenate with new data
                    if isfield(merged_set.(drug),header_pad_spot)
                        merged_set.(drug).(header_pad_spot) = ...
                            vertcat(merged_set.(drug).(header_pad_spot), filtered_set.(drug).(img));
                    else
                        % if field doesn't exist, create it 
                        merged_set.(drug).(header_pad_spot) = filtered_set.(drug).(img);
                    end
                else 
                    % if field doesn't exist, create it
                    merged_set.(drug).(header_pad_spot) = filtered_set.(drug).(img);
                end 
            
                % if there aren't multiple replicate imagese, just go and
                % store it as the image
            else
                merged_set.(drug).(img) = filtered_set.(drug).(img);
            end
        end

    end 
end

