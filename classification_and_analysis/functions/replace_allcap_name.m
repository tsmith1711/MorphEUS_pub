function full_names = replace_allcap_name(abbv_names)
    % converts all drug abbreviations to all cap versions

    abbv_names = upper(abbv_names);
    abbv_names = strrep(abbv_names,'_025x',' L');
    abbv_names = strrep(abbv_names,'_3x', ' H');
    
    full_names = abbv_names;
end 