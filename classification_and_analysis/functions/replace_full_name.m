function full_names = replace_full_name(abbv_names)
    % converts abbreviations into full names for the sake of the plot
    
    abbv_names = strrep(abbv_names, 'Amp','ampicillin');
    abbv_names = strrep(abbv_names, 'Ctax','cefotaxime');
    abbv_names = strrep(abbv_names, 'Mer','meropenem');
    abbv_names = strrep(abbv_names, 'INH','isoniazid');
    abbv_names = strrep(abbv_names, 'EMB','ethambutol');
    abbv_names = strrep(abbv_names, 'ETA','ethionamide');
    abbv_names = strrep(abbv_names, 'IMI','imipenem');
    abbv_names = strrep(abbv_names, 'Van','vancomycin');
    abbv_names = strrep(abbv_names, 'Cyc','cycloserine');
    abbv_names = strrep(abbv_names, 'Carb','carbenicillin');
    abbv_names = strrep(abbv_names, 'Del','delamanid');
    abbv_names = strrep(abbv_names, 'Lev','levofloxacin');
    abbv_names = strrep(abbv_names, 'Mox','moxifloxacin');
    abbv_names = strrep(abbv_names, 'Clz','clofazimine');
    abbv_names = strrep(abbv_names, 'MIT','mitomycin');
    abbv_names = strrep(abbv_names, 'Olf','ofloxacin'); 
    abbv_names = strrep(abbv_names, 'Kan','kanamycin');
    abbv_names = strrep(abbv_names, 'Amk','amikacin');
    abbv_names = strrep(abbv_names, 'Cam','chloramphenicol');
    abbv_names = strrep(abbv_names, 'Cla','clarithromycin');
    abbv_names = strrep(abbv_names, 'Dox','doxycycline'); 
    abbv_names = strrep(abbv_names, 'Gent','gentamicin');
    abbv_names = strrep(abbv_names, 'Strep','streptomycin');
    abbv_names = strrep(abbv_names, 'Tet','tetracycline');
    abbv_names = strrep(abbv_names, 'Tig','tigecycline');
    abbv_names = strrep(abbv_names, 'Lin','linezolid');
    abbv_names = strrep(abbv_names, 'Pre','pretomanid');
    abbv_names = strrep(abbv_names, 'Cer','cerulenin');
    abbv_names = strrep(abbv_names, 'Gra','gramicidin');
    abbv_names = strrep(abbv_names, 'Mon','monensin');
    abbv_names = strrep(abbv_names, 'Nig','nigericin');
    abbv_names = strrep(abbv_names, 'Nis','nisin');
    abbv_names = strrep(abbv_names, 'Thi','thioridazine');
    abbv_names = strrep(abbv_names, 'RifT','rifapentine');
    abbv_names = strrep(abbv_names, 'Sulf','sulfamethizole');
    abbv_names = strrep(abbv_names, 'BDQ','bedaquiline');
    abbv_names = strrep(abbv_names, 'RIF','rifampicin');
    abbv_names = strrep(abbv_names, 'PZA','pyrazinamide');
    
    full_names = abbv_names;
end 