function full_names = replace_allcap_name(abbv_names)
    % converts all drug abbreviations to all cap versions
%     % will also change _025x to " low" and _3x to " high" 
%     
%     abbv_names = strrep(abbv_names, 'Amp','AMP');
%     abbv_names = strrep(abbv_names, 'Ctax','CTAX');
%     abbv_names = strrep(abbv_names, 'Mer','MER');
% %     abbv_names = strrep(abbv_names, 'INH','isoniazid');
% %     abbv_names = strrep(abbv_names, 'EMB','ethambutol');
% %     abbv_names = strrep(abbv_names, 'ETA','ethionamide');
% %     abbv_names = strrep(abbv_names, 'IMI','imipenem');
%     abbv_names = strrep(abbv_names, 'Pip','PIP');
%     abbv_names = strrep(abbv_names, 'Oxa','OXA');
%     abbv_names = strrep(abbv_names, 'PenG','PENG');
%     abbv_names = strrep(abbv_names, 'Van','VAN');
%     abbv_names = strrep(abbv_names, 'Cyc','CYC');
%     abbv_names = strrep(abbv_names, 'Carb','CARB');
%     abbv_names = strrep(abbv_names, 'Del','DEL');
%     abbv_names = strrep(abbv_names, 'Dau','DAU'); 
%     abbv_names = strrep(abbv_names, 'Nal','NAL');
%     abbv_names = strrep(abbv_names, 'Lev','LEV');
%     abbv_names = strrep(abbv_names, 'Mox','MOX');
%     abbv_names = strrep(abbv_names, 'Clz','CLZ');
% %     abbv_names = strrep(abbv_names, 'MIT','MIT');
%     abbv_names = strrep(abbv_names, 'Olf','OLF'); 
%     abbv_names = strrep(abbv_names, 'Kan','KAN');
%     abbv_names = strrep(abbv_names, 'Amk','AMK');
%     abbv_names = strrep(abbv_names, 'Cam','CAM');
% %     abbv_names = strrep(abbv_names, 'AZT','azithromycin');
%     abbv_names = strrep(abbv_names, 'Cla','CLA');
%     abbv_names = strrep(abbv_names, 'Dox','DOX'); 
%     abbv_names = strrep(abbv_names, 'Gent','GENT');
%     abbv_names = strrep(abbv_names, 'Strep','STREP');
%     abbv_names = strrep(abbv_names, 'Tet','TET');
%     abbv_names = strrep(abbv_names, 'Tig','TIG');
%     abbv_names = strrep(abbv_names, 'Lin','LIN');
%     abbv_names = strrep(abbv_names, 'Pre','PRE');
%  %   abbv_names = strrep(abbv_names, 'CCCP','carbonyl cyanide 3-chlorophenylhydrazone');
%     abbv_names = strrep(abbv_names, 'Cer','CER');
%     abbv_names = strrep(abbv_names, 'Gra','GRA');
%     abbv_names = strrep(abbv_names, 'Mon','MON');
%     abbv_names = strrep(abbv_names, 'Nig','NIG');
%     abbv_names = strrep(abbv_names, 'Nis','NIS');
%     abbv_names = strrep(abbv_names, 'Tri','TRI');
%     abbv_names = strrep(abbv_names, 'Ver','VER');
%     abbv_names = strrep(abbv_names, 'Thi','THI');
%     abbv_names = strrep(abbv_names, 'RifT','RIFT');
%     abbv_names = strrep(abbv_names, 'Nit','NIT');
%     abbv_names = strrep(abbv_names, 'Sulf','SULF');
%     abbv_names = strrep(abbv_names, 'BDQ','bedaquiline');
%     abbv_names = strrep(abbv_names, 'RIF','rifampicin');
%     abbv_names = strrep(abbv_names, 'PZA','pyrazinamide');

    abbv_names = upper(abbv_names);
    abbv_names = strrep(abbv_names,'_025x',' L');
    abbv_names = strrep(abbv_names,'_3x', ' H');
    
    full_names = abbv_names;
end 