%% knn color set up
% creates the colormap, all_categories, and color_struct
% need to have set large_group and efflux_in_lipid booleans

%% Set up workspace 
% set up colormap
vec = [100; 83.3; 66.6; 50; 33.3; 16.6; 0;];
hex = ['#006687'; '#00b0e9';'#87e2ff';'#feffff';'#ff7575';'#b00000';'#812020'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
map = interp1(vec,raw,linspace(100,0,N),'pchip');

%% Different drug categories
% put all drug categories in a structure
all_categories = struct;

% want to also auto add 025x and 3x to each array 

% categories that are the same regardless of large or small group 
all_categories.rnap = ["RIF","RifT", "OLDRifT","rifapentine","rifampicin","rif6h","rift6h"];
all_categories.dna = ["Lev","MIT","Mox","Olf","OLDMOX",...
    "levofloxacin","moxifloxacin","mitomycin",...
    "ofloxacin","mox6h"];


all_categories.control = ["water","Untreated","NaOH","MeOH","EtOH","DMSO"];

all_categories.unknown = ["Unknown2019","Unknown2238","Unknown2239",...
    "Unknown3825","Unknown3285"];

% different categories depending on fine grouping 
if large_group
    all_categories.protein = ["Amk","Cam","Cla","Dox","Gent","Kan","Lin","Strep","Tet","Tig",...
        "kanamycin","amikacin","chloramphenicol","clarithromycin",...
        "doxycycline","gentamicin","streptomycin","tetracycline","tigecycline",...
        "linezolid","lin6h"];

    all_categories.lipid = ["Tri","Nis","Nig","Mon","Gra","CCCP","PZA","Clz","BDQ","OLDBDQ", ...
            "OLDPZA","OLDClz","clofazimine","carbonyl cyanide 3-chlorophenylhydrazone",...
            "gramicidin","monensin","nigericin","nisin","triclosan","bedaquiline",...
            "pyrazinamide","bdq6h","pza6h","clz6h","Ver","Thi","verapamil","thioridazine" ];

    % for simplicity, create a blank efflux array so that efflux still
    % exists and won't get upset later -- outdated 
    all_categories.efflux = [""];
    
    all_categories.cell_wall = ["Amp","Ctax","Cyc","Del","EMB","ETA",...
        "IMI","INH","Mer","Pre","Van","THL","Cer",...
        "OLDPre","OLDEMB","OLDINH",...
        "ampicillin","cefotaxime","meropenem","isoniazid","ethambutol",...
        "ethionamide","imipenem",...
        "vancomycin","cycloserine","delamanid","pretomanid",...
        "cerulenin","emb6h","inh6h","pre6h"];
    %others = ["PZA", "Ver", "Thi", "Clz","BDQ"];
else 
    % finer settings
    all_categories.peptidoglycan = ["Mer","Amp","Ctax","IMI","Van","Cyc",...
        "meropenem","ampicillin","cefotaxime","imipenem",...
        "vancomycin","cycloserine"];
    all_categories.mycolic_acid = ["INH","ETA","Del","EMB","Pre","THL","Cer",...
        "isoniazid","ethionamide","ethambutol","delamanid","pretomanid",...
        "cerulenin"];
    all_categories.s50_subunit = ["Cla","Cam","Lin"];
    all_categories.s30_subunit = ["Gent","Kan","Amk","Strep","Tet","Dox","Tig"];
    all_categories.atp_synthesis = ["BDQ","CCCP","Clz","Gra","Mon","Nig","Nis","PZA"];
    %efflux = ["Ver","Thi"];
end 

%% loop through and add 025x of each

all_fields = fieldnames(all_categories)';

for pop = all_fields
    current_field = pop{1};
    
    all_drugs_in_category = all_categories.(current_field);
    
    % if not empty 
    if ~(all_drugs_in_category == "")
        for zip = all_drugs_in_category
            % add a 3x version
            x3_ver = strcat(zip,"_3x");
            x025_ver = strcat(zip,"_025x");
            H_ver = strcat(zip," H");
            L_ver = strcat(zip," L");
            all_categories.(current_field) = [all_categories.(current_field) x3_ver x025_ver H_ver L_ver];
        end 
    end 
    
    all_drugs_in_category = all_categories.(current_field);
        
    % same but all uppercase
    if ~(all_drugs_in_category == "")
        for zip = all_drugs_in_category
            % add a 3x version
            upper_ver = upper(zip);
            all_categories.(current_field) = [all_categories.(current_field) upper_ver];
        end 
    end 
    
end 


%%  create additional color struct
% color structure for plotting knn
color_struct = struct;

color_struct.protein = [0 200 0]/255;
color_struct.lipid = [2 0 200]/255;
color_struct.dna = [220 0 0]/255;
color_struct.control = [0 0 0]/255;
color_struct.rnap = [240 240 0]/255;
color_struct.cell_wall = [255 51 255]/255;
color_struct.other = [0 153 153]/255;
color_struct.unknown = [0 153 153]/255;

% smaller settings colors
color_struct.peptidoglycan = [127 26 127]/255;
color_struct.mycolic_acid = [255 51 255]/255;
color_struct.s50_subunit = [0 200 0]/255;
color_struct.s30_subunit = [0 110 0]/255;
color_struct.atp_synthesis = [0 255 255]/255;
color_struct.efflux = [0 150 255]/255; 