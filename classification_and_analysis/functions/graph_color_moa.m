function  graph_color_moa(fig)
%less loopy than before, this bad boy ensures everything is exactly the
% same every time 
% jk ended up being p loopy 
% color_maker and graph_helper2 needed to run this

    %now to the scatter plot colors
    ca = get(fig,'CurrentAxes');
    
    %make so that it works with scatter and gscatter
    graph_class = class(get(ca,'Children'));
  
    scat = get(ca,'Children'); %array of all scatters in graph
    
    %if a gscatter object, need to change marker styles
    if strcmp(graph_class,'matlab.graphics.chart.primitive.Line') 
        for row = 1:size(scat,1)
            set(scat(row,1),'Marker','o');
            set(scat(row,1),'MarkerSize',6);
            ecolor = get(scat(row,1),'MarkerEdgeColor');
            set(scat(row,1),'MarkerFaceColor',ecolor);
            set(scat(row,1),'MarkerEdgeColor','none');
        end
    end
    
    %give a grey background color to graph
    set(ca,'Color',[0.8, 0.8, 0.8]);
    
    % array of marker options to rotate through
    markers = ['o'];
    
    
    % sort drug by MoA
    RNAP = ["RifT","RIF","RNAP"];
    lipid = ["PZA","CCCP","Gra","Mon","Nig","Nis","Sulf","BDQ","Clz","Cerold"];
    lipid_2 = ["Tri","Ver","Thi"];
    protein = ["Kan","Amk","Cam","AZT","Cla","Dox","Gent","Strep","Tet","Tig","Lin","protein"];
    dna = ["Dau","Nal","Lev","Mox","MIT","Olf","dna","MOX","Olfold","Nit"];
    cell_wall = ["Mer","Amp","Amox","Ctax","Clex","INH","ETA","IMI","Pip","Cfox", ...
        "Mec","Oxa","PenG","Van","Pre","Cyc","A22","Carb","Del","EMB","cell_wall","Cer","THL"];
    controls = ["MeOH","EtOH","NaOH","water","DMSO","Untreated","control"];
    unknowns = ["Unknown2019", "Unknown2239", "Unknown2911","Unknown3285","Unknown4050"];
    extras = ["INH_control","INH_control_3x"];

    stonybrooks = ["high_692","high_701","high_702","low_691","low_692","low_701","low_702"];
    
    %Set up arrays of similar colors for different MoA
    RNAPcol = color_maker(RNAP, [255 255 0], [0 0 50]);
    lipidcol = color_maker(lipid,[0 83 166],[8 13 35]);
    lipid_2col = color_maker(lipid_2,[0, 211, 211],[100 10 5]);
    proteincol = color_maker(protein,[55 92 0],[10 20 4]);
    dnacol = color_maker(dna,[255 234 196],[0,-38,-62]);
    cell_wallcol = color_maker(cell_wall,[52 0 109],[9 8 14]);
    controlscol = color_maker(controls,[33 33 33],[25 25 25]);
    unknownscol = color_maker(unknowns,[150 20 170],[20 40 -20]);
    extrascol = color_maker(extras,[255 255 255],[0 -15 -50]);
    stonybrookscol = color_maker(stonybrooks,[240 120 21], [-28 20 32]);

    %failsafe if a drug shows up that is not listed
    nogroup = 0;
    ct = 0;
    
    %going through each item in the scatter/line array
    for row = 1:size(scat,1)
        name = get(scat(row,1),'DisplayName');

        
        % want to see if it is dose response data--if so, chop off ending
        % so that it can sort into an exisitng category
        find_underscore_w_dose = '([_][0-9])\w+';
        % that regexp will look to see if there is an underscore followed
        % by a number, indicating a dose         
        % if this is not empty, then it is dose data
        is_dose = ~isempty(regexp(name,find_underscore_w_dose));
        
        if is_dose
            underscore_loc = regexp(name,find_underscore_w_dose);
            % find what the dose is for coloring purposes
            dose = name(underscore_loc+1:end);
            % want to extract before that underscore and continue with that
            name = name(1:underscore_loc-1);
            % we also want to add a border or something to differentiate
            % plotting 3x and 025x, etc. 
            if strcmp(dose,'025x')
                set(scat(row,1),'MarkerEdgeColor',[.12, .12, .12])
            elseif strcmp(dose,'05x')
                set(scat(row,1),'MarkerEdgeColor',[.5, .5, .5])
            elseif strcmp(dose,'075x')
                set(scat(row,1),'MarkerEdgeColor',[1 .75, .6])
            elseif strcmp(dose,'1x')
                set(scat(row,1),'MarkerEdgeColor',[1, 1, 1])
            end
            % no border for 3x
            
        end
        
        % ... if they belong to a group, use the correct color array to
        % change the color/marker
        set(scat(row,1),'Marker','o')
        if ismember(name,RNAP) 
            set(scat(row,1),'MarkerFaceColor',[240 240 0]/255)
            
        elseif ismember(name, lipid)           
            set(scat(row,1),'MarkerFaceColor',[2 0 200]/255)
            
        elseif ismember(name, lipid_2)
            set(scat(row,1),'MarkerFaceColor',[2 0 200]/255)
            
        elseif ismember(name, protein)
            set(scat(row,1),'MarkerFaceColor',[0 200 0]/255)
            
        elseif ismember(name, dna)
            set(scat(row,1),'MarkerFaceColor',[220 0 0]/255)
            
        elseif ismember(name, cell_wall)
            set(scat(row,1),'MarkerFaceColor',[255 51 255]/255)
            
        elseif ismember(name,stonybrooks)
            set(scat(row,1),'MarkerFaceColor',[240 240 0]/255)
            
        elseif ismember(name, controls)
            set(scat(row,1),'MarkerFaceColor',[.09, .09, .09])
            
        elseif ismember(name,unknowns)
            set(scat(row,1),'MarkerFaceColor',[240 240 0]/255)
            
        elseif ismember(name,extras)
            ind = find(strcmpi(extras,name));
            graph_helper(scat(row,1),extrascol,ind,markers)
            set(scat(row,1),'MarkerEdgeColor',[.09, .09, .09])
            
        else 
             set(scat(row,1),'MarkerFaceColor',[1, 1, 1])
             nogroup = nogroup+1;
             disp(name)
             set(scat(row,1),'Marker',markers(nogroup));
        end
        
  
    end
    
   
end



function graph_helper(scat, rgb, ind,markers)
% helper function for graph_color_s
% inputs: scat: a line or scatter object
% rgb: array of rgb values from 0-255
% ind: the index of where we currently are in the rgb array
% rgb(ind,1) corresponds to the r value at this index
% markers: array of markers to swap through

    % obtain correct rgb values
    r = rgb(ind,1);
    g = rgb(ind,2);
    b = rgb(ind,3);
 
    % will cause to loop through markers array when it gets too big
    while ind > size(markers,2)
        ind = ind - size(markers,2);
    end
    
    
    % switch from 0 - 255 scale to 0 - 1 scale and set value
    set(scat,'MarkerFaceColor',[r/255, g/255, b/255])
    
    % p is smaller than the rest by deault, need to make it bigger
    % size works slightly differently for scatter vs line objects
    % need to check the class before changing the size
    if strcmp(markers(ind),'p')
        classy = class(scat);
        if strcmp(classy,'matlab.graphics.chart.primitive.Scatter')
            set(scat,'SizeData',60)
        else
            set(scat,'MarkerSize',10);
        end
    end
    
    %looping through markers array
    set(scat,'Marker',markers(ind));
    
end


function color_ary = color_maker(ary, rgb, rgbadd)
% helper function for graph_color_s
    a = [];
    r = rgb(1);
    g = rgb(2);
    b = rgb(3);
    
    radd = rgbadd(1);
    gadd = rgbadd(2);
    badd = rgbadd(3);
    
    for i = 1:size(ary,2)

        r = r + radd;
        if r > 255 
            r = 255;
        elseif r < 0
            r = 0;
        end
        
        g = g + gadd;
        if g > 255 
            g = 255;
        elseif g < 0
            g = 0;
        end
        
        b = b + badd;
        if b > 255 
            b = 255;
        elseif b < 0
            b = 0;
        end
        a = [a; [r g b]];
        
    end
    
    color_ary = a;
end


