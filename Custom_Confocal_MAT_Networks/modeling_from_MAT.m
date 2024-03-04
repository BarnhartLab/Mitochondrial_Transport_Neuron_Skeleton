%add network tools to path
addpath('/Users/jordankalai/Documents/Barnhart_Lab_Code/Architecture_analysis/networktools-master')
addpath('/Users/jordankalai/Documents/Barnhart_Lab_Code/Architecture_analysis/dendriticTrees-master')


%% directory name containing multiple .mat files
% each .mat file needs to contain the following saved variables:
% NT = a network object, already arranged to be a directed tree, rootnode set
% umperpx = conversion factor, micrometers per pixel
dirname = '/Users/jordankalai/Documents/Barnhart_Lab/MCFO/HS/cleanedNetWorkspaces/';
%dirname = '/Users/jordankalai/Documents/Barnhart_Lab_Code/Architecture_analysis/architecture_pythonCodes/pythonDemoScripts/asym_sim_skel/';
%% load in all .mat files within the directory
% and save data
files = dir([dirname '*.mat'])
filenames = {files.name};

clear allNetworks
for fc = 1:length(filenames)
        
    fname = [dirname filenames{fc}];
    disp(fname)
    load(fname);
    
    allNetworks(fc) = NT;    
    allumperpx(fc) = umperpx;
end

%% variables needed for exporting
cell_type = []; 
number = [];
alpha_fin = []; %try 2, 3/2, and 3
beta_fin = [];
splitting_rule = [];
sister_subtree_asym = [];
distal_enrichment = [];
%% loop to get all measurements for each cell, starting with converting to net file, then running model predictions for each splitting pattern
for nc = 1:length(allNetworks)    
    NT = allNetworks(nc);

    trunkedge = NT.nodeedges(NT.rootnode);
    %% define which are the distal edges
    % distance of each edge from root node
    [edgedist] = getEdgeDist(NT,trunkedge);
    edgedist = edgedist*allumperpx(nc);
    
    % distal edges are defined as > this prefactor times max distance to root
    distalprefactor= 0.75;

    % distal is > 75% of max distance   
    distcutoff(1) = max(edgedist)*distalprefactor;
    distcutoff(2) = inf;
    distind = find(edgedist>distcutoff(1) & edgedist<distcutoff(2));  
    %% run model prediction for each splitting rule for different alphas and for different betas
    % this will be a triple nested loop, first layer is splitting rule,
    % second layer is alpha
    %third layer is beta
    % how to split radii. Possible values:
    % LV = split so that subtree L ~ V; does *not* use rm
    % LD = split so that branch radius ~ subtree length / depth
    % L = split so that branch radius ~ subtree length
    % equal = split equally
    splittype = {'LV', 'LD', 'equal', 'L'};
    alpha = [2.2, 2, 3/2, 3]; % alpha value relating parent and daughter width: r0^a = r1^a + r2^a
    rm = 0; % `minimal radius' (from Liao paper) such that r0^a + rm^a = r1^a + r2^a
    rtrunk = 3; % trunk radius in um
    %don't know what to set image radius too
    beta = [0, 0.8, 1.5, 2, 3]; % beta is the exponent relating stopping rate and trunk width: ks ~ 1/r^beta
    for s = 1: length(splittype)
        for a = 1: length(alpha)
            for b = 1: length(beta)
                number = [number nc];
                cell_type = [cell_type "HS"]; %can change based on cell type
                alpha_fin = [alpha_fin alpha(a)]; %try 2, 3/2, and 3
                beta_fin = [beta_fin beta(b)];
                split_rule_add = string(splittype{s});
                splitting_rule = [splitting_rule split_rule_add];
                %set empty stL, stV, and stD parameters
                stL = []; stV = []; stD = [];
                % define radius of trunk
                rtrunk = NT.edgewidth{trunkedge}(1,1); %set as the initial radius in the trunk
                % calculate some basic statistics
                [stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,alpha(a),'LV');        
                switch(splittype{s}) %calculate volume and radii based on each case in loop
                    case('LV')
                        radii = setRadiiFromSubtreeInfo(NT,trunkedge,alpha(a),rtrunk,muvals);   
                        stV = stEta.*radii.^2;
                    case('LD')
                        [radii,stV] = setRadiiWithRm(NT,trunkedge,alpha(a),rm,rtrunk,stL./stD);
                    case('equal')       
                        [radii,stV] = setRadiiWithRm(NT,trunkedge,alpha(a),rm,rtrunk,ones(size(stL)));
                    case('L')        
                        [radii,stV] = setRadiiWithRm(NT,trunkedge,alpha(a),rm,rtrunk,stL);
                end 
                %% for the network with measured edge widths, calculate linear mitochondria concentrations
                %Ask Erin what she wants for these values
                % ratio of stopped to walking rates
                kskw = 1000;
                % density of mitos in trunk (arbitrary units)
                rhoWtrunk = 1;
                
                trunkedge = NT.nodeedges(NT.rootnode,1);
                % rhoSvals = stationary density on each edge
                % rhoWvals = mobile density on each edge
                % stM = average mito density in the subtree below (and including) each edge
                [rhoSvals,rhoWvals,stM] = setMitoConcFromRadii(NT,trunkedge,beta(b),kskw,rhoWtrunk,radii);
                rhovals = rhoSvals + rhoWvals;
                
                %% Calculate asymmetry in mito concentrations below each junction
                % calculate asymmetry at each junction
                % indexed by edge leading to that junction
                [asymmetry,nontermjunc] = getAsymmetry(NT,stM, stV);
                
                % get the average asymmetry
                avgasym = sqrt(sum(asymmetry(nontermjunc).^2)/length(nontermjunc));
                sister_subtree_asym = [sister_subtree_asym avgasym];

                %% Calculate distal enrichment of mitochondrial density
                Mdist = sum(rhovals(distind).*NT.edgelens(distind));
                Vdist = sum(NT.edgelens(distind).*radii(distind).^2);
                distenrich= (Mdist/Vdist) / (rhovals(trunkedge)/radii(trunkedge).^2);
                distal_enrichment = [distal_enrichment distenrich];
                [alpha(a) beta(b) splittype{s}]
            end
        end
    end

    
end
%% after looping through all .swcs then add calculations to a singular array as columns to then export as csv
length(cell_type) %right value
length(number) %forgot to add numbers
length(alpha_fin)
length(beta_fin)
length(splitting_rule)
length(sister_subtree_asym)
output_table = table(cell_type(:), number(:), alpha_fin(:), beta_fin(:), splitting_rule(:), sister_subtree_asym(:), distal_enrichment(:));

%set column names
output_table.Properties.VariableNames = ["cell", "n", "alpha", "beta", "splitting_rule", "sister_subtree_asym", "mito_distal_enrichment"];
writetable(output_table,'/Users/jordankalai/Documents/Barnhart_Lab/Paper Figures/simulated_model_predictions.csv','Delimiter',',')