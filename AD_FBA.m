% tutorial1
%   This is a short introduction that shows how to load a genome-scale
%   metabolic model (GEM), set reaction constraints, objective function and
%   perform an optimization through flux balance analysis (FBA). The
%   resulting fluxes are visualized and exported to a PDF file.
%   For a more detailed description of the individual functions, see
%   [raven_directory]/doc/index.html.
%   A GEM for the filamentous fungus Penicillium chrysogenum is used in
%   this tutorial. The model can be found in a Microsoft Excel file under
%   the name "iAL1006 v1.00.xlsx" and in SBML file "iAL1006 v1.00.xml".
%   See Tutorial 1 in "RAVEN tutorials.docx" for more details.

%Import the model from Excel. This function performs a number of checks
%regarding the model structure (such as for incorrectly written equations
%or illegal characters). In this structure there is only one warning; that
%the formula for the metabolite LPE could not be parsed. The "false" flag
%imports a model with exchange reactions in their "closed" form. This makes
%the model unsuited for modelling, but it is useful for some quality
%control steps.

model_full_O2_intake=importModel('/home/andrei/Downloads/12918_2012_1050_MOESM2_ESM/hippocampus.xml',false);
model_half_O2_intake=importModel('/home/andrei/Downloads/12918_2012_1050_MOESM2_ESM/hippocampus.xml',false);

%The Excel interface is supposed to work in all the syDownloadsstems (Windows, Unix,
%macOS), but upon any problems, the model can be imported from SBML format
%instead. However, in such case the user would not be able to run Tutorials
%2-4, which involve the editing of Excel files. Run the command below
%(remove "%" sign) instead, if having such problem:
%model=importModel('iAL1006 v1.00.xml',false);

%The following function prints some properties of the model. The two "true"
%flags say that it should also list potential problems such as dead-end
%reactions or unconnected metabolites.
% printModelStats(model,true,true);

%As can be seen the model contains 1129 reactions, 920 metabolites and
%1045 genes

%Most modelling approaches using GEMs are based on the mass balancing
%around the internal metabolites in the system. However, in order for the
%system to uptake or excrete metabolites, some metabolites have been
%defined as "unconstrained". In order to simulate something, those
%metabolites have to be removed from the model. The function simplifyModel
%is a general-purpose function for making models smaller. This includes the
%options such as grouping linear reactions and deleting reactions which
%cannot carry flux. Here it is chosen to delete the exchange metabolites,
%all reactions that are constrained to zero (mainly uptake of non-standard
%carbon sources), and all reactions that cannot carry flux (mainly
%reactions that were dependent on any of those non-standard carbons
%sources).

% model=simplifyModel(model,true,false,true,true);

% for hippocampus model, no element is removed 

% explore Hippocampus model properties

%%% Getting reactions and metabolite lists

% display(model.compNames);
% display(model_full_O2_intake.metNames);

% all_mets_hippocampus = model_full_O2_intake.metNames
% save('all_mets_hippocampus.mat', 'all_mets_hippocampus');


% all_reactions_hippocampus = model_full_O2_intake.rxns


% save('all_reactions_hippocampus.mat', 'all_reactions_hippocampus');

% display(model_full_O2_intake.comps);
% display(model_full_O2_intake.compNames);
% display(model.mets)

% mets_extracellular = model.mets(endsWith(model.mets, '_e'));

% display(mets_extracellular);

% find reactions containing >=1 from the above list

% display(model.rxns(1).metabolite);

%As a first way of validating the model, calculate the theoretical yield of
%carbon dioxide from glucose. The supplied model already allows for uptake
%of phosphate, sulfate, NH3, O2 and the co-factor precursors thiamin and
%pimelate. The setParam function is used for setting constraints,
%reversibility and objective function coefficients.
%Set the uptake of glucose to be no more than 1 mmol/gDW/h and no uptake of
%ethanol.
%model=setParam(model,'ub',{'glcIN' 'etohIN'},[1 0]);


model_full_O2_intake=setParam(model_full_O2_intake,'eq',{'O2t'},[1]);

model_half_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.5]);

model_0125_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.125]);
model_025_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.25]);
model_0375_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.375]);

model_0625_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.625]);
model_075_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.75]);
model_0875_O2_intake=setParam(model_half_O2_intake,'eq',{'O2t'},[0.875]);

list_models = {model_0125_O2_intake,model_025_O2_intake, model_0375_O2_intake, model_half_O2_intake, model_0625_O2_intake,model_075_O2_intake, model_0875_O2_intake, model_full_O2_intake};

% list_metabolites_of_interest = {"atp_m", "glu_DASH_DASH_DASH_L_e", "aacoa_c" }

list_indeces_reactions_of_interest = {209, 961, 34, 153, 189, 193, 541, 542, 543, 544, 545, 546, 790, 791, 792, 835, 1036, 114, 115};

O2_values = {0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1};

n = length(list_indeces_reactions_of_interest); % Number of rows
m = length(O2_values); % Number of columns


% display(list_models(1));
% aux_model = setParam(model_half_O2_intake,'eq',{'O2t'},[0.875]);

% sol_aux=solveLP(aux_model);

% aux_O2value = O2_values{1};

% aux_model = setParam(model_full_O2_intake,'eq',{'O2t'},aux_O2value);
% sol_aux=solveLP(aux_model);
% aux_reaction_idx = list_indeces_reactions_of_interest(1);
% flux = sol_aux.x(aux_reaction_idx);
% display(flux)


matrix_fluxes_reactionsByModels = zeros(n, m);

% Nested for loop to iterate over each cell
for i = 1:n  % Loop over rows
    for j = 1:m  % Loop over columns
        % Perform some operation (e.g., sum of indices)        

        aux_O2value = O2_values{j};
        aux_model = setParam(model_full_O2_intake,'eq',{'O2t'},aux_O2value);
        sol_aux=solveLP(aux_model);
        aux_reaction_idx = list_indeces_reactions_of_interest{i};
        flux = sol_aux.x(aux_reaction_idx);
        
        % Assign the value to the corresponding cell
        matrix_fluxes_reactionsByModels(i, j) = flux; 
    end
end



save('matrix_fluxes_reactionsByModels.mat', 'matrix_fluxes_reactionsByModels');


%Set the objective for the simulation to maximize ATP synthase in cytosol
model_full_O2_intake=setParam(model_full_O2_intake,'obj',{'ATPS4m'},1);
model_half_O2_intake=setParam(model_half_O2_intake,'obj',{'ATPS4m'},1);

%The problem can now be solved using linear programming. The solveLP
%function takes a model and solves the linear programming problem defined
%by the constraints and the objective value coefficients.
sol_full=solveLP(model_full_O2_intake);
sol_half=solveLP(model_half_O2_intake);

% LEFT HERE: figure out 
display(sol_full.x(4))

display(length(sol_full.x))

display(length(model_0125_O2_intake.rxnNames))



%The results show that the growth rate is 0.084/h and that the system now
%also requires sulfate, phosphate, NH3, thiamin and pimelate. Compare this
%to the growth on ethanol instead of glucose. Use three times the molar
%flux of ethanol since it contains 2 carbons rather than 6. 'eq' means
%'equal to' and sets the lower and upper bound to the same value.

% modelETH=setParam(model,'eq',{'glcIN' 'etohIN'},[0 3]);
% solETH=solveLP(modelETH,1);
% printFluxes(modelETH, solETH.x, true, 10^-7);

%Investigate how metabolism changes between the two carbon sources.
%followChanged takes two flux distributions and lets the user select which
%reactions should be printed. Here the reactions are shown that differ with
%more than 50%, has a flux higher than 0.5 mmol/gDW/h and an absolute
%difference higher than 0.5 mmol/gDW/h.
followChanged(model_half_O2_intake,sol_full.x,sol_half.x, 50, 0.5, 0.5);

%There are 65 such reactions. By studying them one can start to get an idea
%about where the key changes occur. Visualization can help a lot in this
%regard. One can investigate how ATP metabolism changes by running the
%following command:
followChanged(model_half_O2_intake,sol_full.x,sol_half.x, 30, 0.4, 0.4,{'ATP'});

%See that on glucose ATP is generated in glycolysis but on ethanol it seems
%to have to do with acetate and so on. This allows the user to look further
%and further down until one understands the underlying flux redistributions
%that give rise to different phenotypes.

%The fluxes can also be visualized on a metabolic map. Green corresponds to
%reactions which are more used for growth on glucose, and red are reactions
%which are more used for growth on ethanol. Open the "GLCvsETH.pdf" file to
%be able to zoom in on individual reactions.
load 'pcPathway.mat' pathway;

% TODO: if no colors displayed, play with cutoffs (10^)
drawMap('Full vs Half Oxygen intake',pathway,model_full_O2_intake,sol_full.x,sol_half.x,model_half_O2_intake,'fullVsHalfO2Intake.pdf',10^-20);
