model_HMR=importExcelModel('/home/andrei/Downloads/Human-GEM.xlsx',false);

printModelStats(model_HMR,true,true);

% base model looks good


% integrate a patient’s *expression data* ([[RNA seq]]) to in the above generic model

% to run ftINIT, we need:
% 1. to load prepData (for compuation rasons) 
%  for this we need the taskStruct file with the metabolic tasks





% making prepData object
    % making metaboic tasks object

    % I took generic essential human metabolic tasks from RAVEN

taskStruct = parseTaskList('/home/andrei/Downloads/RAVEN_human_metabolicTasks_Essential.txt');


origRefModel = model_HMR;


prepData = prepINITModel(origRefModel, taskStruct);    % step 3 takes ~10 min, step 4 takes ~1h



% making transcrData object


% Load the dataset from a .tsv file
filename = '/home/andrei/Downloads/dataset_RNAseq_Patient1_LasrssonEtAl/dataset_RNAseq_Patient1_LasrssonEtAl.tsv'; % Replace with your actual filename
data = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');

% Extract the gene names for the "genes" component
genes = data.gene_name;

% Create the "tissues" and "celltypes" components
% Set all tissues to "brain" and all cell types to "neuron"
tissues = repmat({'brain'}, height(data), 1); % 1D cell array with "brain"
celltypes = repmat({'neuron'}, height(data), 1); % 1D cell array with "neuron"

% Extract the "levels" component
levels = data.unstranded;

% Create the object with the desired components
geneExpression = struct();
geneExpression.genes = genes;        % 1D cell array of gene names
geneExpression.tissues = tissues;    % 1D cell array of "brain"
geneExpression.celltypes = celltypes; % 1D cell array of "neuron"
geneExpression.levels = levels;      % 1D array of expression levels

% Save the object to a .mat file
% save('geneExpression.mat', 'geneExpression');

%disp('The MATLAB object has been created and saved as "geneExpression.mat".');

disp(geneExpression.genes(69));

transcrData = geneExpression





% final piece: building model

[model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, ...
               fullMipRes] = ...
               ftINIT(prepData, transcrData.tissue, transcrData.celltype, transcrData)



