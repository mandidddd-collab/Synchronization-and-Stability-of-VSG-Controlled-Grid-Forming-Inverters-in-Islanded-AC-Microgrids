function selected = prepare(cfg)
%PREPARE Use only the bundled MATPOWER case files.
% Create all output folders before any calculation writes a file.
folders = {cfg.resultsFolder,cfg.figureFolder,cfg.caseDataFolder};
for k = 1:numel(folders)
    if ~isfolder(folders{k}), mkdir(folders{k}); end
end
% Read the bundled catalog and select the cases requested in cfg.caseKeys.
catalog = normbench.case_catalog();
selected = repmat(catalog(1),0,1);
for k = 1:numel(cfg.caseKeys)
    index = find(strcmp({catalog.key},cfg.caseKeys{k}),1);
    assert(~isempty(index),'Unknown case: %s',cfg.caseKeys{k});
    spec = catalog(index);
    assert(isfile(fullfile(cfg.dataFolder,[spec.functionName,'.m'])), ...
        'Missing bundled case file: %s',spec.functionName);
    selected(end+1,1) = spec; %#ok<AGROW>
end
% Put the bundled data folder first on the MATLAB search path.
addpath(cfg.dataFolder,'-begin');
% Check that the local case file, not another installation, is selected.
for k = 1:numel(selected)
    assert(strcmpi(which(selected(k).functionName), ...
        fullfile(cfg.dataFolder,[selected(k).functionName,'.m'])), ...
        'A different %s shadows the bundled data.',selected(k).functionName);
end
end
