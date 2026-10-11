function export_case_data(data,cfg)
%EXPORT_CASE_DATA Export original/converted bus and branch data to CSV.

% Create one audit folder for the current converted system.
caseFolder = fullfile(cfg.caseDataFolder,data.key);
if exist(caseFolder,'dir') ~= 7
    mkdir(caseFolder);
end
% Export the converted bus and branch tables for independent inspection.
writetable(data.busTable,fullfile(caseFolder,'bus_data.csv'));
writetable(data.branchTable,fullfile(caseFolder,'branch_data.csv'));

% Write a text summary of the conversion assumptions and selected buses.
fid = fopen(fullfile(caseFolder,'conversion_notes.txt'),'w');
if fid < 0
    error('Could not create conversion notes for %s.',data.displayName);
end
% Always close the file, including when a later fprintf operation fails.
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'System: %s\n',data.displayName);
fprintf(fid,'Official case function: %s\n',data.caseFunction);
fprintf(fid,'Official source: %s\n',data.sourceURL);
fprintf(fid,'baseMVA: %.12g\n',data.baseMVA);
fprintf(fid,'AC power-flow source: %s\n',data.powerFlowMessage);
fprintf(fid,'VSG buses (internal): %s\n',mat2str(data.vsgBuses.'));
fprintf(fid,'CPL buses (internal): %s\n',mat2str(data.cplBuses.'));
fprintf(fid,'Zero-injection algebraic buses: %s\n', ...
    mat2str(data.zeroInjectionBuses.'));
fprintf(fid,'Loss surplus removed from positive injections [p.u.]: %.12g\n', ...
    data.acLossSurplusPU);
fprintf(fid,'Positive-injection rescaling: %.12g\n',data.sourceRescaling);
fprintf(fid,'Balanced injection residual: %.3e\n',sum(data.varpiNominal));
fprintf(fid,'Conversion: fixed voltage, R=0, Bsh=0, phase shifts ignored.\n');
fprintf(fid,'Line coefficient: a_e=V_i V_j/(x_e |tap_e|).\n');
end
