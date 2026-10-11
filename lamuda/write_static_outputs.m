function write_static_outputs(resultsTable,cfg)
%WRITE_STATIC_OUTPUTS Write diagnostics and the manuscript loading-factor table.
writetable(resultsTable,fullfile(cfg.resultsFolder,'static_benchmark_results.csv'));

% The paper compares the four loading limits directly.  Keep this table
% separate from the full CSV so it can be copied into the manuscript.
tablePath = fullfile(cfg.resultsFolder,'critical_loading_factors.tex');
fid = fopen(tablePath,'w');
assert(fid >= 0,'Cannot create the LaTeX table.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
rowEnd = [char(92),char(92),newline];

fprintf(fid,'\\begin{table}[!t]\n');
fprintf(fid,'\\caption{Critical loading factors of the existing and proposed tests.}\n');
fprintf(fid,'\\label{tab:critical_loading_factors}\n');
fprintf(fid,'\\centering\\footnotesize\n');
fprintf(fid,'\\setlength{\\tabcolsep}{2pt}\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.20}\n');
fprintf(fid,'\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}lcccc@{}}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'& \\multicolumn{2}{c}{2-norm tests} & \\multicolumn{2}{c}{$\\infty$-norm tests} %s',rowEnd);
fprintf(fid,'\\cmidrule(lr){2-3}\n');
fprintf(fid,'\\cmidrule(lr){4-5}\n');
fprintf(fid,'System & $\\lambda_1$ & $\\lambda_3$ & $\\lambda_2$ & $\\lambda_4$ %s',rowEnd);
fprintf(fid,'\\midrule\n');
for k = 1:height(resultsTable)
    name = tex_escape(resultsTable.System(k));
    fprintf(fid,'%s & %.4f & %.4f & %.4f & %.4f %s',name, ...
        resultsTable.lambda_1(k),resultsTable.lambda_3(k), ...
        resultsTable.lambda_2(k),resultsTable.lambda_4(k),rowEnd);
end
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular*}\n');
fprintf(fid,'\\end{table}\n');
clear cleanup;

% Preserve the old filename as a compatibility copy for scripts that used
% the earlier output name, but keep its contents consistent with the final
% four-limit comparison.
copyfile(tablePath,fullfile(cfg.resultsFolder,'critical_ratio_comparison.tex'),'f');
end

function value = tex_escape(value)
value = char(string(value));
value = strrep(value,'_','\\_');
value = strrep(value,'&','\\&');
end

