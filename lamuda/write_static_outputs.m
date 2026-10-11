function write_static_outputs(resultsTable,cfg)
%WRITE_STATIC_OUTPUTS Write the numerical benchmark table.
writetable(resultsTable,fullfile(cfg.resultsFolder,'static_benchmark_results.csv'));
end

