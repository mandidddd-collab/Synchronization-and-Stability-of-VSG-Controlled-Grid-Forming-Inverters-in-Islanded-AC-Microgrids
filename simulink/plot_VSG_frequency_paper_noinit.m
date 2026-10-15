% Plot the VSG30_38_Frequency simulation results.
% Run Simulink first, then run plot_VSG_frequency_paper_noinit.

%% 1. Plot settings.
cfg.DataName = 'f_vsg30_38_Hz';
cfg.YLim = [49.90 50.40];          % Frequency limits in Hz; alternative: [49.6 51.0].
cfg.YTickStep = 0.10;              % Frequency tick spacing in Hz; alternative: 0.20.
cfg.XLim = [];                    % [] uses the full recorded time; alternatively, set [0 18.1].
cfg.XTickStep = [];               % [] uses automatic ticks; use 3 for a three-second spacing.
cfg.FigureSizeCM = [16.8 8.6];     % Figure width and height, in cm.
cfg.FontSize = 12;
cfg.LineWidth = 0.75;
cfg.Channels = 1:9;               % Channel order: Bus 30,31,...,38.
cfg.ShowLegend = false;           % Hide the legend.
cfg.StageStart = [0 0.05 3 6 9 12 15];
cfg.StageRatio = [0.2 0.4 0.6 0.8 1.0 1.2 1.35];
cfg.OutputDir = fullfile(fileparts(mfilename('fullpath')), 'VSG_paper_figures');
cfg.FileStem = 'VSG_frequency_paper';

% Optional: override settings using a userOptions structure.
if exist('userOptions','var') && ~isempty(userOptions)
    assert(isstruct(userOptions) && isscalar(userOptions), ...
        'userOptions必须是标量结构体。');
    keys = fieldnames(userOptions);
    for k = 1:numel(keys)
        assert(isfield(cfg, keys{k}), '未知设置：%s', keys{k});
        cfg.(keys{k}) = userOptions.(keys{k});
    end
end

%% 2. Read the original log directly without helper functions.
% Prefer manually supplied data; otherwise find a workspace log or simulation output.
% Keep automatically detected logs in frequencyLog to avoid stale input on later runs.
if exist('data','var') && ~isempty(data)
    frequencyLog = data;
    dataSource = '工作区 data';
elseif exist(cfg.DataName,'var')
    frequencyLog = eval(cfg.DataName);
    dataSource = cfg.DataName;
else
    frequencyLog = [];
    candidates = {'out','simOut','simout','simulationOutput'};
    for k = 1:numel(candidates)
        name = candidates{k};
        if ~exist(name,'var'), continue; end
        obj = eval(name);
        if isa(obj,'Simulink.SimulationOutput') && any(strcmp(who(obj),cfg.DataName))
            frequencyLog = get(obj,cfg.DataName);
            dataSource = [name '.' cfg.DataName];
            break;
        elseif isstruct(obj) && isscalar(obj) && isfield(obj,cfg.DataName)
            frequencyLog = obj.(cfg.DataName);
            dataSource = [name '.' cfg.DataName];
            break;
        end
    end
    assert(~isempty(frequencyLog), ...
        '未找到%s。请先运行模型，或设置data=你的日志后再运行脚本。',cfg.DataName);
end
if isa(frequencyLog,'Simulink.SimulationOutput')
    assert(any(strcmp(who(frequencyLog),cfg.DataName)), ...
        'SimulationOutput中没有变量%s。',cfg.DataName);
    frequencyLog = get(frequencyLog,cfg.DataName);
end
if isa(frequencyLog,'Simulink.SimulationData.Signal')
    frequencyLog = frequencyLog.Values;
end

% Convert Timeseries or Structure With Time logs to a time-sample-by-channel matrix.
if isa(frequencyLog,'timeseries')
    t = double(frequencyLog.Time(:));
    values = {frequencyLog.Data};
elseif isstruct(frequencyLog) && isscalar(frequencyLog) && ...
        isfield(frequencyLog,'time') && isfield(frequencyLog,'signals')
    t = double(frequencyLog.time(:));
    assert(~isempty(t),'请使用Timeseries或Structure With Time格式保存数据。');
    values = {frequencyLog.signals.values};
else
    error('数据类型%s不支持，请使用timeseries或含time的日志结构体。',class(frequencyLog));
end
cols = cell(1,numel(values));
for k = 1:numel(values)
    samples = squeeze(double(values{k}));
    if size(samples,1)==numel(t)
        cols{k} = reshape(samples,numel(t),[]);
    elseif ismatrix(samples) && size(samples,2)==numel(t)
        cols{k} = samples.';
    else
        error('频率数组维度与时间长度不对应。');
    end
end
f_Hz = cat(2,cols{:});
assert(~isempty(t) && numel(t) > 1, '日志为空或样本不足；请先完成仿真。');
assert(all(isfinite(t)) && all(diff(t) >= 0) && t(end) > t(1), ...
    '时间必须是有效、非递减的秒数，且至少含两个不同时刻。');
assert(isreal(f_Hz), '频率数据应为实数。');
assert(all(cfg.Channels >= 1 & cfg.Channels <= size(f_Hz,2) & ...
    cfg.Channels == floor(cfg.Channels)), '所选通道超出日志范围。');
assert(numel(cfg.YLim) == 2 && diff(cfg.YLim) > 0, 'YLim 必须递增。');
assert(isscalar(cfg.YTickStep) && cfg.YTickStep > 0, 'YTickStep 必须为正数。');
assert(numel(cfg.FigureSizeCM) == 2 && all(cfg.FigureSizeCM > 0), ...
    'FigureSizeCM 必须为两个正数。');

% The logged frequency is already in Hz; do not multiply it by 50 again.
y0 = f_Hz(1:min(100,size(f_Hz,1)), :);
y0 = y0(isfinite(y0));
if ~isempty(y0) && median(abs(y0)) < 5
    error(['输入似乎不是 Hz 频率。请使用 f_vsg30_38_Hz，', ...
        '不要直接传入 w2、w3 等标幺频率。']);
end
if any(~isfinite(f_Hz(:)))
    warning('日志包含 NaN/Inf；不填补或平滑这些点，请检查仿真是否正常结束。');
end
fprintf('数据来源：%s；%d 个样本，%d 个通道；时间 %.6g～%.6g s。\n', ...
    dataSource, numel(t), size(f_Hz,2), t(1), t(end));

%% 3. Use the loading times and ratios specified in this script.
stageStart = cfg.StageStart;
stageRatio = cfg.StageRatio;
scheduleSource = 'StageStart / StageRatio in this script';
assert(numel(stageStart) == numel(stageRatio) && ...
    all(isfinite(stageStart)) && all(isfinite(stageRatio)) && ...
    all(diff(stageStart) >= 0), '加载时刻与倍率不匹配。');
stageStart = stageStart(:).';
stageRatio = stageRatio(:).';
fprintf('阶段标注依据：%s。\n', scheduleSource);
fprintf('请确保日志与该加载时序来自同一次设置；脚本不改动加载时刻。\n');

if isempty(cfg.XLim)
    xLim = [t(1), t(end)];
else
    xLim = cfg.XLim;
    assert(numel(xLim) == 2 && diff(xLim) > 0, 'XLim 必须递增。');
    assert(xLim(1) < t(end) && xLim(2) > t(1), 'XLim 与数据时间不重叠。');
end
xSpan = diff(xLim);
ySpan = diff(cfg.YLim);

%% 4. White background, Times New Roman fonts, thin lines and a light grid.
fig = figure('Color','w', 'Name','VSG frequency - paper style', ...
    'NumberTitle','off', 'Units','centimeters', ...
    'Position',[3 3 cfg.FigureSizeCM], 'Renderer','painters');
ax = axes('Parent',fig, 'Units','normalized', ...
    'Position',[0.12 0.17 0.855 0.79]);
hold(ax,'on');
set(ax, 'Color','w', 'Box','on', 'LineWidth',0.70, ...
    'FontName','Times New Roman', 'FontSize',cfg.FontSize, ...
    'XColor','k', 'YColor','k', 'TickDir','in', 'TickLength',[0.012 0.012], ...
    'XMinorTick','on', 'YMinorTick','on', 'Layer','top', ...
    'XGrid','on', 'YGrid','on', 'GridColor',[0.65 0.65 0.65], ...
    'GridAlpha',0.24, 'GridLineStyle','-');
xlim(ax,xLim);
ylim(ax,cfg.YLim);
if ~isempty(cfg.XTickStep)
    assert(cfg.XTickStep > 0, 'XTickStep 必须为正数。');
    xticks(ax, ceil(xLim(1)/cfg.XTickStep)*cfg.XTickStep : ...
        cfg.XTickStep : xLim(2));
end
yTicks = cfg.YLim(1):cfg.YTickStep:cfg.YLim(2);
yticks(ax,yTicks);
yticklabels(ax,arrayfun(@(v)sprintf('%g',v), yTicks, 'UniformOutput',false));

% Nominal-frequency reference line, excluded from the legend.
plot(ax,xLim,[50 50], ':', 'Color',[0.48 0.48 0.48], ...
    'LineWidth',0.60, 'HandleVisibility','off');

% Draw step markers behind the curves; use red dashes for 1.35lambda4.
for k = 2:numel(stageStart)
    if stageStart(k) > xLim(1) && stageStart(k) < xLim(2)
        color = [0.55 0.55 0.55];
        if abs(stageRatio(k)-1.35) < 1e-8
            color = [0.85 0.10 0.10];
        end
        plot(ax,[stageStart(k) stageStart(k)],cfg.YLim, '--', ...
            'Color',color, 'LineWidth',0.65, 'HandleVisibility','off');
    end
end

% Nine distinct colors; retain every curve even when synchronization is lost.
colors = [0.000 0.447 0.741;
          0.850 0.325 0.098;
          0.780 0.580 0.120;
          0.494 0.184 0.556;
          0.466 0.674 0.188;
          0.301 0.650 0.833;
          0.635 0.078 0.184;
          0.000 0.550 0.500;
          0.750 0.300 0.650];
h = gobjects(numel(cfg.Channels),1);
for k = 1:numel(cfg.Channels)
    c = cfg.Channels(k);
    colorIndex = mod(c-1,size(colors,1))+1;
    h(k) = plot(ax,t,f_Hz(:,c), 'Color',colors(colorIndex,:), ...
        'LineWidth',cfg.LineWidth, 'DisplayName',sprintf('Bus %d',29+c));
end
xlabel(ax,'Time (s)', 'FontName','Times New Roman', ...
    'FontSize',cfg.FontSize+1, 'FontWeight','bold');
ylabel(ax,'Frequency (Hz)', 'FontName','Times New Roman', ...
    'FontSize',cfg.FontSize+1, 'FontWeight','bold');

%% 5. Label the six load increases; describe the initial 0.2lambda4 state in the paper.
% Retain all data and step times; omit only the initial-stage text label.
stageEnd = [stageStart(2:end), max(t(end),xLim(2))];
for k = 2:numel(stageStart)
    left = max(stageStart(k),xLim(1));
    right = min(stageEnd(k),xLim(2));
    if right <= left
        continue;
    end
    xText = (left+right)/2;
    yText = cfg.YLim(2)-0.025*ySpan;
    if abs(stageRatio(k)-1) < 1e-8
        label = '\lambda=\lambda_4';
    else
        label = sprintf('\\lambda=%g\\lambda_4',stageRatio(k));
    end
    textColor = [0 0 0];
    if abs(stageRatio(k)-1.35) < 1e-8
        textColor = [0.85 0.10 0.10];
    end
    text(ax,xText,yText,label, 'Interpreter','tex', ...
        'FontName','Times New Roman', 'FontSize',cfg.FontSize, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'Color',textColor, 'Clipping','on');
end
if cfg.ShowLegend
    legend(ax,h,'Location','southoutside','NumColumns',3, ...
        'Box','off','FontName','Times New Roman','FontSize',cfg.FontSize-1);
end

% Set ticks, axis labels, stage labels and legend text to bold before export.
set(findall(fig, '-property', 'FontWeight'), 'FontWeight', 'bold');
drawnow;

%% 6. Export vector PDF, 600-dpi PNG, editable FIG and unscaled data.
if ~exist(cfg.OutputDir,'dir')
    [ok,msg] = mkdir(cfg.OutputDir);
    assert(ok,'无法建立输出文件夹：%s',msg);
end
base = fullfile(cfg.OutputDir,cfg.FileStem);
save([base '_data.mat'],'t','f_Hz','cfg','stageStart','stageRatio', ...
    'dataSource','scheduleSource');
savefig(fig,[base '.fig']);
drawnow;
if exist('exportgraphics','file') ~= 0
    exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,[base '.png'],'Resolution',600,'BackgroundColor','white');
else
    set(fig,'PaperUnits','centimeters','PaperPositionMode','manual', ...
        'PaperPosition',[0 0 cfg.FigureSizeCM],'PaperSize',cfg.FigureSizeCM, ...
        'InvertHardcopy','off');
    print(fig,[base '.pdf'],'-dpdf','-painters');
    print(fig,[base '.png'],'-dpng','-r600');
end
fprintf('已导出到：%s\n',cfg.OutputDir);
yShown = f_Hz(t>=xLim(1) & t<=xLim(2),cfg.Channels);
if any(yShown(:)<cfg.YLim(1) | yShown(:)>cfg.YLim(2))
    fprintf(['提示：部分曲线超出所选纵轴范围，仅在图面外被裁剪；', ...
        'FIG和MAT仍保存完整数据，未限幅。\n']);
end
