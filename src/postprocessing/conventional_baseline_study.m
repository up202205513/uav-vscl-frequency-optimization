%% conventional_baseline_study.m
%
% Sweeps conventional fixed-angle laminates as a reference baseline for the
% four-gene VSCL results. Setting all four genes equal (T0 = T1 = theta) gives
% a uniform-angle laminate, so this is a strict subset of the GA search space.
%
% Run BEFORE the main GA to populate the baseline table (~15-25 min, 13 solves).
%
% Outputs (results/):
%   baseline_results_<timestamp>.json
%   results/plots/conventional_baseline_<timestamp>.pdf/.png

clc;

%% Path Setup
thisScript = mfilename('fullpath');
srcPostDir = fileparts(thisScript);
srcDir     = fileparts(srcPostDir);
projectDir = fileparts(srcDir);

addpath(fullfile(srcDir, 'optimization'));
addpath(fullfile(srcDir, 'preprocessing'));
addpath(srcPostDir);

resultsDir = fullfile(projectDir, 'results');
plotsDir   = fullfile(resultsDir, 'plots');
if ~isfolder(plotsDir), mkdir(plotsDir); end

%% Sweep
panelHalfWidth = 200;

fprintf('Conventional Fixed-Angle Baseline Study\n');
fprintf('=========================================\n');
fprintf('%-10s  %-10s  %-25s  %s\n', 'Theta [deg]', 'Freq [Hz]', 'Effective laminate', 'Note');
fprintf('%s\n', repmat('-', 1, 70));

angles = -90:15:90;
freqs  = NaN(size(angles));

for k = 1:numel(angles)
    ang = angles(k);

    chromosome = [ang, ang, ang, ang];
    freqs(k)   = evaluate_panel(chromosome, panelHalfWidth);

    ply1 = ang;
    ply2 = mod(ang + 90, 180);
    laminateStr = sprintf('[%d / %d / %d / %d]', ply1, ply2, ply2, ply1);

    note = '';
    if ang == 0,  note = '<-- original baseline'; end
    if ang == 90, note = 'equivalent to [90/0/0/90]'; end

    fprintf('%-11.0f  %-10.4f  %-25s  %s\n', ang, freqs(k), laminateStr, note);
end

%% Summary
[bestConvFreq, bestIdx] = max(freqs);
baselineFreq  = 39.5443;
vsclOptFreq   = 41.2425;   % two-zone (4-gene) VSCL optimum <90|35> [Hz]

fprintf('\nBest conventional laminate: theta = %.0f deg  ->  %.4f Hz\n', ...
        angles(bestIdx), bestConvFreq);
fprintf('Original [0/90/90/0] baseline:          %.4f Hz\n', baselineFreq);
fprintf('Two-zone VSCL optimum <90|35>:          %.4f Hz\n', vsclOptFreq);

%% Plot
fig = figure('Name', 'Conventional Baseline Study', 'NumberTitle', 'off');
plot(angles, freqs, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
yline(baselineFreq, 'b--', 'Baseline [0/90/90/0]',  'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left');
yline(vsclOptFreq,  'r--', 'Two-zone VSCL optimum', 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left');

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11);
xlabel('Uniform Ply Rotation Angle \theta [deg]', 'Color', 'k', 'FontSize', 12);
ylabel('Mode 1 Natural Frequency [Hz]',           'Color', 'k', 'FontSize', 12);
t = title('Mode 1 Frequency vs. Constant Fiber Angle - Conventional Laminates');
t.Color = 'k'; t.FontSize = 13;
set(gca, 'GridLineStyle', ':', 'GridColor', [0.2 0.2 0.2], 'GridAlpha', 0.2);
xticks(angles); grid on;

%% Save Results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% JSON results record
runResult.timestamp             = datestr(now, 'yyyy-mm-dd HH:MM:SS');
runResult.angles_deg            = angles;
runResult.frequencies_Hz        = freqs;
runResult.best.angle_deg        = angles(bestIdx);
runResult.best.frequency_Hz     = bestConvFreq;
runResult.reference.baseline_Hz = baselineFreq;
runResult.reference.vscl2g_Hz   = vsclOptFreq;

jsonFile = fullfile(resultsDir, sprintf('baseline_results_%s.json', timestamp));
fid = fopen(jsonFile, 'w');
fprintf(fid, '%s\n', jsonencode(runResult));
fclose(fid);
fprintf('Results saved : %s\n', jsonFile);

exportgraphics(fig, fullfile(plotsDir, sprintf('conventional_baseline_%s.pdf', timestamp)), 'ContentType', 'vector');
exportgraphics(fig, fullfile(plotsDir, sprintf('conventional_baseline_%s.png', timestamp)), 'Resolution', 300);
fprintf('Figures saved : %s\n', plotsDir);
