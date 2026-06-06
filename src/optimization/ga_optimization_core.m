%% ga_optimization_core.m  -  Multi-start driver for two-zone VSCL optimisation
%
% Runs a multi-start GA campaign (numRuns independent restarts from logged
% seeds, see run_ga.m) to maximise the Mode-1 natural frequency of a two-zone
% VSCL panel, and reports the global best.
%
%   chromosome = [T0_A, T1_A, T0_B, T1_B]
%   Zone A (x <  x_c): theta_A(x) = T0_A + ((T1_A - T0_A)/b)*(x_c - x)
%   Zone B (x >= x_c): theta_B(x) = T0_B + ((T1_B - T0_B)/b)*(x - x_c)
%
% Outputs (results/):
%   ga_results_<timestamp>.json                 - config, seeds, per-run + best
%   results/plots/ga_multistart_<timestamp>.*   - convergence + diversity plots

clc; clear; close all;

%% 0. Path Setup
thisScript = mfilename('fullpath');
srcOptDir  = fileparts(thisScript);          % .../src/optimization
srcDir     = fileparts(srcOptDir);           % .../src
projectDir = fileparts(srcDir);              % project root

addpath(srcOptDir);
addpath(fullfile(srcDir, 'preprocessing'));
addpath(fullfile(srcDir, 'postprocessing'));

resultsDir = fullfile(projectDir, 'results');
plotsDir   = fullfile(resultsDir, 'plots');
if ~isfolder(plotsDir), mkdir(plotsDir); end

%% 1. Algorithm Configuration
config.popSize           = 40;
config.numGen            = 50;
config.mutationRateStart = 0.30;   % high exploration early
config.mutationRateEnd   = 0.05;   % fine exploitation late
config.panelHalfWidth    = 200;    % [mm]
config.K_limit           = 3.28;   % AFP curvature manufacturing limit [1/m]
config.stallLimit        = 15;     % per-run early stop
config.allowedAngles     = -90:5:90;
config.numGenes          = 4;      % [T0_A, T1_A, T0_B, T1_B]

% Multi-start settings. Set numRuns = 1 for a single fixed-seed run (quick
% check); raise it for a proper campaign. Seeds are logged for reproducibility.
numRuns  = 5;
baseSeed = 12345;
seeds    = baseSeed + (0:numRuns - 1);

baselineFreq = 39.5443;   % [0/90/90/0] reference [Hz]

%% 2. Shared Fitness Cache (handle object - persists across all restarts)
cache = containers.Map('KeyType', 'char', 'ValueType', 'double');

fprintf('Two-Zone VSCL Optimisation - Multi-Start 4-Gene GA\n');
fprintf('Restarts: %d   Population: %d   Generations: %d\n', ...
        numRuns, config.popSize, config.numGen);
fprintf('Base seed: %d   Seeds: [%s]\n', baseSeed, strtrim(sprintf('%d ', seeds)));
fprintf('Mutation rate schedule: %.2f -> %.2f\n\n', ...
        config.mutationRateStart, config.mutationRateEnd);

%% 3. Multi-Start Loop
bestFreqAll  = NaN(1, numRuns);
bestGenAll   = NaN(1, numRuns);
bestChromAll = NaN(numRuns, config.numGenes);
fitHistAll   = cell(1, numRuns);
divHistAll   = cell(1, numRuns);

campaignTimer = tic;
for r = 1:numRuns
    fprintf('===== Restart %d/%d  (seed = %d) =====\n', r, numRuns, seeds(r));
    res = run_ga(config, seeds(r), cache);

    bestFreqAll(r)    = res.bestFreq;
    bestGenAll(r)     = res.bestGen;
    bestChromAll(r,:) = res.bestChromosome;
    fitHistAll{r}     = res.fitnessHistory;
    divHistAll{r}     = res.diversityHistory;

    fprintf('  -> restart best: %.4f Hz  [%3.0f %3.0f %3.0f %3.0f]\n\n', ...
            res.bestFreq, res.bestChromosome);
end
elapsed = toc(campaignTimer);

%% 4. Aggregate Across Restarts
if all(isnan(bestFreqAll))
    fprintf('--- Campaign Complete ---\n');
    fprintf('No restart produced a feasible, valid design.\n');
    fprintf('Check the solver path in evaluate_panel.m or the constraint settings.\n');
    return;
end

[overallBestFreq, bestRunIdx] = max(bestFreqAll);   % max() ignores NaN
overallBestChrom = bestChromAll(bestRunIdx, :);
improvementPct   = 100 * (overallBestFreq - baselineFreq) / baselineFreq;
validBest        = bestFreqAll(~isnan(bestFreqAll));   % NaN-free for statistics

fprintf('--- Campaign Complete (%.1f s) ---\n', elapsed);
fprintf('Restart best frequencies [Hz]:\n');
for r = 1:numRuns
    marker = ''; if r == bestRunIdx, marker = '  <-- global best'; end
    fprintf('  seed %d : %.4f Hz%s\n', seeds(r), bestFreqAll(r), marker);
end
fprintf('\nGlobal best : %.4f Hz  (seed %d, generation %d)\n', ...
        overallBestFreq, seeds(bestRunIdx), bestGenAll(bestRunIdx));
fprintf('Chromosome  : T0_A=%3.0f  T1_A=%3.0f  T0_B=%3.0f  T1_B=%3.0f\n', overallBestChrom);
fprintf('Mean +/- std: %.4f +/- %.4f Hz   (min %.4f, max %.4f)\n', ...
        mean(validBest), std(validBest), min(validBest), max(validBest));
fprintf('Baseline    : %.4f Hz   Improvement: +%.2f%%\n', baselineFreq, improvementPct);
fprintf('Unique CalculiX evaluations (cache size): %d\n', cache.Count);

%% 5. Plots - convergence and diversity across restarts
fig = figure('Name', 'GA Multi-Start - Two-Zone VSCL', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [100 100 1150 460]);

% (a) Convergence
subplot(1, 2, 1); hold on;
for r = 1:numRuns
    fh = fitHistAll{r};  vg = find(~isnan(fh));
    if r == bestRunIdx
        plot(vg, fh(vg), '-o', 'Color', [0 0.2 0.8], 'LineWidth', 2, 'MarkerSize', 4, ...
             'DisplayName', sprintf('best (seed %d)', seeds(r)));
    else
        plot(vg, fh(vg), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1, ...
             'HandleVisibility', 'off');
    end
end
yline(baselineFreq, 'r--', 'Baseline [0/90/90/0]', 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11, ...
    'GridLineStyle', ':', 'GridColor', [0.2 0.2 0.2], 'GridAlpha', 0.2);
xlabel('Generation', 'Color', 'k', 'FontSize', 12);
ylabel('Best Frequency (Mode 1) [Hz]', 'Color', 'k', 'FontSize', 12);
title('Convergence Across Restarts', 'Color', 'k', 'FontSize', 13);
legend('Location', 'southeast'); grid on; hold off;

% (b) Diversity
subplot(1, 2, 2); hold on;
for r = 1:numRuns
    dh = divHistAll{r};  vg = find(~isnan(dh));
    if r == bestRunIdx
        plot(vg, dh(vg), '-o', 'Color', [0 0.2 0.8], 'LineWidth', 2, 'MarkerSize', 4);
    else
        plot(vg, dh(vg), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1);
    end
end
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11, ...
    'GridLineStyle', ':', 'GridColor', [0.2 0.2 0.2], 'GridAlpha', 0.2);
xlabel('Generation', 'Color', 'k', 'FontSize', 12);
ylabel('Unique Individuals in Population', 'Color', 'k', 'FontSize', 12);
title('Population Diversity', 'Color', 'k', 'FontSize', 13);
ylim([0 config.popSize + 1]); grid on; hold off;

%% 6. Save Results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

out.timestamp                  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
out.config                     = config;
out.config.numRuns             = numRuns;
out.config.baseSeed            = baseSeed;
out.seeds                      = seeds;
out.elapsed_seconds            = elapsed;
out.perRun.seeds               = seeds;
out.perRun.bestFreq_Hz         = bestFreqAll;
out.perRun.bestGeneration      = bestGenAll;
out.perRun.bestChromosome      = bestChromAll;
out.overallBest.frequency_Hz   = overallBestFreq;
out.overallBest.chromosome     = overallBestChrom;
out.overallBest.seed           = seeds(bestRunIdx);
out.overallBest.run            = bestRunIdx;
out.overallBest.generation     = bestGenAll(bestRunIdx);
out.stats.mean_Hz              = mean(validBest);
out.stats.std_Hz               = std(validBest);
out.stats.min_Hz               = min(validBest);
out.stats.max_Hz               = max(validBest);
out.totalUniqueEvaluations     = cache.Count;
out.baseline_frequency_Hz      = baselineFreq;
out.improvement_pct            = improvementPct;
out.fitnessHistories           = fitHistAll;
out.diversityHistories         = divHistAll;

jsonFile = fullfile(resultsDir, sprintf('ga_results_%s.json', timestamp));
fid = fopen(jsonFile, 'w');
fprintf(fid, '%s\n', jsonencode(out));
fclose(fid);
fprintf('\nResults saved : %s\n', jsonFile);

exportgraphics(fig, fullfile(plotsDir, sprintf('ga_multistart_%s.pdf', timestamp)), 'ContentType', 'vector');
exportgraphics(fig, fullfile(plotsDir, sprintf('ga_multistart_%s.png', timestamp)), 'Resolution', 300);
fprintf('Figures saved : %s\n', plotsDir);
