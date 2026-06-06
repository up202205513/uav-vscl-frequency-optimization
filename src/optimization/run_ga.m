function result = run_ga(config, seed, cache)
% run_ga  One seeded run of the two-zone VSCL genetic algorithm.
%
% Executes a complete GA optimisation for a given seed and returns a result
% struct. Called repeatedly by the multi-start driver (ga_optimization_core.m),
% sharing a fitness cache. Reproducible (rng(seed)); discrete index-space
% BLX-alpha crossover + power mutation on the 5-deg grid; mutation rate and
% strength annealed over the run; per-generation diversity logged.
%
% Inputs:
%   config - struct: popSize, numGen, mutationRateStart, mutationRateEnd,
%            panelHalfWidth, K_limit, stallLimit, allowedAngles, numGenes
%   seed   - integer RNG seed for this run
%   cache  - containers.Map shared across runs (char -> double), mutated in place
% Output:
%   result - struct: seed, bestFreq, bestChromosome, bestGen,
%            fitnessHistory, diversityHistory

    rng(seed);

    popSize        = config.popSize;
    numGen         = config.numGen;
    panelHalfWidth = config.panelHalfWidth;
    K_limit        = config.K_limit;
    stallLimit     = config.stallLimit;
    allowedAngles  = config.allowedAngles;
    numGenes       = config.numGenes;

    % ------------------------------------------------------------------
    % Feasibility-aware initialisation (rejection sampling)
    % Only ~12% of random 4-gene individuals satisfy the AFP curvature limit
    % on both zones, so each starting individual is resampled until feasible.
    % ------------------------------------------------------------------
    population = zeros(popSize, numGenes);
    maxDraws   = 2000;
    for i = 1:popSize
        candidate = allowedAngles(randi(numel(allowedAngles), 1, numGenes));
        for attempt = 1:maxDraws
            violA = check_manufacturing_constraint(candidate(1), candidate(2), panelHalfWidth);
            violB = check_manufacturing_constraint(candidate(3), candidate(4), panelHalfWidth);
            if ~violA && ~violB, break; end
            candidate = allowedAngles(randi(numel(allowedAngles), 1, numGenes));
        end
        population(i, :) = candidate;
    end

    fitnessHistory   = NaN(numGen, 1);
    diversityHistory = NaN(numGen, 1);
    bestChromosome   = NaN(1, numGenes);
    bestFreqSoFar    = -Inf;

    % ------------------------------------------------------------------
    % Main evolutionary loop
    % ------------------------------------------------------------------
    for gen = 1:numGen
        % Annealed schedules: high exploration early, exploitation late
        mutationRate = config.mutationRateEnd + ...
                       (config.mutationRateStart - config.mutationRateEnd) * (1 - gen / numGen);
        powerIndex   = max(2, 10 * (1 - gen / numGen));

        rawFitness = zeros(popSize, 1);
        violations = zeros(popSize, 1);

        % --- Phase 1: constraint screening + (cached) fitness evaluation ---
        for i = 1:popSize
            T0_A = population(i, 1);  T1_A = population(i, 2);
            T0_B = population(i, 3);  T1_B = population(i, 4);

            [violA, curvA] = check_manufacturing_constraint(T0_A, T1_A, panelHalfWidth);
            [violB, curvB] = check_manufacturing_constraint(T0_B, T1_B, panelHalfWidth);

            if violA || violB
                violations(i) = max(curvA - K_limit, 0) + max(curvB - K_limit, 0);
                rawFitness(i) = -1.0;
            else
                violations(i) = 0;
                key = sprintf('%d_%d_%d_%d', round(population(i, :)));
                if isKey(cache, key)
                    rawFitness(i) = cache(key);
                else
                    rawFitness(i) = evaluate_panel(population(i, :), panelHalfWidth);
                    cache(key)    = rawFitness(i);
                end
            end
        end

        % Diversity: number of distinct individuals in the population
        diversityHistory(gen) = size(unique(population, 'rows'), 1);

        % Record best feasible fitness / chromosome
        feasibleMask = (violations == 0) & (rawFitness > 0);
        if any(feasibleMask)
            [genBestFreq, genBestIdx] = max(rawFitness .* feasibleMask - 1e6 * ~feasibleMask);
            fitnessHistory(gen) = genBestFreq;
            if genBestFreq > bestFreqSoFar
                bestFreqSoFar  = genBestFreq;
                bestChromosome = population(genBestIdx, :);
            end
            fprintf('    Gen %2d/%d | best %8.4f Hz | diversity %2d/%d | mutRate %.3f\n', ...
                    gen, numGen, genBestFreq, diversityHistory(gen), popSize, mutationRate);
        else
            fprintf('    Gen %2d/%d | no feasible individual | diversity %2d/%d\n', ...
                    gen, numGen, diversityHistory(gen), popSize);
        end

        % --- Phase 2: feasibility-based binary tournament selection ---
        matingPool = zeros(popSize, numGenes);
        for i = 1:popSize
            idx1 = randi(popSize);
            idx2 = randi(popSize);
            while idx2 == idx1, idx2 = randi(popSize); end

            feas1 = (violations(idx1) == 0);
            feas2 = (violations(idx2) == 0);
            if feas1 && feas2
                winnerIdx = ifelse(rawFitness(idx1) >= rawFitness(idx2), idx1, idx2);
            elseif feas1
                winnerIdx = idx1;
            elseif feas2
                winnerIdx = idx2;
            else
                winnerIdx = ifelse(violations(idx1) <= violations(idx2), idx1, idx2);
            end
            matingPool(i, :) = population(winnerIdx, :);
        end

        % --- Phase 3: elitism, crossover, mutation ---
        nextGeneration = zeros(popSize, numGenes);
        if sum(violations == 0) >= 2
            eliteScores = rawFitness;
            eliteScores(violations > 0) = -Inf;
            [~, sortedIdx] = sort(eliteScores, 'descend');
        else
            [~, sortedIdx] = sort(violations, 'ascend');
        end
        nextGeneration(1:2, :) = population(sortedIdx(1:2), :);

        for i = 3:2:popSize
            parent1 = matingPool(i, :);
            parent2 = matingPool(min(i + 1, popSize), :);

            [offspring1, offspring2] = discrete_crossover(parent1, parent2, allowedAngles);
            offspring1 = discrete_power_mutation(offspring1, mutationRate, powerIndex, allowedAngles);
            offspring2 = discrete_power_mutation(offspring2, mutationRate, powerIndex, allowedAngles);

            nextGeneration(i, :) = offspring1;
            if i + 1 <= popSize
                nextGeneration(i + 1, :) = offspring2;
            end
        end
        population = nextGeneration;

        % --- Stall convergence check ---
        if gen >= stallLimit
            recentHistory = fitnessHistory(gen - stallLimit + 1 : gen);
            if all(~isnan(recentHistory)) && all(recentHistory == recentHistory(1))
                fprintf('    [stall] best unchanged for %d generations; stopping at gen %d\n', ...
                        stallLimit, gen);
                fitnessHistory   = fitnessHistory(1:gen);
                diversityHistory = diversityHistory(1:gen);
                break;
            end
        end
    end

    % ------------------------------------------------------------------
    % Assemble result
    % ------------------------------------------------------------------
    validGens = find(~isnan(fitnessHistory));
    if isempty(validGens)
        bestFreq = NaN;  bestGen = NaN;
    else
        [bestFreq, bestGen] = max(fitnessHistory);
    end

    result.seed             = seed;
    result.bestFreq         = bestFreq;
    result.bestChromosome   = bestChromosome;
    result.bestGen          = bestGen;
    result.fitnessHistory   = fitnessHistory(:)';
    result.diversityHistory = diversityHistory(:)';
end

%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------

function res = ifelse(cond, a, b)
    if cond, res = a; else, res = b; end
end

function [c1, c2] = discrete_crossover(p1, p2, allowedAngles)
% BLX-alpha blend crossover in gene-index space, so offspring land exactly
% on the discrete grid (no sub-grid rounding that would collapse diversity).
    n     = numel(allowedAngles);
    step  = allowedAngles(2) - allowedAngles(1);
    i1    = round((p1 - allowedAngles(1)) / step) + 1;
    i2    = round((p2 - allowedAngles(1)) / step) + 1;
    alpha = 0.5;

    c1 = zeros(size(p1));
    c2 = zeros(size(p2));
    for j = 1:numel(p1)
        lo = min(i1(j), i2(j));
        hi = max(i1(j), i2(j));
        d  = hi - lo;
        low  = lo - alpha * d;
        high = hi + alpha * d;

        ni1 = min(max(round(low + rand() * (high - low)), 1), n);
        ni2 = min(max(round(low + rand() * (high - low)), 1), n);
        c1(j) = allowedAngles(ni1);
        c2(j) = allowedAngles(ni2);
    end
end

function y = discrete_power_mutation(x, mutRate, powerIndex, allowedAngles)
% Power mutation in gene-index space; larger powerIndex => finer local moves,
% smaller => larger exploratory jumps across the grid.
    n    = numel(allowedAngles);
    step = allowedAngles(2) - allowedAngles(1);
    idx  = round((x - allowedAngles(1)) / step) + 1;
    y    = x;
    for j = 1:numel(x)
        if rand() < mutRate
            s = rand() ^ powerIndex;
            if rand() < 0.5
                newIdx = idx(j) - round(s * (idx(j) - 1));      % toward low edge
            else
                newIdx = idx(j) + round(s * (n - idx(j)));       % toward high edge
            end
            newIdx = min(max(newIdx, 1), n);
            y(j)   = allowedAngles(newIdx);
        end
    end
end
