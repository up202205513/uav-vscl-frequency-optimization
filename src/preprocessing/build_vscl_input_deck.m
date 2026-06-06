function build_vscl_input_deck(inputFileName, outputFileName, T0_A, T1_A, T0_B, T1_B)
% build_vscl_input_deck  Two-zone steered CalculiX input-deck generator.
%
% Reads a baseline unsteered mesh and applies an independent linear steering
% law to each chordwise zone, writing per-element *Orientation / *Shell
% Section cards (CalculiX RECTANGULAR + type-3 rotation, PrePoMax convention).
%
%   Zone A (x <  x_c): theta_A(x) = T0_A + ((T1_A - T0_A)/b)*(x_c - x)
%   Zone B (x >= x_c): theta_B(x) = T0_B + ((T1_B - T0_B)/b)*(x - x_c)
%   b = panel half-width; T0 = centre angle, T1 = edge angle (zones independent).
%
% Inputs:
%   inputFileName  - Path to the baseline CalculiX input deck (.inp)
%   outputFileName - Path for the generated steered input deck (.inp)
%   T0_A, T1_A     - Zone A centre / edge angles [deg]
%   T0_B, T1_B     - Zone B centre / edge angles [deg]
%
% With no arguments, runs in standalone mode and writes a fiber-map figure
% to results/plots/.

    % Resolve project root for standalone default paths
    thisFile   = mfilename('fullpath');
    srcPrepDir = fileparts(thisFile);
    srcDir     = fileparts(srcPrepDir);
    projectDir = fileparts(srcDir);

    if nargin == 0
        inputFileName  = fullfile(projectDir, 'fea', 'mesh', 'Baseline_Trial.inp');
        outputFileName = fullfile(projectDir, 'results', 'current_run', 'Steered_Trial.inp');
        T0_A =  0;  T1_A = 45;
        T0_B = 30;  T1_B = -30;
        isStandalone = true;
    else
        isStandalone = false;
    end

    if ~isfile(inputFileName)
        error('Pre-processing Aborted: Input file "%s" not found.', inputFileName);
    end

    % Ensure output directory exists
    outDir = fileparts(outputFileName);
    if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end

    [plateLength, ~] = parseInpMesh(inputFileName);
    panelHalfWidth   = plateLength / 2;
    fileLines        = readlines(inputFileName);

    if isStandalone
        fprintf('File Read Complete: %d lines from "%s".\n', length(fileLines), inputFileName);
        fprintf('Panel length: %.1f mm   Half-width: %.1f mm\n\n', plateLength, panelHalfWidth);
    end

    % ------------------------------------------------------------------
    % Pass 1: Locate keyword block boundaries
    % ------------------------------------------------------------------
    nodeStartLine    = 0;  nodeEndLine    = 0;
    elementStartLine = 0;  elementEndLine = 0;
    shellSectionLine = 0;  shellSectionEnd = 0;

    for i = 1:length(fileLines)
        currentLine = strtrim(fileLines(i));

        if nodeStartLine == 0 && ...
                startsWith(currentLine, '*Node', 'IgnoreCase', true) && ...
                ~startsWith(currentLine, '*Nset', 'IgnoreCase', true)
            nodeStartLine = i + 1;
        end

        if elementStartLine == 0 && ...
                startsWith(currentLine, '*Element', 'IgnoreCase', true) && ...
                contains(currentLine, 'S8R', 'IgnoreCase', true)
            elementStartLine = i + 1;
            if nodeEndLine == 0, nodeEndLine = i - 1; end
        end

        if elementStartLine > 0 && elementEndLine == 0 && ...
                i > elementStartLine && startsWith(currentLine, '*')
            elementEndLine = i - 1;
        end

        if shellSectionLine == 0 && ...
                startsWith(currentLine, '*Shell Section', 'IgnoreCase', true)
            shellSectionLine = i;
        end

        if shellSectionLine > 0 && shellSectionEnd == 0 && ...
                i > shellSectionLine && startsWith(currentLine, '*')
            shellSectionEnd = i - 1;
        end
    end

    if nodeStartLine == 0 || nodeEndLine == 0
        error('build_vscl_input_deck: *Node block not found in "%s".', inputFileName);
    end
    if elementStartLine == 0 || elementEndLine == 0
        error('build_vscl_input_deck: *Element (S8R) block not found in "%s".', inputFileName);
    end
    if shellSectionLine == 0 || shellSectionEnd == 0
        error('build_vscl_input_deck: *Shell Section block not found in "%s".', inputFileName);
    end

    % ------------------------------------------------------------------
    % Pass 2: Build nodal coordinate lookup table
    % ------------------------------------------------------------------
    maxNodeID = 0;
    for k = nodeStartLine:nodeEndLine
        rawLine = strtrim(fileLines(k));
        if isempty(rawLine) || startsWith(rawLine, '*'), continue; end
        parsedRow = str2double(strsplit(rawLine, ','));
        if ~isnan(parsedRow(1)), maxNodeID = max(maxNodeID, parsedRow(1)); end
    end

    nodeLookupTable = NaN(maxNodeID, 3);
    parsedNodeCount = 0;
    for k = nodeStartLine:nodeEndLine
        rawLine = strtrim(fileLines(k));
        if isempty(rawLine) || startsWith(rawLine, '*'), continue; end
        parsedRow = str2double(strsplit(rawLine, ','));
        if isnan(parsedRow(1)), continue; end
        nodeLookupTable(parsedRow(1), :) = parsedRow(2:4);
        parsedNodeCount = parsedNodeCount + 1;
    end

    validX       = nodeLookupTable(~isnan(nodeLookupTable(:,1)), 1);
    panelCenterX = min(validX) + panelHalfWidth;

    % ------------------------------------------------------------------
    % Pass 3: Parse element connectivity, assign zone, evaluate steering law
    % ------------------------------------------------------------------
    elementIDs       = [];
    elementCentroids = [];
    elementAngles    = [];
    elementZones     = [];

    for k = elementStartLine:elementEndLine
        rawLine = strtrim(fileLines(k));
        if isempty(rawLine) || startsWith(rawLine, '*') || startsWith(rawLine, '**')
            continue;
        end
        parsedRow = str2double(strsplit(rawLine, ','));
        if isnan(parsedRow(1)) || numel(parsedRow) < 9, continue; end

        eID          = parsedRow(1);
        connectivity = parsedRow(2:9);
        nodalX       = nodeLookupTable(connectivity, 1);
        nodalY       = nodeLookupTable(connectivity, 2);
        Xc           = mean(nodalX);
        Yc           = mean(nodalY);

        if Xc < panelCenterX
            theta = T0_A + ((T1_A - T0_A) / panelHalfWidth) * (panelCenterX - Xc);
            zone  = 'A';
        else
            theta = T0_B + ((T1_B - T0_B) / panelHalfWidth) * (Xc - panelCenterX);
            zone  = 'B';
        end

        elementIDs       = [elementIDs;       eID];       %#ok<AGROW>
        elementCentroids = [elementCentroids; Xc, Yc];    %#ok<AGROW>
        elementAngles    = [elementAngles;    theta];      %#ok<AGROW>
        elementZones     = [elementZones;     zone];       %#ok<AGROW>
    end

    totalElements = numel(elementIDs);
    if totalElements == 0
        error('build_vscl_input_deck: No valid S8R element records could be parsed.');
    end

    % ------------------------------------------------------------------
    % Orientation field visualisation (standalone execution only)
    % ------------------------------------------------------------------
    if isStandalone
        plotsDir = fullfile(projectDir, 'results', 'plots');
        if ~isfolder(plotsDir), mkdir(plotsDir); end

        fig = figure('Name', 'Two-Zone Variable Stiffness Orientation Map', 'NumberTitle', 'off');
        scatter(elementCentroids(:,1), elementCentroids(:,2), 18, elementAngles, 'filled');
        colormap(jet);
        set(gcf, 'Color', 'w');
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        cb = colorbar; cb.Color = 'k';
        cb.Label.String = 'Ply Orientation Angle, \theta [deg]';
        cb.Label.Color  = 'k';
        xline(panelCenterX, 'w--', 'Zone boundary', 'LineWidth', 1.2, ...
              'LabelVerticalAlignment', 'bottom');
        xlabel('X Coordinate [mm]', 'Color', 'k');
        ylabel('Y Coordinate [mm]', 'Color', 'k');
        t = title(sprintf(['Two-Zone VSCL  |  Zone A: \\langle%.0f|%.0f\\rangle   ' ...
                           'Zone B: \\langle%.0f|%.0f\\rangle'], T0_A, T1_A, T0_B, T1_B));
        t.Color = 'k';
        set(gca, 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.2);
        axis equal tight; grid on;

        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        exportgraphics(fig, fullfile(plotsDir, sprintf('fiber_map_%s.pdf', timestamp)), 'ContentType', 'vector');
        exportgraphics(fig, fullfile(plotsDir, sprintf('fiber_map_%s.png', timestamp)), 'Resolution', 300);
        fprintf('Fiber map saved to: %s\n', plotsDir);
    end

    % ------------------------------------------------------------------
    % Ply Parser: Extract baseline stacking sequence
    % ------------------------------------------------------------------
    plyData = [];
    for compositeLine = (shellSectionLine + 1):shellSectionEnd
        rawPly = strtrim(fileLines(compositeLine));
        if isempty(rawPly) || startsWith(rawPly, '**'), continue; end
        parts = strsplit(rawPly, ',');
        if length(parts) >= 3
            thk    = strtrim(parts{1});
            oriStr = strtrim(parts{2});
            mat    = strtrim(parts{3});
            if strcmpi(oriStr, 'ORI_90')
                baseAng = 90.0;
            else
                baseAng = str2double(oriStr);
                if isnan(baseAng), baseAng = 0.0; end
            end
            plyData = [plyData; struct('thk', thk, 'baseAng', baseAng, 'mat', mat)]; %#ok<AGROW>
        end
    end

    % ------------------------------------------------------------------
    % File Writer
    % ------------------------------------------------------------------
    fid = fopen(outputFileName, 'w');
    if fid == -1
        error('File Write Error: Cannot open "%s" for writing.', outputFileName);
    end

    for lineIdx = 1:(shellSectionLine - 1)
        fprintf(fid, '%s\n', fileLines(lineIdx));
    end

    % Write per-element blocks: *Elset -> *Orientation -> *Shell Section
    for e = 1:totalElements
        eID    = elementIDs(e);
        eAngle = elementAngles(e);

        fprintf(fid, '*Elset, elset=Set-Elem-%d\n', eID);
        fprintf(fid, ' %d\n', eID);

        for p = 1:length(plyData)
            finalAng = plyData(p).baseAng + eAngle;
            fprintf(fid, '*Orientation, name=Ori-E%d-P%d, system=Rectangular\n', eID, p);
            fprintf(fid, ' 1.0, 0.0, 0.0, 0.0, 1.0, 0.0\n');
            fprintf(fid, ' 3, %.4f\n', finalAng);
        end

        fprintf(fid, '*Shell Section, elset=Set-Elem-%d, composite\n', eID);
        for p = 1:length(plyData)
            fprintf(fid, '%s, , %s, Ori-E%d-P%d\n', plyData(p).thk, plyData(p).mat, eID, p);
        end
    end

    % Transcribe remainder, neutralising any duplicate *Shell Section blocks
    skipMode = false;
    for lineIdx = (shellSectionEnd + 1):length(fileLines)
        currentLine = strtrim(fileLines(lineIdx));
        if startsWith(currentLine, '*Shell Section', 'IgnoreCase', true)
            skipMode = true;
            fprintf(fid, '** Neutralised duplicate block: %s\n', fileLines(lineIdx));
            continue;
        end
        if skipMode
            if startsWith(currentLine, '*') && ~startsWith(currentLine, '**')
                skipMode = false;
            else
                fprintf(fid, '** Neutralised duplicate layer: %s\n', fileLines(lineIdx));
                continue;
            end
        end
        fprintf(fid, '%s\n', fileLines(lineIdx));
    end

    fclose(fid);

    if isStandalone
        nZoneA = sum(elementZones == 'A');
        nZoneB = sum(elementZones == 'B');
        fprintf('\nExecution Summary:\n');
        fprintf('  Nodes resolved     : %d\n', parsedNodeCount);
        fprintf('  Elements total     : %d  (Zone A: %d, Zone B: %d)\n', ...
                totalElements, nZoneA, nZoneB);
        fprintf('  Angle range Zone A : [%.2f, %.2f] deg\n', ...
                min(elementAngles(elementZones == 'A')), max(elementAngles(elementZones == 'A')));
        fprintf('  Angle range Zone B : [%.2f, %.2f] deg\n', ...
                min(elementAngles(elementZones == 'B')), max(elementAngles(elementZones == 'B')));
        fprintf('  *Orientation cards : %d\n', totalElements * length(plyData));
    end
end
