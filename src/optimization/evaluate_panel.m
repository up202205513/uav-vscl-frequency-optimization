function fitness = evaluate_panel(chromosome, panelHalfWidth) %#ok<INUSD>
% evaluate_panel  GA fitness evaluator: builds the steered deck, runs
% CalculiX, and returns the Mode-1 frequency (-1.0 on failure).
%
% Inputs:
%   chromosome     - [T0_A, T1_A, T0_B, T1_B] gene vector [deg]
%   panelHalfWidth - Panel half-width [mm] (reserved for future extensions)
% Output:
%   fitness - Mode 1 natural frequency [Hz], or -1.0 on evaluation failure

    % Resolve project root from this file's location
    thisFile   = mfilename('fullpath');
    srcOptDir  = fileparts(thisFile);           % .../src/optimization
    srcDir     = fileparts(srcOptDir);          % .../src
    projectDir = fileparts(srcDir);             % project root

    meshFile   = fullfile(projectDir, 'fea', 'mesh', 'Baseline_Trial.inp');
    runDir     = fullfile(projectDir, 'results', 'current_run');
    solverDir  = fullfile(projectDir, 'tools', 'PrePoMax v2.5.1 dev', ...
                          'PrePoMax v2.5.1 dev', 'Solver');

    if ~isfolder(runDir), mkdir(runDir); end

    outputJobName = 'Steered_Trial';
    inpFileName   = fullfile(runDir, [outputJobName, '.inp']);
    datFileName   = fullfile(runDir, [outputJobName, '.dat']);

    % Purge legacy solver output files to prevent stale readings
    extensionsToPurge = {'.dat', '.inp', '.frd', '.sta', '.cvg', '.log'};
    for extIdx = 1:length(extensionsToPurge)
        targetFile = fullfile(runDir, [outputJobName, extensionsToPurge{extIdx}]);
        if isfile(targetFile)
            try delete(targetFile); catch; end
        end
    end

    T0_A = chromosome(1);  T1_A = chromosome(2);
    T0_B = chromosome(3);  T1_B = chromosome(4);

    originalPath = getenv('PATH');
    origDir      = pwd();

    try
        % Step 1: Generate the two-zone steered finite element input deck
        build_vscl_input_deck(meshFile, inpFileName, T0_A, T1_A, T0_B, T1_B);

        % Step 2: Invoke the CalculiX solver from the run directory so that
        %         all solver output files are written there automatically
        cd(runDir);
        setenv('PATH', [originalPath pathsep solverDir]);
        [status, solverLog] = system(sprintf('ccx_dynamic %s', outputJobName));
        setenv('PATH', originalPath);
        cd(origDir);

        if status ~= 0
            warning('evaluate_panel: CalculiX returned status %d for [%.0f,%.0f,%.0f,%.0f].\n%s', ...
                    status, T0_A, T1_A, T0_B, T1_B, solverLog);
            fitness = -1.0;
            return;
        end

        % Step 3: Parse Mode 1 natural frequency from the solver output
        if isfile(datFileName)
            fitness = extract_calculix_frequency(datFileName);
        else
            warning('evaluate_panel: Solver finished but "%s" was not created.', datFileName);
            fitness = -1.0;
        end

    catch ME
        setenv('PATH', originalPath);
        if ~strcmp(pwd(), origDir), cd(origDir); end
        warning('evaluate_panel: Exception for [%.0f,%.0f,%.0f,%.0f] - %s', ...
                T0_A, T1_A, T0_B, T1_B, ME.message);
        fitness = -1.0;
    end
end
