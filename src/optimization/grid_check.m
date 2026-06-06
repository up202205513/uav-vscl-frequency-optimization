%% grid_check.m  -  Grid Search Utility (does not overwrite GA results)
%
% Manual parameter sweep over T0/T1 combinations that pass the AFP
% manufacturing constraint, printing the Mode 1 frequency for each.

clc;

%% Path Setup
thisScript = mfilename('fullpath');
srcOptDir  = fileparts(thisScript);
srcDir     = fileparts(srcOptDir);
projectDir = fileparts(srcDir);

addpath(srcOptDir);
addpath(fullfile(srcDir, 'preprocessing'));
addpath(fullfile(srcDir, 'postprocessing'));

%% Parameters
panelHalfWidth = 200;
K_limit        = 3.28;  %#ok<NASGU>

T0_range = -45:15:45;
T1_range = -45:15:45;

for T0_test = T0_range
    for T1_test = T1_range
        [viol, ~] = check_manufacturing_constraint(T0_test, T1_test, panelHalfWidth);
        if ~viol
            f = evaluate_panel([T0_test, T1_test], panelHalfWidth);
            fprintf('T0 = %4.0f  T1 = %4.0f  ->  %.4f Hz\n', T0_test, T1_test, f);
        end
    end
end
