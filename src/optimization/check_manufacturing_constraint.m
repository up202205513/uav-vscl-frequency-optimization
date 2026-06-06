function [isViolated, maxCurvature] = check_manufacturing_constraint(T0, T1, panelHalfWidth)
% check_manufacturing_constraint  Analytical AFP curvature pre-filter.
%
% O(1) check that rejects unmanufacturable steering before any FEA call.
% For a linear law the peak curvature is kappa = (d_theta/dx)*cos^2(theta),
% maximised at the panel centre (or at a zero-crossing if T0, T1 differ in
% sign). The limit K_limit = 3.28 m^-1 is the ~305 mm (12 in) AFP minimum
% bend radius, below which tow deposition wrinkles.
%
% Inputs:
%   T0             - Fiber angle at panel centre [deg]
%   T1             - Fiber angle at panel edge   [deg]
%   panelHalfWidth - Distance from panel centre to edge [mm]
% Outputs:
%   isViolated   - true if the curvature limit is exceeded
%   maxCurvature - Peak curvature of the steered fiber path [1/m]

    K_limit = 3.28;   % AFP minimum-radius manufacturing limit [1/m]

    T0_rad      = deg2rad(T0);
    T1_rad      = deg2rad(T1);
    halfWidth_m = panelHalfWidth / 1000;   % Convert mm to m

    angleGradient = abs(T1_rad - T0_rad) / halfWidth_m;   % [rad/m]

    if T0 * T1 <= 0
        cosAtPeak = 1.0;
    else
        cosAtPeak = cos(min(abs(T0_rad), abs(T1_rad)));
    end

    maxCurvature = angleGradient * cosAtPeak^2;   % [1/m]
    isViolated   = (maxCurvature > K_limit);
end
