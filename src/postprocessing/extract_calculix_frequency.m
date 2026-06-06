function frequency = extract_calculix_frequency(datFileName)
% extract_calculix_frequency  Parse the Mode-1 frequency from CalculiX output.
%
% Reads the EIGENVALUE OUTPUT table in a .dat file and returns column 4
% (frequency in Hz) of the first data row, validated as finite and positive.
%
% Inputs:
%   datFileName - Path to the CalculiX output file (*.dat)
% Outputs:
%   frequency   - Mode 1 natural frequency [Hz]

    fid = fopen(datFileName, 'r');
    if fid == -1
        error('extract_calculix_frequency: Failed to open "%s".', datFileName);
    end

    frequency         = NaN;
    inEigenvalueBlock = false;

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end

        if contains(line, 'E I G E N V A L U E')
            inEigenvalueBlock = true;
            continue;
        end

        if inEigenvalueBlock
            tokens = strsplit(line);

            if numel(tokens) >= 4 && strcmp(tokens{1}, '1')
                candidate = str2double(tokens{4});
                if ~isnan(candidate) && candidate > 0
                    frequency = candidate;
                    break;
                end
            end
        end
    end

    fclose(fid);

    if isnan(frequency)
        error(['extract_calculix_frequency: Mode 1 frequency not found in "%s". ' ...
               'Verify solver execution status and confirm the EIGENVALUE OUTPUT ' ...
               'section is present in the .dat file.'], datFileName);
    end
end
