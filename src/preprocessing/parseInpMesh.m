function [plateLength, totalElements] = parseInpMesh(filename)
% parseInpMesh  Mesh geometry extractor (keeps the pipeline mesh-independent).
%
% Single-pass scan of a CalculiX/Abaqus .inp deck for the X-span and the
% S8R element count, so no geometry is hardcoded.
%
% Inputs:
%   filename      - Path to the baseline CalculiX input deck (.inp)
% Outputs:
%   plateLength   - Panel span along X: max(X_nodes) - min(X_nodes) [mm]
%   totalElements - Number of valid S8R element records in the mesh

fid = fopen(filename, 'r');
if fid == -1
    error('parseInpMesh: Input deck "%s" could not be opened.', filename);
end

allX           = [];
elementCount   = 0;
inNodeBlock    = false;
inElementBlock = false;

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line), continue; end

    if startsWith(line, '*', 'IgnoreCase', false)
        if startsWith(line, '*Node', 'IgnoreCase', true) && ...
                ~startsWith(line, '*Nset', 'IgnoreCase', true)
            inNodeBlock    = true;
            inElementBlock = false;
        elseif startsWith(line, '*Element', 'IgnoreCase', true) && ...
                contains(line, 'S8R', 'IgnoreCase', true)
            inNodeBlock    = false;
            inElementBlock = true;
        else
            inNodeBlock    = false;
            inElementBlock = false;
        end
        continue;
    end

    if inNodeBlock
        parsedRow = str2double(strsplit(line, ','));
        if numel(parsedRow) >= 2 && ~isnan(parsedRow(2))
            allX(end+1) = parsedRow(2); %#ok<AGROW>
        end
    end

    % Count element records: first token is a positive integer element ID.
    % Handles S8R definitions that span continuation lines.
    if inElementBlock
        tokens     = strsplit(line, ',');
        firstToken = str2double(strtrim(tokens{1}));
        if ~isnan(firstToken) && firstToken == floor(firstToken) && firstToken > 0
            elementCount = elementCount + 1;
        end
    end
end
fclose(fid);

if isempty(allX)
    error('parseInpMesh: No nodal X-coordinate data found in "%s".', filename);
end
if elementCount == 0
    error('parseInpMesh: No valid S8R element records detected in "%s".', filename);
end

plateLength   = max(allX) - min(allX);
totalElements = elementCount;
end
