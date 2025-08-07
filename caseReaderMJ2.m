function [faceMx, ptCoordMx, dia, BC, np, nf, nt] = caseReaderMJ2(filename)
    % Read face matrix
    faceMx = load(strcat(filename, '.fMx'));
    
    % Read point coordinates
    ptCoordMx = load(strcat(filename, '.pMx'));
    
    % Read diameters
    dia = load(strcat(filename, '.dia'));
    
    % Optional: Read boundary conditions
    bcFile = strcat(filename, '.BC');
    if isfile(bcFile)
        BC = load(bcFile);
    else
        BC = []; % Return empty if not present
        warning('No .BC file found. Returning empty BC matrix.');
    end
    
    % Basic stats
    np = length(ptCoordMx); 
    nf = size(faceMx, 1); 
    nt = np + nf;
end