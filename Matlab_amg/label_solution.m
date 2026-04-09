function label_solution(sz)
% LABEL_SOLUTION  Read p_interior_<sz>.txt and add node indices.
%
%   label_solution('100x1000x1000')
%
%   Reads the BC file to recover interior node indices, then writes
%   p_interior_<sz>_labeled.txt with columns: [node_index, pressure].

if nargin < 1, sz = '100x100x10'; end
if iscell(sz), sz = sz{1}; end

dname = 'four_large_Cartesian_networks';

%% Read BC to get boundary nodes and total node count
bcfile = fullfile(dname, ['nwk_Cartesian_' sz '.BC']);
fprintf('Reading %s ...\n', bcfile);
BC = load(bcfile);
boundary_nodes = BC(:,1);

%% Total node count from pMx file (number of lines)
pmxfile = fullfile(dname, ['nwk_Cartesian_' sz '.pMx']);
fid = fopen(pmxfile, 'r');
np = 0;
while ~feof(fid)
    fgetl(fid);
    np = np + 1;
end
fclose(fid);
fprintf('  np = %d, n_boundary = %d\n', np, length(boundary_nodes));

%% Interior nodes (same logic as read_nwk.m)
interior_nodes = setdiff((1:np)', boundary_nodes);
n_int = length(interior_nodes);
fprintf('  n_interior = %d\n', n_int);

%% Read pressure values
solfile = sprintf('p_interior_%s.txt', sz);
fprintf('Reading %s ...\n', solfile);
p = load(solfile);
assert(length(p) == n_int, 'Mismatch: %d values vs %d interior nodes', length(p), n_int);

%% Write labeled output
outfile = sprintf('p_interior_%s_labeled.txt', sz);
fprintf('Writing %s ...\n', outfile);
fid = fopen(outfile, 'w');
fprintf(fid, '%d %.16e\n', [interior_nodes(:)'; p(:)']);
fclose(fid);
fprintf('Done. %d lines written.\n', n_int);

end
