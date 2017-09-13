function block = read_block(filename)
%
% USAGE: RING = read_ring(FILENAME)
%
% INPUT ARGUMENTS:
%
% FILENAME:
%  The filename of the ring-parameter file.
%
% OUTPUTS:
%
% BLOCK
%  Structure of block parameters. Fields include:
%      FIELD       DESCRIPTION
%      fn          Name of block parameter file
%
%      ref         3D coordinate of block reference point used to
%                  position block in ring
%
%      min         3-element vector of lower bound on block volume (X,Y,Z)
%                  given with respect to "ref" point
%
%      max         3-element vector of upper bound on block volume (X,Y,Z)
%                  given with respect to "ref" point
%
%      x           Array of X coordinates for layer divisions (w.r.t. ref).
%                  For Nx layers, "x" should be a (Nx-1)-element vector.
%                  Values should be monotonically increasing.
%
%      y           List of 1D-arrays (one array for each layer) of the form:
%                     {[...],[...],...}, where [...] is a 1D-array.
%                  The array for a given layer gives the y coordinates for
%                  voxel boundaries/material changes. For (Ny+1) segments,
%                  there are Ny divisions.  Note that the divisions for
%                  each layer can be distinct.
%
%      z           List of 1D-arrays (one array for each layer) of the form:
%                     {[...],[...],...}, where [...] is a 1D-array.
%                  The array for a given layer gives the z coordinates for
%                  voxel boundaries/material changes. For Nz segments,
%                  there are (Nz-1) divisions.  Note that the divisions for
%                  each layer can be distinct.
%
%      idx         List of 2D-arrays (one array for each layer) of the form:
%                     {[...],[...],...}, where [...] is a ((Ny+1)*(Nz+1))
%                  2D-array. The array for a given layer gives the material
%                  index for each of the ((Ny+1)*(Nz+1)) elements of that
%                  layer voxel boundaries/material changes.
%
%      act         List of 2D-arrays (one array for each layer) of the form:
%                     {[...],[...],...}, where [...] is a ((Ny+1)*(Nz+1))
%                  2D-array.  The array for a given layer flags if each of
%                  the ((Ny+1)*(Nz+1)) elements of that layer is active
%                  (0=inactive, 1=active).

if (exist(filename, 'file') == 2)
    blk_params = fileread(filename);
else
    error('The block-parameter file "%s" cannot be found.', filename);
end

block.fn = filename;

dims = ['x', 'y', 'z'];
ref_str = {'REAL\s+block_reference_', '\s+=\s+'};
min_str = {'REAL\s+block_', '_minimum\s+=\s+'};
max_str = {'REAL\s+block_', '_maximum\s+=\s+'};

block.ref = get_coordinates(ref_str, blk_params);
block.min = get_coordinates(min_str, blk_params);
block.max = get_coordinates(max_str, blk_params);



end

function coor = get_coordinates(str, params)
% Parses params, the block-parameter file, for the given pattern which
% defines the coordinates of a point of interest of the block.
%
% USAGE: COOR = get_coordinates(PATTERN, PARAMS)
%
% INPUT ARGUMENTS:
%
% PATTERN
%  A 1x2 cell-array where the first entry is the line before the
%  character x, y, z and the second entry is the rest of the line.
%
% PARAMS
%  A character vector of the block-parameter file.
%
% OUTPUT: COOR
%  Returns a 1x3 vector with the x, y, z coordinates.

dims = ['x', 'y', 'z'];
coor = zeros(1, 3);

% Loop through the x, y, z dimensions
for i = 1:length(dims)
    expr = [str{1}, dims(i), str{2}];
    coor(i) = get_val(expr, params);
end

end

function val = get_val(pattern, params)
% Parses params for the given pattern which is a line which ends with a value.
%
% USAGE: VAL = get_val(PATTERN, PARAMS)
%
% INPUT ARGUMENTS:
%
% PATTERN
%  A character array which contains the beginning of the line which ends
%  with the value which is desired.
%
% PARAMS
%  A character vector of the block-parameter file.
%
% OUTPUT: VAL
%  A 1xm vector with the values embedded in the pattern, where m is the number
%  of lines that match in the pattern in params.

expr = [pattern, '(-?[0-9]+\.?[0-9]*)'];
out = regexp(params, expr, 'tokens');
if (~isempty(out))
    for m = 1:numel(out)
        val(m) = str2double(out{m});
    end
else
    error([ 'The parameter file does not contain the value ', ...
            'in the proper format.']);
end

end
