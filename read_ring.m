function [ring, blocks] = read_ring(filename)
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
% RING
%  Structure of ring parameters. Fields include:
%      FIELD       DESCRIPTION
%      filename    Ring-parameter filename (overwritten)
%
%      r_min       minimum radial extent of bounding volume.
%
%      r_max       maximum radial extent of bounding volume.
%
%      z_min       minimum axial position of bounding volume.
%
%      z_max       maximum axial position of bounding volume.
%
% BLOCKS
%  Array of structures giving number and arrangement for each
%  type of block. Fields of the I-th block type are stored in
%  BLOCKS(I).FIELD, where FIELD is:
%      FIELD       DESCRIPTION
%      filename    Name of block parameter file
%
%      r_min       Inner radius at which to distribute blocks
%
%      count       Scalar set to one since the the block-parameter
%                  doesn't organize each block by block-type.
%
%      azimuth     Azmuthal position of one or all blocks with
%                  respect to block reference point (given in
%                  block param file).  If only one angle is given
%                  then this value is take as azmuthal position of
%                  first block.  In this case, the rest of the
%                  (count-1) blocks are uniformly distributed about
%                  this position.  An array of azmuthal positions
%                  for all blocks is returned in this field in
%                  either case.
%      z_pos       Axial position blocks (one value for all blocks).
%
%      tilt        Transaxial orientation of blocks from block
%                  reference point.  A scalar value indicates that
%                  all blocks have the same tilt.  A vector must
%                  used to provide different tilt for each block.
%
% Deniz Alpay, 2017-09-06

if (exist(filename, 'file') == 2)
    ring_params = fileread(filename);
else
    error('The ring-parameter file "%s" cannot be found.', filename);
end

% RING

% filename
ring.filename = filename;

% r_min and r_max
dims = ['x', 'y'];
innerstr = {'REAL\s+ring_', '_inner_radius\s+=\s+'};
outerstr = {'REAL\s+ring_', '_outer_radius\s+=\s+'};

for i = 1:numel(dims)
    innerpat = [innerstr{1}, dims(i), innerstr{2}];
    inner(i) = get_val(innerpat, ring_params);

    outerpat = [outerstr{1}, dims(i), outerstr{2}];
    outer(i) = get_val(outerpat, ring_params);
end

ring.r_min = min(inner);
ring.r_max = max(outer);

% z_min and z_max
minstr = 'REAL\s+ring_z_minimum\s+=\s+';
maxstr = 'REAL\s+ring_z_maximum\s+=\s+';

ring.z_min = get_val(minstr, ring_params);
ring.z_max = get_val(maxstr, ring_params);

% BLOCKS

% filename
fnexpr = 'STR\s+ring_block_parameter_file\s*=\s*"([a-zA-Z_-.]+)"';
fn = regexp(ring_params, fnexpr, 'tokens');

% r_min
rstr = 'REAL\s+ring_block_radial_position\s*=\s*';
r_min = get_val(rstr, ring_params);

% azimuth
azstr = 'REAL\s+ring_block_angular_position\s*=\s*';
azimuth = get_val(azstr, ring_params);

% z_pos
zstr = 'REAL\s+ring_block_z_position\s*=\s*';
z_pos = get_val(zstr, ring_params);

% tilt
tiltstr = 'REAL\s+ring_block_transaxial_orientation\s*=\s*';
tilt = get_val(tiltstr, ring_params);

% Check that the there's the same number of each parameter
num_blks = numel(fn);
if ( num_blks ~= [numel(r_min), numel(azimuth), numel(z_pos), numel(tilt)])
    error('The parameter file is not formatted correctly.');
end

for i = 1:num_blks
    blocks(i) = struct( 'filename', fn{i},      ...
                        'r_min',    r_min(i),   ...
                        'count',    1,          ...
                        'azimuth',  azimuth(i), ...
                        'z_pos',    z_pos(i),   ...
                        'tilt',     tilt(i) );
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
