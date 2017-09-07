function ring = build_ring(data)
%
% USEAGE: RING = build_ring(DATA);
%
% INPUT ARGUMENTS:
%
% DATA
%  If DATA is a 1x2 struct it is assumed that DATA{1} contains the ring
%  parameters and DATA{2} contains the block arrangement.
%  If DATA is a 1x3 struct, (as is the case when using the output of
%  define_ring) DATA{1} is ignored. Instead, DATA{2} is taken to be
%  the ring parameters and DATA{3} is taken to be the block arrangement.
%
% DATA{1}
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
% DATA{2}
%  Array of structures giving number and arrangement for each
%  type of block. Fields of the I-th block type are stored in
%  BLOCKS(I).FIELD, where FIELD is:
%      FIELD       DESCRIPTION
%      filename    Name of block parameter file
%
%      r_min       Inner radius at which to distribute blocks
%
%      count       Number of this block type in ring
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
%
%      z_pos       Axial position blocks (one value for all blocks).
%
%      tilt        Transaxial orientation of blocks from block
%                  reference point.  A scalar value indicates that
%                  all blocks have the same tilt.  A vector must
%                  used to provide different tilt for each block.
%
% OUTPUTS: RING
%  A 1xm cell array containg arrays of strutures for each block of a block
%  type, where m is the number of block types. Each element of the cell
%  array contains a 1xn array of structs for that block type, where n
%  is the number of blocks of that block type.
%  The structures contains the following fields:
%       FIELD       DESCRIPTION
%       Vertices    A 8x3 matrix where the three columns correspond to
%                   the x, y, and z coordinates of the 8 vertices that
%                   define a block.
%
%       Faces       A 6x4 matrix where each row defines a face of the block
%                   and the values correspond to the rows in the
%                   Vertices field. The four columns connect four
%                   vertices to define a rectangular face.
%
%       FaceColor   A 3x1 vector containing the RGB values of the faces.
%
% Deniz Alpay, 2017-09-06

elems = size(data);
if (elems(1) == 1 && elems(2) == 3)
    ring_data = data{2};
    block_data = data{3};
elseif (elems(1) == 1 && elems(2) == 2)
    ring_data = data{1};
    block_data = data{2};
else
    error('The input should be a 1x2 or a 1x3 cell array. See header comment.');
end

ring = cell(1, size(block_data, 2));

cm = colormap(lines);

% Loop through each block type of the ring
prev_filename = '';
for i = 1:numel(block_data)

    block = block_data(i);

    % Read the block parameter file for the current block type and retrieve
    % the ref, min, and max coordinates using regular expressions, if the
    % current block has a different block-parameter file from the previous.
    if (exist(block.filename, 'file') == 2)
        if (~strcmp(block.filename, prev_filename))

            color = cm(i, :);

            prev_filename = block.filename;

            blockparams = fileread(block.filename);

            refstr = {'REAL\s+block_reference_', '\s+=\s+'};
            minstr = {'REAL\s+block_', '_minimum\s+=\s+'};
            maxstr = {'REAL\s+block_', '_maximum\s+=\s+'};

            ref = get_coordinates(refstr, blockparams);
            min = get_coordinates(minstr, blockparams);
            max = get_coordinates(maxstr, blockparams);
        end
    else
        error(  'The block-parameter file "%s" cannot be found.', ...
                block.filename);
    end

    % Compute the remaining azimuths if they are evenly distributed
    if (isscalar(block.azimuth))
        block.azimuth = block.azimuth + (0:(block.count-1))*(360/block.count);
    end

    % Repeat the tilt for the remaining blocks if they're the same
    if (isscalar(block.tilt))
        block.tilt = repmat(block.tilt, block.count, 1);
    end

    % Loop through each block of a block-type
    for j = 1:block.count

        ring{i}(j) = struct('Vertices', [], 'Faces', [], 'FaceColor', color);

        % Vertices are defined from the top-left going clockwise when
        % looking from the xy-plane which is parallel to the xz plane of
        % the coordinate space where the block dimensions are defined
        vert = [    min(1), max(3);
                    max(1), max(3);
                    max(1), min(3);
                    min(1), min(3)   ];

        % Rotate the block by the sum of the azimuth and tilt by
        % multiplying the vertices in the xy-plane by the rotation matrix
        az = block.azimuth(j)*pi/180;
        tilt = block.tilt(j)*pi/180;
        ang = az + tilt;
        rot = [ cos(ang), -sin(ang);
                sin(ang), cos(ang)  ];
        vert = (rot*vert')';

        % Append the z coordinates to the vertices
        vert = [    vert, min(2)*ones(4, 1);
                    vert, max(2)*ones(4, 1) ];

        % Shift the position of the vertices by the azimuth and radius in
        % the xy-plane and the z-position in the z dimension
        [x_pos, y_pos] = pol2cart(az, block.r_min);
        scale = [ ones(8, 1)*x_pos, ones(8, 1)*y_pos, ones(8, 1)*block.z_pos ];
        vert = vert - ones(8, 3)*ref' + scale;

        ring{i}(j).Vertices = vert;
        ring{i}(j).Faces = [    1 2 3 4;    % Right
                                4 3 7 8;    % Bottom
                                5 6 7 8;    % Left
                                5 6 2 1;    % Top
                                5 1 4 8;    % Front
                                6 2 3 7 ];  % Back
    end
end

% Save the ring visualization data for the tomograph visualization
filename = split(ring_data.filename, '.');
filename = [ char(filename(1)), '.mat' ];
save(filename, 'ring');

end

function coor = get_coordinates(pattern, params)
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
token = '(-?[0-9]+\.?[0-9]*)';
dims = ['x', 'y', 'z'];
coor = zeros(1, 3);

% Loop through the x, y, z dimensions
for i = 1:3
    expr = [pattern{1}, dims(i), pattern{2}, token];
    val = regexp(params, expr, 'tokens');
    if (~isempty(val))
        coor(i) = str2double(val{1, 1});
    else
        error([ 'The block-parameter file does not contain the ', ...
                'block dimensions in the proper format.']);
    end
end

end
