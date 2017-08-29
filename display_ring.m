function ring = display_ring(ring_data, block_data)
%
% USEAGE: BLOCKS = display_ring(RING_DATA, BLOCKS_DATA);
%
% INPUT ARGUMENTS:
%
% RING_DATA
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
% BLOCK_DATA
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
% Deniz Alpay, 2017-07-30

ring = cell(1, size(block_data, 2));
cm = colormap(lines);

% Loop through each block type
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

            dims = ['x', 'y', 'z'];
            refstr = {'REAL	block_reference_', '       =   '};
            minstr = {'REAL	block_', '_minimum         =   '};
            maxstr = {'REAL	block_', '_maximum         =   '};
            token = '(-?[0-9]+\.[0-9]+)';

            ref = zeros(1, 3);
            min = zeros(1, 3);
            max = zeros(1, 3);

            % Loop through the x, y, z dimensions
            for j = 1:3
                % Reference
                ref_expr = [refstr{1}, dims(j), refstr{2}, token];
                ref_token = regexp(blockparams, ref_expr, 'tokens');
                ref(j) = str2double(ref_token{1, 1});

                % Minimum
                min_expr = [minstr{1}, dims(j), minstr{2}, token];
                min_token = regexp(blockparams, min_expr, 'tokens');
                min(j) = str2double(min_token{1, 1});

                % Maxiumum
                max_expr = [maxstr{1}, dims(j), maxstr{2}, token];
                max_token = regexp(blockparams, max_expr, 'tokens');
                max(j) = str2double(max_token{1, 1});
            end
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

    % Loop through each individual block
    for j = 1:block.count

        ring{i}(j) = struct('Vertices', [], 'Faces', [], 'FaceColor', color);

        % Vertices are defined from the top-left going clockwise when
        % looking from the xy-plane which is parallel to the xz plane of
        % the coordinate space where the block dimensions are defined
        vert = [    min(1) max(3);
                    max(1) max(3);
                    max(1) min(3);
                    min(1) min(3)   ];

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

% Display the blocks and the bounding radial cyclinders
n = 50; % Number of faces on the cylinders
[X_inner, Y_inner, ~] = cylinder(ring_data.r_min, n);
[X_outer, Y_outer, ~] = cylinder(ring_data.r_max, n);
Z_inner = [ring_data.z_min*ones(1, n+1); ring_data.z_max*ones(1, n+1)];
Z_outer = [ring_data.z_min*ones(1, n+1); ring_data.z_max*ones(1, n+1)];

view(3);
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
xlim([-ring_data.r_max ring_data.r_max]);
ylim([-ring_data.r_max ring_data.r_max]);
zlim([-ring_data.r_max ring_data.r_max]);
pbaspect(2*ring_data.r_max*ones(1, 3));

hold on;
outer = surf(X_outer, Y_outer, Z_outer);
inner = surf(X_inner, Y_inner, Z_inner);
set(outer, 'FaceAlpha', 0.5, 'FaceColor', 'yellow', 'EdgeColor', 'none');
set(inner, 'FaceAlpha', 0.5, 'FaceColor', 'yellow', 'EdgeColor', 'none');

% Loop through each block-type of the ring
for i = 1:numel(ring)
    % Loop through each block of a block-type
    for j = 1:numel(ring{i})
       patch(ring{i}(j));
    end
end
hold off;

end
