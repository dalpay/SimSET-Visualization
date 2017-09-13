function ring = display_ring(ring_data, block_data)
%
% USAGE: RING = display_ring(DATA);
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

ring_fn = split(ring_data.filename , '.');
ring_vis = [ char(ring_fn(1)), '.mat' ];

if (exist(ring_vis, 'file') == 2)
    load(ring_vis);
elseif (exist(ring_data.filename, 'file') == 2)
    ring = build_ring(ring_data, block_data);
else
    error([ 'Neither the ring-visualization file "%s" or the ', ...
            'ring-parameter file "%s" can be found.' ], ...
            ring_vis, ring_data.filename);
end

% Define the bounding radial cyclinders
n = 50; % Number of faces on the cylinders
[X_inner, Y_inner, ~] = cylinder(ring_data.r_min, n);
[X_outer, Y_outer, ~] = cylinder(ring_data.r_max, n);
Z_inner = [ring_data.z_min*ones(1, n+1); ring_data.z_max*ones(1, n+1)];
Z_outer = [ring_data.z_min*ones(1, n+1); ring_data.z_max*ones(1, n+1)];

% Display the rings and cyclinders
view(3);
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
xlim([-ring_data.r_max, ring_data.r_max]);
ylim([-ring_data.r_max, ring_data.r_max]);
zlim([ring_data.z_min, ring_data.z_max]);
pbaspect([2*ring_data.r_max*ones(1, 2), (ring_data.z_max - ring_data.z_min)]);

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
