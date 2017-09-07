function tomo = display_tomo(tomo_data)
%
% USEAGE: TOMO = display_tomo(TOMO_DATA)
%
% INPUT ARGUMENTS:
%
% TOMO_DATA:
%  A 1xm struct array containg the tomograph parameters. Fields include:
%       FIELD       DESCRIPTION
%       filename    The ring-parameter filename.
%
%       shift       The axial position of the ring.
%
%       rot         The transaxial rotation, a counterclockwise rotation
%                   between 0 and 360 degrees.
%
% OUTPUTS:
%  A 1xm cell array where m is the number of rings in the tomograph. Each entry
%  of the 1xm cell array contains a 1xn cell array where n is the number of
%  block-types in a ring. Each entry of the 1xn cell array contains an array of
%  structs for each block of a block-type.
%  The structs contain the following fields:
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
% Deniz Alpay, 2017-08-27

tomo = cell(1, numel(tomo_data));

% A 3x2 matrix with the x, y, z bounds of the tomograph. The first column
% is the lower bound and the second column is the upper bound. The three rows
% correspond to the x, y, and z dimensions.
lims = zeros(3, 2);

% Loop through each ring in the tomograph
for i = 1:numel(tomo_data)

    % The data describing the rings of the tomograph is retrieved in one of
    % two ways. First, whether a ring visualization file exists is checked.
    % The ring visualization file is the output of display_ring, which is
    % automatically saved when display_ring is called. If display_ring was
    % called on the rings used in the tomograph, the ring data will be aquired
    % from the ring visualization files. In the case that the ring visualization
    % file doesn't exist, whether the ring parameter file exists is checked.
    % If the ring parameter file exists, it is read and parsed using read_ring.
    % The second option is much more computationally expensive.

    ring_parms = tomo_data(i).filename;
    ring_fn = split(ring_parms , '.');
    ring_vis = [ char(ring_fn(1)), '.mat' ];

    if (exist(ring_vis, 'file') == 2)
        load(ring_vis);
    elseif (exist(ring_parms, 'file') == 2)
        [ring, blocks] = read_ring(ring_parms);
        ring = build_ring({ring, blocks});
    else
        error([ 'Neither the ring-visualization file "%s" or the ', ...
                'ring-parameter file "%s" can be found.' ], ...
                ring_vis, ring_parms);
    end

    % Loop through each block-type of the ring
    tomo{i} = cell(1, numel(ring));
    for j = 1:numel(ring)
        tomo{i}{j} = ring{j};
        % Loop through each block of a block-type
        for k = 1:numel(ring{j})

            % Shift the vertices in the z-dimension
            shift = ones(8, 1)*tomo_data(i).shift;
            tomo{i}{j}(k).Vertices(:, 3) = ring{j}(k).Vertices(:, 3) + shift;

            % Rotate the blocks in the xy-plane
            ang = tomo_data(i).rot*pi/180;
            rot = [ cos(ang), -sin(ang);
                    sin(ang), cos(ang)  ];
            tomo{i}{j}(k).Vertices(:, 1:2) = ring{j}(k).Vertices(:, 1:2)*rot;

            % Get the coordinates of the bounds of the ring, and replace
            % the existing bounds if the new bounds contain the old.
            % Loop through the x, y, z dimensions
            for n = 1:3
                verts = tomo{i}{j}(k).Vertices;
                lims(n, 1) = min(lims(n, 1), min(verts(:, n)));
                lims(n, 2) = max(lims(n, 2), max(verts(:, n)));
            end
        end
    end
end

% Display the tomograph
view(3);
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
xlim(lims(1, :));
ylim(lims(2, :));
zlim(lims(3, :));
pbaspect((lims(:, 2) - lims(:, 1))');

% Loop through each ring in the tomograph
hold on;
for i = 1:numel(tomo)
    % Loop through each block-type in the ring
    for j = 1:numel(tomo{i})
        % Loop through each block of a block-type
        for k = 1:numel(tomo{i}{j})
            patch(tomo{i}{j}(k));
        end
    end
end
hold off;

end
