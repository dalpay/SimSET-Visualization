function tomo = display_tomo(tomo_data)
%
% USEAGE: TOMO = display_tomo(TOMO_DATA)
%
% INPUT ARGUMENTS:
%
% TOMO_DATA:
%  A struct array containg the tomograph parameters. Fields include:
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

% Loop through each ring in the tomograph
for i = 1:numel(tomo_data)

    % The data describing the rings of the tomograph are retrieved in one of
    % two ways. First, whether a ring-visualization file is exists is checked.
    % The ring-visualization file is the output of display_ring, which is
    % automatically saved when display_ring is called. If display_ring was
    % called on the rings used in the tomograph, the ring data will be aquired
    % from the ring-visualization files. In a future version of display_tomo,
    % the ring data will be retrieved from the ring-paramter files if the
    % ring-visualization files do not exist.

    ring_parms = tomo_data(i).filename;
    ring_fn = split(ring_parms , '.');
    ring_vis = [ char(ring_fn(1)), '.mat' ];

    if (exist(ring_vis, 'file') == 2)
        load(ring_vis);
        ring
    elseif (exist(ring_parms, 'file') == 2)
        message = [ 'The tomograph cannot be visualized using the ',
                    'ring-parameter file at the moment'  ];
        disp(message);
        %{
        build_ring will be everything in display_ring except the visual part
        ring = build_ring(read_ring(ring_parms));
        %}
    else
        err_str = [ 'Neither the ring-visualization file "%s" ',
                    'nor the ring-parameter file "%s" can be found.'    ];
        error(err_str, ring_vis, ring_parms);
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

            % Rotate the blocks in the xy-plane by multiplying the vertices
            % by the rotation matrix
            ang = tomo_data(i).rot*pi/180;
            rot = [ cos(ang), -sin(ang);
                    sin(ang), cos(ang)  ];
            tomo{i}{j}(k).Vertices(:, 1:2) = ring{j}(k).Vertices(:, 1:2)*rot;
        end
    end
end

% Display the ring-visualization
view(3);
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
c = 25;     % TODO: Fix the plotbox limits and aspect ratio
xlim([-c c]);
ylim([-c c]);
zlim([-c c]);
pbaspect(2*c*[ones(1, 3)]);

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
