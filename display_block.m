function block = display_block(blk)
%
% USAGE: display_block(BLK);
%
% INPUT ARGUMENTS:
%
% BLK
%  Structure of block parameters. Fields include:
%      FIELD       DESCRIPTION
%      fn    Name of block parameter file
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
%
% OUTPUT: BLOCK
%  A 1xm cell array containing scalar structs, where m is the number of
%  layers in the block.
%  The stuctures contain the following fields:
%      XData       A 4xn matrix where n is two times the number of faces for
%                  the current layer plus the four faces connecting the front
%                  to the back. Each column describes the x-coordinates of a
%                  vertex, with four coordinates/rows to describe a rectangle.
%                  The last four columns connect the front face to the back.
%                  If the layer is nonscintillating, it is a 4x6 matrix.
%                  If the layer is scintillating, the number of faces depends
%                  on the complexity of the block.
%
%      YData       Same as XData for the y-coordinates.
%
%      ZData       Same as XData for the z-coordinates.
%
%      FaceVertexCData
%                  A nx3 matrix where n is the two times the number of faces
%                  plus the four faces connecting the front to the back.
%                  Each column describes the RGB value of the face that is
%                  described the columns of XData, YData, ZData that has the
%                  same index as the row of FaceVertexCData.
%
%      EdgeColor   A character array containing 'none'.
%
%      FaceColor   A character array containing 'flat'.
%
% Deniz Alpay, 2017-08-15

x = [ blk.min(1), blk.x, blk.max(1) ];

block = cell(1, numel(x) - 1);

cm = colormap(lines);

prev = logical(0);
curr = logical(0);

% Loop through each layer
for i = 1:numel(x) - 1

    faces = (length(blk.y{i})+1)*(length(blk.z{i})+1);

    block{i} = struct(  'XData', [], ...
                        'YData', [], ...
                        'ZData', [], ...
                        'FaceVertexCData', [], ...
                        'EdgeColor', 'none', ...
                        'FaceColor', 'flat');

    % prev is true if the current layer is nonscintillating
    prev = curr;
    % curr is true if the current layer is nonscintillating
    curr = isscalar(blk.act{i}) && (blk.act{i} == 0);

    % Nonscintillating layers are transparent
    if (curr)
        block{i}.FaceAlpha = 0.3;
    end

    y = [ blk.min(2), blk.y{i}, blk.max(2) ];
    z = [ blk.min(3), blk.z{i}, blk.max(3) ];
    idx = reshape(blk.idx{i}, length(y)-1, length(z)-1);

    % Define the cross-section of the layer
    ind = 1;
    for j = 1:length(y) - 1
        for k = 1:length(z) - 1
            block{i}.YData(:, ind) = [ y(j+1) y(j+1) y(j) y(j) ]';
            block{i}.ZData(:, ind) = [ z(k+1) z(k) z(k) z(k+1) ]';
            block{i}.FaceVertexCData(ind, :) = cm(idx(j, k), :);
            ind = ind + 1;
        end
    end

    % To prevent interference between the faces of the adjacent layers whether
    % the front or back of a nonscintillating layer is displayed depends on
    % the previous and next layer.
    if ( curr )             % The current layer is nonscintillating
        if ( ~prev )        % The previous layer is scintillating
            % Only display the back face
            block{i}.XData = x(i+1)*ones(4, faces);
        %{
        else if ( ~next )   % The next layer is scintillating
            % Only display the front face
            block{i}.XData = x(i)*ones(4, faces);
        else                % Both the front and back

        %}
        end
    else                    % The current layer is scintillating
        % Display both the front and back faces
        block{i}.XData = [ x(i)*ones(4, faces), x(i+1)*ones(4, faces) ];
        block{i}.YData = repmat(block{i}(1).YData, 1, 2);
        block{i}.ZData = repmat(block{i}(1).ZData, 1, 2);
        block{i}.FaceVertexCData = repmat(block{i}.FaceVertexCData, 2, 1);
    end

    % Connect the front and back with four faces
    x_sides = [ x(i)*ones(2, 4); x(i+1)*ones(2, 4) ];

    y_sides = [ y(end)*ones(4, 1), ...
                [ y(end), y(1), y(1), y(end) ]', ...
                y(1)*ones(4, 1), ...
                [ y(end), y(1), y(1), y(end) ]' ];

    z_sides = [ [ z(end), z(1), z(1), z(end) ]', ...
                z(1)*ones(4, 1), ...
                [ z(1), z(end), z(end), z(1) ]', ...
                z(end)*ones(4, 1) ];

    c_sides = repmat(cm(idx(1, 1), :), 4, 1);

    block{i}.XData = [ block{i}.XData, x_sides ];
    block{i}.YData = [ block{i}.YData, y_sides ];
    block{i}.ZData = [ block{i}.ZData, z_sides ];
    block{i}.FaceVertexCData = [ block{i}.FaceVertexCData; c_sides ];
end

% Display the block detector
view(3);
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
xlim([blk.min(1) blk.max(1)]);
ylim([blk.min(2) blk.max(2)]);
zlim([blk.min(3) blk.max(3)]);
pbaspect([	blk.max(1) - blk.min(1),
			blk.max(2) - blk.min(2),
			blk.max(3) - blk.min(3)	]);

hold on;
% Loop through each layer
for i = 1:numel(block)
    patch(block{i});
end
hold off;

end
