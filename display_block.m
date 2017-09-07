function block = display_block(blk_data)
%
% USAGE: BLOCK = display_block(blk_data);
%
% INPUT ARGUMENTS:
%
% BLK_DATA
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

x = [ blk_data.min(1), blk_data.x, blk_data.max(1) ];

block = cell(1, numel(x) - 1);

cm = colormap(lines);

% prev, curr, and next are true if the previous, current, next layer is
% nonscintillating and false if they are scintillating
prev = logical(1);
curr = logical(1);
next = isscalar(blk_data.act{1}) && (blk_data.act{1} == 0);

% Loop through each layer of the block
for i = 1:numel(x) - 1

    faces = (length(blk_data.y{i})+1)*(length(blk_data.z{i})+1);

    block{i} = struct(  'XData', [], ...
                        'YData', [], ...
                        'ZData', [], ...
                        'FaceVertexCData', [], ...
                        'EdgeColor', 'none', ...
                        'FaceColor', 'flat');

    % Update prev, curr, next
    prev = curr;
    curr = next;
    if (i < numel(x) - 1)
        next = isscalar(blk_data.act{i+1}) && (blk_data.act{i+1} == 0);
    else
        next = logical(1);
    end

    % The current layer is transparent if it is nonscintillating
    if (curr)
        block{i}.FaceAlpha = 0.3;
    end

    % Append the boundaries to the y and z vectors
    y = [ blk_data.min(2), blk_data.y{i}, blk_data.max(2) ];
    z = [ blk_data.min(3), blk_data.z{i}, blk_data.max(3) ];

    % Turn the idx vector into a matrix where the indices of the rows and
    % columns correspond to the y and z coordinates of that segment
    idx = reshape(blk_data.idx{i}, length(y)-1, length(z)-1);

    % If current layer is nonscintillating and both the previous and next
    % layers are scintillating the front and back faces of the current layer
    % aren't generated.
    % Using DeMorgan's law ~(curr && ~prev && ~next) = ~curr || prev || next
    if (~curr || prev || next)
        % Define the yz planar cross-section of the layer
        ind = 1;
        for j = 1:length(y) - 1
            for k = 1:length(z) - 1
                block{i}.YData(:, ind) = [ y(j+1) y(j+1) y(j) y(j) ]';
                block{i}.ZData(:, ind) = [ z(k+1) z(k) z(k) z(k+1) ]';
                block{i}.FaceVertexCData(ind, :) = cm(idx(j, k), :);
                ind = ind + 1;
            end
        end
    end

    % Generate the front and back layers only if the current layer is
    % scintillating or the current layer is nonscintillating and previous and
    % next layers are also nonscintillating.
    % By DeMorgan's law,
    % ~(curr || (curr && prev && next)) = curr && ~(prev && next)
    if (curr && (~prev || ~next))
        if (~prev)          % The previous layer is scintillating
            % Only display the back face
            block{i}.XData = x(i+1)*ones(4, faces);
        elseif (~next)     % The next layer is scintillating
            % Only display the front face
            block{i}.XData = x(i)*ones(4, faces);
        end
    else
        % Display both the front and back faces
        block{i}.XData = [ x(i)*ones(4, faces), x(i+1)*ones(4, faces) ];
        block{i}.YData = repmat(block{i}(1).YData, 1, 2);
        block{i}.ZData = repmat(block{i}(1).ZData, 1, 2);
        block{i}.FaceVertexCData = repmat(block{i}.FaceVertexCData, 2, 1);
    end

    % Connect the front and back with four faces
    x_sides = [ x(i)*ones(2, 4); x(i+1)*ones(2, 4) ];

    y_sides = [ [ y(end), y(end), y(end), y(end) ]', ...
                [ y(end), y(1),   y(1),   y(end) ]', ...
                [ y(1),   y(1),   y(1),   y(1)   ]', ...
                [ y(end), y(1),   y(1),   y(end) ]' ];

    z_sides = [ [ z(end), z(1),   z(1),   z(end) ]', ...
                [ z(1),   z(1),   z(1),   z(1)   ]', ...
                [ z(1),   z(end), z(end), z(1)   ]', ...
                [ z(end), z(end), z(end), z(end) ]' ];

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
xlim([blk_data.min(1), blk_data.max(1)]);
ylim([blk_data.min(2), blk_data.max(2)]);
zlim([blk_data.min(3), blk_data.max(3)]);
pbaspect([	blk_data.max(1) - blk_data.min(1),
			blk_data.max(2) - blk_data.min(2),
			blk_data.max(3) - blk_data.min(3)	]);

hold on;
% Loop through each layer of the block
for i = 1:numel(block)
    patch(block{i});
end
hold off;

end
