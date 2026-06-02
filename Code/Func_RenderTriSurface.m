% This function renders a reconstructed depth map as a shaded triangular
% surface. It builds a structured image-grid triangulation and removes invalid
% faces whose vertices contain NaN values or whose longest 3D edge is too long.
%
% Input parameters:
% hAxes          -- target axes handle
% mX             -- reconstructed X coordinate map, format: H-by-W
% mY             -- reconstructed Y coordinate map, format: H-by-W
% mZ             -- reconstructed Z coordinate map, format: H-by-W
% vClim          -- color-axis range, format: [cMin, cMax]
% vXlim          -- X-axis range, format: [xMin, xMax]
% vYlim          -- Y-axis range, format: [yMin, yMax]
% vZlim          -- Z-axis range, format: [zMin, zMax]
% fMaxEdgeLength -- maximum valid triangle edge length in 3D space
%
% Output:
% hSurface       -- handle of the rendered trisurf object
%
% Author: Anonymous Authors, 2026/06/01
function hSurface = Func_RenderTriSurface(hAxes, mX, mY, mZ, vClim, vXlim, vYlim, vZlim, fMaxEdgeLength)
if nargin < 9 || isempty(fMaxEdgeLength)
    fMaxEdgeLength = 10;
end

[H, W] = size(mZ);
mVertexIndex = reshape(1:H*W, H, W);
mTopLeft = mVertexIndex(1:end-1, 1:end-1);
mTopRight = mVertexIndex(1:end-1, 2:end);
mBottomLeft = mVertexIndex(2:end, 1:end-1);
mBottomRight = mVertexIndex(2:end, 2:end);
mFaces = [mTopLeft(:), mTopRight(:), mBottomLeft(:);
          mBottomLeft(:), mTopRight(:), mBottomRight(:)];

mXYZ = [mX(:), mY(:), mZ(:)];
mValidVertex = all(isfinite(mXYZ), 2);
mValidFace = all(mValidVertex(mFaces), 2);

if any(mValidFace)
    mFacesCandidate = mFaces(mValidFace, :);
    mP1 = mXYZ(mFacesCandidate(:, 1), :);
    mP2 = mXYZ(mFacesCandidate(:, 2), :);
    mP3 = mXYZ(mFacesCandidate(:, 3), :);

    vEdge12 = sqrt(sum((mP1 - mP2).^2, 2));
    vEdge23 = sqrt(sum((mP2 - mP3).^2, 2));
    vEdge31 = sqrt(sum((mP3 - mP1).^2, 2));
    mFaces = mFacesCandidate(max([vEdge12, vEdge23, vEdge31], [], 2) <= fMaxEdgeLength, :);
else
    mFaces = zeros(0, 3);
end

axes(hAxes);
cla(hAxes);
hSurface = trisurf(mFaces, mXYZ(:, 1), mXYZ(:, 2), mXYZ(:, 3), mXYZ(:, 3), ...
    'Parent', hAxes, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.35, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.1);

axis(hAxes, 'equal');
xlim(hAxes, vXlim);
ylim(hAxes, vYlim);
zlim(hAxes, vZlim);
view(hAxes, 0, -90);
set(hAxes, 'Box', 'off', ...
    'XColor', 'black', ...
    'YColor', 'black', ...
    'ZColor', 'black', ...
    'TickDir', 'out');
colormap(hAxes, jet);
caxis(hAxes, vClim);
% material(hAxes, 'dull');
camlight;
end
