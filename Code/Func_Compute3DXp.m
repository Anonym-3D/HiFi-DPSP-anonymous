% This function triangulates 3D points from a main-camera pixel grid and the
% corresponding horizontal coordinate in a second camera or projector.
%
% Input parameters:
% mCamera1 -- 3-by-4 projection matrix of the main camera
% mCamera2 -- 3-by-4 projection matrix of the auxiliary camera or projector
% mXp      -- matched horizontal coordinate in the auxiliary camera or projector, format: H-by-W
% dX       -- optional horizontal offset added to the main-camera pixel grid
% dY       -- optional vertical offset added to the main-camera pixel grid
%
% Output:
% mX       -- reconstructed X coordinate
% mY       -- reconstructed Y coordinate
% mZ       -- reconstructed Z coordinate
%
% Author: Anonymous Authors, 2026/06/01
function [mX, mY, mZ] = Func_Compute3DXp(mCamera1, mCamera2, mXp, dX, dY)
if nargin < 4 || isempty(dX)
    dX = 0;
end
if nargin < 5 || isempty(dY)
    dY = 0;
end

[H, W] = size(mXp);
[mx, my] = meshgrid(1:W, 1:H);
mx = mx + dX;
my = my + dY;

% Build a three-equation linear system for each pixel:
% camera x, camera y, and auxiliary/projector x.
a11 = mCamera1(1, 1) - mx.*mCamera1(3, 1);
a12 = mCamera1(1, 2) - mx.*mCamera1(3, 2);
a13 = mCamera1(1, 3) - mx.*mCamera1(3, 3);
b1  = mx.*mCamera1(3, 4) - mCamera1(1, 4);

a21 = mCamera1(2, 1) - my.*mCamera1(3, 1);
a22 = mCamera1(2, 2) - my.*mCamera1(3, 2);
a23 = mCamera1(2, 3) - my.*mCamera1(3, 3);
b2  = my.*mCamera1(3, 4) - mCamera1(2, 4);

a31 = mCamera2(1, 1) - mXp.*mCamera2(3, 1);
a32 = mCamera2(1, 2) - mXp.*mCamera2(3, 2);
a33 = mCamera2(1, 3) - mXp.*mCamera2(3, 3);
b3  = mXp.*mCamera2(3, 4) - mCamera2(1, 4);

mDet = a11.*(a22.*a33 - a23.*a32) - ...
       a12.*(a21.*a33 - a23.*a31) + ...
       a13.*(a21.*a32 - a22.*a31);

mX = (b1 .*(a22.*a33 - a23.*a32) - ...
      a12.*(b2 .*a33 - a23.*b3 ) + ...
      a13.*(b2 .*a32 - a22.*b3 )) ./ mDet;
mY = (a11.*(b2 .*a33 - a23.*b3 ) - ...
      b1 .*(a21.*a33 - a23.*a31) + ...
      a13.*(a21.*b3  - b2 .*a31)) ./ mDet;
mZ = (a11.*(a22.*b3  - b2 .*a32) - ...
      a12.*(a21.*b3  - b2 .*a31) + ...
      b1 .*(a21.*a32 - a22.*a31)) ./ mDet;

mInvalid = abs(mDet) < eps | ~isfinite(mXp);
mX(mInvalid) = nan;
mY(mInvalid) = nan;
mZ(mInvalid) = nan;
end
