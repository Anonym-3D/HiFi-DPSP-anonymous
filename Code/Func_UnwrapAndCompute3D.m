% This function unwraps a wrapped phase map with a stereo-derived reference
% phase, then computes the corresponding 3D point cloud.
%
% Input parameters:
% mPhaseWrapLeft    -- wrapped phase map of the main camera
% mDisparity        -- disparity map from the auxiliary camera to the main camera
% mCamera1Rectified -- 3-by-4 projection matrix of the main camera
% mCamera2Rectified -- 3-by-4 projection matrix of the auxiliary camera
% mProjector        -- 3-by-4 projection matrix of the projector
% FSet              -- fringe frequency
%
% Output:
% mX,mY,mZ          -- reconstructed 3D point cloud
% mXRaw,mYRaw,mZRaw -- raw reconstructed 3D point cloud before optional rectification
% mPhase            -- unwrapped phase map of the main camera
%
% Author: Anonymous Authors, 2026/06/01
function [mX, mY, mZ, mXRaw, mYRaw, mZRaw, mPhase] = Func_UnwrapAndCompute3D(mPhaseWrapLeft, mDisparity, mCamera1Rectified, mCamera2Rectified, mProjector, FSet)
iProjectorWidth = 1280;
[iCameraHeight, iCameraWidth] = size(mPhaseWrapLeft);
[mx, ~] = meshgrid(1:iCameraWidth, 1:iCameraHeight);
mXcMatch = mDisparity + mx;

% Use stereo matching to obtain a reference depth for unwrapping.
[mXMatch, mYMatch, mZMatch] = Func_Compute3DXp(mCamera1Rectified, mCamera2Rectified, mXcMatch, 0, 0);
vProjector = mProjector * [mXMatch(:), mYMatch(:), mZMatch(:), ones(numel(mZMatch), 1)]';
vProjector = vProjector ./ vProjector(3, :);
mProjectorX = reshape(vProjector(1, :), size(mPhaseWrapLeft));
mPhiReference = mProjectorX ./ iProjectorWidth .* (2*pi);

% Unwrap the phase and triangulate the final 3D point cloud.
mPhaseOrder = round((FSet(end)*mPhiReference - mPhaseWrapLeft) ./ (2*pi));
mPhase = (mPhaseOrder.*(2*pi) + mPhaseWrapLeft) ./ FSet(end);
mXpUnwrap = mPhase ./ (2*pi) .* iProjectorWidth;

[mXRaw, mYRaw, mZRaw] = Func_Compute3DXp(mCamera1Rectified, mProjector, mXpUnwrap, 0, 0);
[mX, mY, mZ] = Func_Compute3DXp(mCamera1Rectified, mProjector, mXpUnwrap, 0, 0);
end
