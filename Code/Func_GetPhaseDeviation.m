% This function estimates the phase deviation caused by optical flow. It
% interpolates the reference projector phase at the flow-shifted position and
% subtracts the original reference phase.
%
% Input parameters:
% mdX   -- horizontal optical-flow field, format: H-by-W
% mdY   -- vertical optical-flow field, format: H-by-W
% mMask -- valid-pixel mask; invalid pixels are set to NaN in the output
% mC    -- 3-by-4 projection matrix of the camera
% mP    -- 3-by-4 projection matrix of the projector
% H     -- image height
% W     -- image width
% FSet  -- fringe frequency used to scale the reference phase
%
% Output:
% mDeltaPhase -- optical-flow-induced phase deviation
%
% Author: Anonymous Authors, 2026/06/01
function mDeltaPhase = Func_GetPhaseDeviation(mdX, mdY, mMask, mC, mP, H, W, FSet)
mPhiRef = FSet .* Func_GenerateRefPhase(mC, mP, H, W);
[mx, my] = meshgrid(1:W, 1:H);

mPhiRefInterp = griddedInterpolant(my, mx, mPhiRef, 'linear', 'none');
mDeltaPhase = mPhiRefInterp(my + mdY, mx + mdX) - mPhiRef;
mDeltaPhase(~mMask) = nan;
end

function mPhiRef = Func_GenerateRefPhase(mC, mP, iCameraHeight, iCameraWidth)
% Generate the projector reference phase on the camera image plane.
iProjectorWidth = 1280;
mCameraPlane = [mC(:, 1:2), mC(:, 4)];
[mxc, myc] = meshgrid(1:iCameraWidth, 1:iCameraHeight);

vWorld = mCameraPlane \ [mxc(:), myc(:), ones(numel(mxc), 1)]';
vWorld = vWorld ./ vWorld(3, :);

vProjector = mP * [vWorld(1, :); vWorld(2, :); zeros(1, numel(mxc)); ones(1, numel(mxc))];
vProjector = vProjector ./ vProjector(3, :);
mProjectorX = reshape(vProjector(1, :), size(mxc));
mPhiRef = mProjectorX ./ iProjectorWidth .* (2*pi);
end
