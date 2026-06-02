% This function reconstructs a 3D point cloud using optical-flow-assisted
% three-step image-sequential binomial self-compensation.
%
% Input parameters:
% vmIL_Raw          -- raw left-camera image sequence, format: H-by-W-by-(N+2)
% mDisparity        -- disparity map from the auxiliary camera to the main camera
% mFlowX            -- horizontal optical-flow field between uniformly illuminated images
% mFlowY            -- vertical optical-flow field between uniformly illuminated images
% mCamera1Rectified -- 3-by-4 projection matrix of the main camera
% mCamera2Rectified -- 3-by-4 projection matrix of the auxiliary camera
% mProjector        -- 3-by-4 projection matrix of the projector
% FSet              -- fringe frequency
% iBcThresh         -- pixels with modulation below this threshold will be disabled
%
% Output:
% mX,mY,mZ          -- reconstructed 3D point cloud
% mPhase_BSC_Unwrap_L -- unwrapped BSC phase map of the main camera
% vmIL              -- optical-flow-aligned fringe image sequence
%
% Author: Anonymous Authors, 2026/06/01
function [mX, mY, mZ, mPhase_BSC_Unwrap_L, vmIL] = Func_OF3IBSC(vmIL_Raw, mDisparity, mFlowX, mFlowY, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, iBcThresh)
[iCameraHeight, iCameraWidth, iRawImageNum] = size(vmIL_Raw);
iImageNum = iRawImageNum - 2;
mIL1 = 2*vmIL_Raw(:, :, 1);
mValidIllumination = mIL1 > 2*iBcThresh;

% Optical-flow-based image alignment.
mDeltaPhaseL_FPM = Func_GetPhaseDeviation(mFlowX, mFlowY, mValidIllumination, mCamera1Rectified, mProjector, iCameraHeight, iCameraWidth, FSet);
vmIL = zeros(iCameraHeight, iCameraWidth, iImageNum - 2, 'like', vmIL_Raw);
for m = 3:iImageNum
    vmIL(:, :, m - 2) = Func_ImageAlign(vmIL_Raw(:, :, m), mFlowX, mFlowY, m - 3, iImageNum) .* mValidIllumination;
end

% Three-step IBSC.
[mPhase_BSC, ~] = Func_3StepIBSC(vmIL, iBcThresh);
mTempL = 9 / (2*iImageNum) * mDeltaPhaseL_FPM;
mPhase_BSC = mod(mPhase_BSC - mTempL, 2*pi);

% Phase unwrapping and 3D reconstruction.
[mX, mY, mZ, ~, ~, ~, mPhase_BSC_Unwrap_L] = Func_UnwrapAndCompute3D(mPhase_BSC, mDisparity, mCamera1Rectified, mCamera2Rectified, mProjector, FSet);
mInvalid = (mZ > 250) | isnan(mZ);
mX(mInvalid) = nan;
mY(mInvalid) = nan;
mZ(mInvalid) = nan;
end
