% This script decodes the captured image sequences and visualizes the 3D
% reconstruction results of four-step PSP and HiFi-DPSP under three motion
% conditions.
%
% Output:
% A 2-by-3 figure showing shaded surface reconstruction results. The first
% row shows four-step PSP, and the second row shows HiFi-DPSP.
%
% Notes:
% The optical-flow and disparity maps loaded by this script are precomputed
% results from external algorithms. Optical flow is estimated with RAFT, and
% stereo disparity is estimated with RAFT-Stereo. Their implementations are
% not included in this MATLAB demo.
%
% References:
% Z. Teed and J. Deng, "RAFT: Recurrent All-Pairs Field Transforms for Optical Flow," ECCV, 2020.
% L. Lipson, Z. Teed, and J. Deng, "RAFT-Stereo: Multilevel Recurrent Field Transforms for Stereo Matching," 3DV, 2021.
%
% Author: Anonymous Authors, 2026/06/01
clc;
clearvars;
close all;

%% Parameter settings
sCodeDir = fileparts(mfilename('fullpath'));
sDataRoot = fullfile(sCodeDir, '..', 'Data');
sCalibRoot = fullfile(sCodeDir, '..', 'CalibMat');

vsFolder = {fullfile(sDataRoot, 'X-trans');
            fullfile(sDataRoot, 'Z-trans');
            fullfile(sDataRoot, 'Y-rot')};
vsMotionLabel = {'X Translation', 'Z Translation', 'Y Rotation'};

NSet = 8;
FSet = 912/20;
iCameraWidth = 800;
iCameraHeight = 600;
iBinomialOrder = 5;
iImageNum = iBinomialOrder + 3 + 2;
iBcThresh = 8;
fMaxTriangleEdge = 10;
vXlim = [-100, 100];
vYlim = [-150, 50];

%% Load calibration matrices
sCamera1 = load(fullfile(sCalibRoot, 'mCamera1Rectified.mat'), 'mCamera1Rectified');
sCamera2 = load(fullfile(sCalibRoot, 'mCamera2Rectified.mat'), 'mCamera2Rectified');
sProjector = load(fullfile(sCalibRoot, 'mProjector.mat'), 'mProjector');
mCamera1Rectified = sCamera1.mCamera1Rectified;
mCamera2Rectified = sCamera2.mCamera2Rectified;
mProjector = sProjector.mProjector;

%% Compare HiFi-DPSP and traditional four-step phase shifting for dynamic 3D scanning
figResult = figure;
set(figResult, 'Position', [0 0 1400 850]);

for iExp = 1:numel(vsFolder)
    vmIL_Raw = zeros(iCameraHeight, iCameraWidth, iImageNum + 2);
    vmIL_Raw_Next = zeros(iCameraHeight, iCameraWidth, iImageNum + 2);

    for iIdx = 1:iImageNum + 2
        vmIL_Raw(:, :, iIdx) = double(imread(fullfile(vsFolder{iExp}, '1_Rectified', sprintf('%04d.bmp', iIdx - 1))));
    end
    for iIdx = 1:iImageNum + 2
        iImageIndex = iImageNum + iIdx - 1;
        vmIL_Raw_Next(:, :, iIdx) = double(imread(fullfile(vsFolder{iExp}, '1_Rectified', sprintf('%04d.bmp', iImageIndex))));
    end
    
    % External RAFT-Stereo and RAFT results. The disparity maps are computed
    % from the left/right camera speckle image pairs, and the optical-flow
    % maps are computed between the first two uniformly illuminated frames.
    mDisparity1 = load(fullfile(vsFolder{iExp}, 'mDisparity_1.mat'), 'mDisparity_1');
    mDisparity2 = load(fullfile(vsFolder{iExp}, 'mDisparity_2.mat'), 'mDisparity_2');
    mFlowX1 = load(fullfile(vsFolder{iExp}, 'mFlowX_1.mat'), 'mFlowX_1');
    mFlowY1 = load(fullfile(vsFolder{iExp}, 'mFlowY_1.mat'), 'mFlowY_1');
    mFlowX2 = load(fullfile(vsFolder{iExp}, 'mFlowX_2.mat'), 'mFlowX_2');
    mFlowY2 = load(fullfile(vsFolder{iExp}, 'mFlowY_2.mat'), 'mFlowY_2');

    mDisparity_1 = mDisparity1.mDisparity_1;
    mDisparity_2 = mDisparity2.mDisparity_2;
    mFlowX_1 = mFlowX1.mFlowX_1;
    mFlowY_1 = mFlowY1.mFlowY_1;
    mFlowX_2 = mFlowX2.mFlowX_2;
    mFlowY_2 = mFlowY2.mFlowY_2;

    % Four-step PSP.
    [mPhase_FourStep, ~] = Func_FourStepPhaseShifting(vmIL_Raw(:, :, 3:6), iBcThresh);
    [mX_FourStep, mY_FourStep, mZ_FourStep, ~, ~, ~, ~] = Func_UnwrapAndCompute3D(mPhase_FourStep, mDisparity_1, mCamera1Rectified, mCamera2Rectified, mProjector, FSet);

    % HiFi-DPSP.
    [~, ~, ~, mPhase_BSC_Unwrap_1, vmIL] = Func_OF3IBSC(vmIL_Raw, mDisparity_1, mFlowX_1, mFlowY_1, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, iBcThresh);
    [~, ~, ~, mPhase_BSC_Unwrap_2] = Func_OF3IBSC(vmIL_Raw_Next, mDisparity_2, mFlowX_2, mFlowY_2, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, iBcThresh);
    mPhase_BSC_Unwrap_2_Aligned = Func_ImageAlign(mPhase_BSC_Unwrap_2, mFlowX_1, mFlowY_1, iImageNum, iImageNum);
    mPhaseOffset = FSet * (mPhase_BSC_Unwrap_2_Aligned - mPhase_BSC_Unwrap_1);
    [mPhase_3StepIBSC, ~] = Func_POS(vmIL(:, :, 1:NSet), mPhaseOffset./iImageNum, NSet, iBcThresh);
    [mX_HiFiDPSP, mY_HiFiDPSP, mZ_HiFiDPSP, ~, ~, ~, ~] = Func_UnwrapAndCompute3D(mPhase_3StepIBSC, mDisparity_1, mCamera1Rectified, mCamera2Rectified, mProjector, FSet);

    % Draw shaded triangular surfaces.
    fMeanDepth = mean(mZ_HiFiDPSP(:), 'omitnan');
    cMin = fMeanDepth - 60;
    cMax = fMeanDepth + 60;
    vClim = [cMin, cMax];
    vZlim = [cMin, cMax];

    figure(figResult);
    axFourStep = subplot(2, 3, iExp);
    Func_RenderTriSurface(axFourStep, mX_FourStep, mY_FourStep, mZ_FourStep, vClim, vXlim, vYlim, vZlim, fMaxTriangleEdge);
    title(vsMotionLabel{iExp}, 'FontSize', 20, 'FontWeight', 'bold');
    if iExp == 1
        ylabel(axFourStep, '4-Step PSP', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'black');
    end

    axHiFiDPSP = subplot(2, 3, iExp + numel(vsFolder));
    Func_RenderTriSurface(axHiFiDPSP, mX_HiFiDPSP, mY_HiFiDPSP, mZ_HiFiDPSP, vClim, vXlim, vYlim, vZlim, fMaxTriangleEdge);
    if iExp == 1
        ylabel(axHiFiDPSP, 'HiFi-DPSP', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'black');
    end
    drawnow;
end

hColorbarAxis = axes('Visible', 'off');
colormap(hColorbarAxis, jet);
cbDepth = colorbar(hColorbarAxis, 'Position', [0.92 0.18 0.016 0.68], 'FontSize', 20, 'FontWeight', 'bold');
set(get(cbDepth, 'Title'), 'String', 'mm', 'FontWeight', 'bold', 'FontSize', 20, 'Color', 'black');
set(cbDepth, 'FontWeight', 'bold', 'FontSize', 20, 'Color', 'black');
caxis(hColorbarAxis, [cMin, cMax]);
