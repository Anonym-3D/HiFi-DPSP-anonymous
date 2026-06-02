% This function conducts three-step image-sequential binomial self-compensation
% (I-BSC) for motion-error suppression. It first accumulates the binomially
% compensated sine and cosine terms, then computes the motion-error-free phase.
%
% Input parameters:
% vmI       -- captured image sequence, format: H-by-W-by-(K+3), where H-by-W is the image resolution and K is the binomial order
% iBcThresh -- pixels with modulation below this threshold will be disabled
%
% Output:
% mPhi      -- motion-error-free wrapped phase obtained through I-BSC
% mBc       -- modulation of the compensated fringe sequence
%
% Author: Anonymous Authors, 2026/06/01
function [mPhi, mBc] = Func_3StepIBSC(vmI, iBcThresh)
[~, ~, iImageNum] = size(vmI);
iBinomialOrder = iImageNum - 3;

% Compute compensated sine and cosine terms.
mSin = 0;
mCos = 0;
for k = 0:iBinomialOrder
    fBinomial = nchoosek(iBinomialOrder, k);
    fCosShift = cos(k*pi/2);
    fSinShift = sin(k*pi/2);

    mSecondDifference = -vmI(:, :, k + 1) + 2*vmI(:, :, k + 2) - vmI(:, :, k + 3);
    mFringeDifference =  vmI(:, :, k + 1) - vmI(:, :, k + 3);

    mSin = mSin + fBinomial*( fCosShift*mSecondDifference + fSinShift*mFringeDifference);
    mCos = mCos + fBinomial*(-fSinShift*mSecondDifference + fCosShift*mFringeDifference);
end

% Compute phase and modulation.
mPhi = pi + atan2(-mSin, -mCos);
mBc = sqrt(mSin.^2 + mCos.^2) ./ 2^(1 + iBinomialOrder);

% Reject low-modulation pixels.
mPhi(mBc < iBcThresh) = nan;
end
