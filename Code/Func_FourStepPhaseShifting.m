% This function computes wrapped phase using the standard four-step
% phase-shifting algorithm.
%
% Input parameters:
% vmI       -- four phase-shifted fringe images, format: H-by-W-by-4
% iBcThresh -- pixels with modulation below this threshold will be disabled
%
% Output:
% mPhi      -- wrapped phase map in the range [0, 2*pi)
% mBc       -- fringe modulation
%
% Author: Anonymous Authors, 2026/06/01
function [mPhi, mBc] = Func_FourStepPhaseShifting(vmI, iBcThresh)
mI0 = vmI(:, :, 1);
mI1 = vmI(:, :, 2);
mI2 = vmI(:, :, 3);
mI3 = vmI(:, :, 4);

mSin = mI1 - mI3;
mCos = mI0 - mI2;
mPhi = mod(atan2(mSin, mCos), 2*pi);
mBc = 0.5*hypot(mSin, mCos);
mPhi(mBc < iBcThresh) = nan;
end
