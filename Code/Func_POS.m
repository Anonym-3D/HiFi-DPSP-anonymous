% This function conducts phase optimization by solving a per-pixel sinusoidal
% least-squares model with a phase-offset sequence.
%
% Input parameters:
% vmI         -- fringe image sequence, format: H-by-W-by-NSet
% mDeltaPhase -- phase increment map between adjacent frames, format: H-by-W
% NSet        -- number of images used in the sinusoidal fitting
% iBcThresh   -- pixels with modulation below this threshold will be disabled
%
% Output:
% mPhase      -- optimized wrapped phase map
% mBc         -- fitted fringe modulation
%
% Author: Anonymous Authors, 2026/06/01
function [mPhase, mBc] = Func_POS(vmI, mDeltaPhase, NSet, iBcThresh)
[H, W] = size(mDeltaPhase);
mPhase = nan(H, W);
mBc = nan(H, W);

vStep = (0:NSet - 1)';
vValidPixel = find(isfinite(mDeltaPhase));
for iPixel = reshape(vValidPixel, 1, [])
    vDelta = vStep * (mDeltaPhase(iPixel) - pi/2);
    mA = [ones(NSet, 1), cos(vDelta), -sin(vDelta)];
    vI = reshape(vmI(iPixel + (0:NSet - 1)*H*W), [], 1);
    vSolve = mA \ vI;

    mPhase(iPixel) = pi + atan2(-vSolve(3), -vSolve(2));
    mBc(iPixel) = hypot(vSolve(3), vSolve(2));
end

mPhase(mBc < iBcThresh) = nan;
end
