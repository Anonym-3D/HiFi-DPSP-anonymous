% This function aligns an image according to the fractional motion between
% two frames. It keeps static pixels unchanged, interpolates pixels that remain
% dynamic, and suppresses static-to-dynamic transition regions.
%
% Input parameters:
% mImg        -- image or phase map to be aligned, format: H-by-W
% mdX         -- horizontal optical-flow field over the full frame interval
% mdY         -- vertical optical-flow field over the full frame interval
% iFrameNum   -- target frame index within the interval
% iFrameTotal -- total number of frame steps in the interval
%
% Output:
% mImg_Align  -- aligned image or phase map
%
% Author: Anonymous Authors, 2026/06/01
function mImg_Align = Func_ImageAlign(mImg, mdX, mdY, iFrameNum, iFrameTotal)
mOffsetX = mdX ./ iFrameTotal .* iFrameNum;
mOffsetY = mdY ./ iFrameTotal .* iFrameNum;
[H, W] = size(mImg);
[mx, my] = meshgrid(1:W, 1:H);
mImg_Align = zeros(H, W, 'like', mImg);

% Estimate the current dynamic mask by forward splatting the initial mask.
mMotion_Init = hypot(mdX, mdY) > 0.1;
mxFull = mx + mdX;
myFull = my + mdY;
x1 = floor(mxFull);
x2 = ceil(mxFull);
y1 = floor(myFull);
y2 = ceil(myFull);
mValidFlow = mMotion_Init & x1 >= 1 & x2 <= W & y1 >= 1 & y2 <= H;

mMotion_Now = zeros(H, W);
if any(mValidFlow(:))
    dx = mxFull(mValidFlow) - x1(mValidFlow);
    dy = myFull(mValidFlow) - y1(mValidFlow);
    vY1 = y1(mValidFlow);
    vY2 = y2(mValidFlow);
    vX1 = x1(mValidFlow);
    vX2 = x2(mValidFlow);

    vIndex = [sub2ind([H, W], vY1, vX1);
              sub2ind([H, W], vY1, vX2);
              sub2ind([H, W], vY2, vX1);
              sub2ind([H, W], vY2, vX2)];
    vWeight = [(1 - dx).*(1 - dy);
                dx     .*(1 - dy);
               (1 - dx).* dy;
                dx     .* dy];
    mMotion_Now = reshape(accumarray(vIndex, vWeight, [H*W, 1], @sum, 0), H, W);
end
mMotion_Now = mMotion_Now > 0.4;

% Classify motion transitions.
mValidFull = x1 >= 1 & x2 <= W & y1 >= 1 & y2 <= H;
mDynamicToDynamic = false(H, W);
mStaticToDynamic = false(H, W);
if any(mValidFull(:))
    vRoundX = round(mxFull(mValidFull));
    vRoundY = round(myFull(mValidFull));
    vCurrentMotion = mMotion_Now(sub2ind([H, W], vRoundY, vRoundX));
    vInitialMotion = mMotion_Init(mValidFull);
    vValidFullIndex = find(mValidFull);

    mDynamicToDynamic(vValidFullIndex(vInitialMotion &  vCurrentMotion)) = true;
    mStaticToDynamic(vValidFullIndex(~vInitialMotion & vCurrentMotion)) = true;
end
mStaticToDynamic = imopen(mStaticToDynamic, strel('disk', 3));

% Align the image at the requested fractional frame position.
mxAlign = mx + mOffsetX;
myAlign = my + mOffsetY;
mValidAlign = floor(mxAlign) >= 1 & ceil(mxAlign) <= W & ...
              floor(myAlign) >= 1 & ceil(myAlign) <= H;
mInterpolated = interp2(mx, my, mImg, mxAlign, myAlign, 'linear', 0);

mCopyOriginal = mValidAlign & ~mDynamicToDynamic & ~mStaticToDynamic;
mImg_Align(mCopyOriginal) = mImg(mCopyOriginal);
mImg_Align(mValidAlign & mDynamicToDynamic) = mInterpolated(mValidAlign & mDynamicToDynamic);
end
