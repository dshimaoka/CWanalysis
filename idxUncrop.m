function [xidx_uncrop, yidx_uncrop] = idxUncrop(xidx, yidx, rowRegion, colRegion, binning)
% [xidx_uncrop, yidx_uncrop] = idxUncrop(xidx, yidx, rowRegion, colRegion, binning)
% convert x and y idx of cropped image to idx in the uncropped image
%
% rowRegion: range in row of cropped image (before binning)
% colRegion: range in column of cropped image (before binning)
% binning: binning size (1,2,4,8), specified when recording

% assume images are not flipped/rotated
% assume resizefactor==1

%these should be identical
%imshowpair(cam2_crop, cam2_uncrop(yidx_uncrop, xidx_uncrop));

%maybe round is not accurate 
rowRegion_binned = round(rowRegion/binning);
colRegion_binned = round(colRegion/binning);

xidx_uncrop = xidx + rowRegion_binned(1);%suppose images are not flipped
yidx_uncrop = yidx + colRegion_binned(1);%suppose images are not flipped
