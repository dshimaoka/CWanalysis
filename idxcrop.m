function [xidx_crop, yidx_crop] = idxcrop(xidx, yidx, rowRegion, colRegion, binning)
% [xidx_crop, yidx_crop] = idxcrop(xidx, yidx, rowRegion, colRegion, binning)
% convert x and y idx of uncropped image to idx in the cropped image, so
% outcome can be negative value
%
% rowRegion: range in row of cropped image (before binning)
% colRegion: range in column of cropped image (before binning)
% binning: binning size (1,2,4,8), specified when recording

% assume images are not flipped/rotated
% assume resizefactor==1

% filledImage = zeros(size(cam2_uncrop,1), size(cam2_uncrop,2));
% filledImage(find(yidx_crop==1):find(yidx_crop==size(cam2_crop,1)), ...
%     find(xidx_crop==1):find(xidx_crop==size(cam2_crop,2)))...
%     = cam2_crop;

%these should be identical
%imshowpair(filledImage, cam2_crop);

%maybe round is not accurate 
rowRegion_binned = round(rowRegion/binning);
colRegion_binned = round(colRegion/binning);

xidx_crop = xidx - rowRegion_binned(1);%suppose images are not flipped
yidx_crop = yidx - colRegion_binned(1);%suppose images are not flipped
