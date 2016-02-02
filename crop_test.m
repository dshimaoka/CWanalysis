rowRegion = [350 1811];
colRegion = [1 2080];
binning = 4;

% rowRegion_binned = round(rowRegion/binning);
% colRegion_binned = round(colRegion/binning);

%load('\\zserver3\Data\Stacks\Affogato52_ratio\20150219\005\Affogato52_20150219_5rotate_flip_cams.mat')
load('\\ZSERVER3\Data\Stacks\M150108_SD2_ratio\001\002\M150108_SD2_1_2rotate_flip_cams.mat')
cam2_crop = cam2;

S = tools.LoadMyStacks('\\ZSERVER3\Data\Stacks\M150108_SD2_cam2\001\003\Resize 100 Stim 001 Repeat All', 'Stack0001');
cam2_max = squeeze(mean(S.Values(:,:,50:60),3));
cam2_max = cam2_max ./ mean(cam2_max(:));

%% convert cropped to max coordinate
xidxcrop = 1:size(cam2_crop,2);
yidxcrop = 1:size(cam2_crop,1);

% coordinates in max image ... length should be identical to cropped one
[xidxcrop_inmax, yidxcrop_inmax] = idxUncrop(xidxcrop, yidxcrop, rowRegion, colRegion,binning);

%these should be identical
imshowpair(cam2_crop, cam2_max(yidxcrop_inmax, xidxcrop_inmax));
hold on
contour(edge(cam2_crop,'canny'),1,'r');
contour(edge(cam2_max(yidxcrop_inmax, xidxcrop_inmax),'canny'),1,'g');

%% convert max to cropped coordinate
xidxmax = 1:size(cam2_max,2);
yidxmax = 1:size(cam2_max,1);

% coordinates in cropped image ... length should be identical to max one
[xidxmax_incrop, yidxmax_incrop] = idxcrop(xidxmax, yidxmax, rowRegion, colRegion, binning);
% xidxmax - rowRegion_binned(1);%suppose images are not flipped
% yidxmax_incrop = yidxmax - colRegion_binned(1);%suppose images are not flipped


filledImage = zeros(size(cam2_max,1), size(cam2_max,2));
filledImage(find(yidxmax_incrop==1):find(yidxmax_incrop==yidxcrop(end)), ...
    find(xidxmax_incrop==1):find(xidxmax_incrop==xidxcrop(end)))...
    = cam2_crop;

%these should be identical
%imshowpair(filledImage, cam2_max);


% imagesc(cam2_max);hold on
% rectangle('position',[xidxcrop_inmax(1) yidxcrop_inmax(1) xidxcrop_inmax(end)-xidxcrop_inmax(1) yidxcrop_inmax(end)-yidxcrop_inmax(1)])

% input_points(:,2) = input_points(:,2) + colRegion_binned(1);
% base_points(:,2) = base_points(:,2) + colRegion_binned(1);
% 
% input_points(:,1) = input_points(:,1) + rowRegion_binned(1);
% base_points(:,1) = base_points(:,1) + rowRegion_binned(1);




