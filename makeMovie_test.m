ROI = [45.3821138211382 93.6750296409214 9.92805724932257 11.1796451558267;...
    30.5 80.25 9.75 11;...
6.33333333333334 78.1666666666667 8.91666666666667 12.25;...
18.9166666666667 56.3333333333334 9.33333333333334 12.6666666666667];

movstack = Savg{6,3};
movstack = movstack.Trim([-0.1 0.6]);

load('\\zserver3\Data\Stacks\M151015_SD_ratio\001\001\M151015_SD_1_1maskInfo.mat')
%mask = imresize(maskInfo.maskIdx2D,[120 92]);
mask_c = maskInfo.signal;
mask = imresize(mask_c./max(mask_c(:)), [120 92]);
movstack.PlayCondition(1,[-0.5 1.5],mask,'jet');
