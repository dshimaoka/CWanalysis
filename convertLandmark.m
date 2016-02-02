function [input_points_converted, base_points_converted] = ...
    convertLandmark(input_points, base_points, rowRegion, colRegion, binning,...
    fliplr_input, flipud_input, fliplr_base, flipud_base)
% [input_points_converted, base_points_converted] = ...
%    convertLandmark(input_points, base_points, rowRegion, colRegion, binning,...
%    fliplr_input, flipud_input, fliplr_base, flipud_base)
%converts input and base points in UNCROPPED coordinates, to CROPPED
%coordinates

%rowRegion and colRegion: [beginIdx_input endIdx_input; beginIdx_base endIdx_base]
%binning: binning factor, assuming CROPPED and NONCROPPED images have same
%binning factor

%2015-2-27 DS created

%TODO: also convert from CROPPED to UNCROPPED coordinates

if size(rowRegion,1) == 1
    rowRegion(2,:) = rowRegion(1,:);
    colRegion(2,:) = colRegion(1,:);
end


rowRegion_binned = round(rowRegion/binning);
colRegion_binned = round(colRegion/binning);


%column and row size of uncropped images
col_input = round(2560/binning); %corrsponding x in camware, and y in matlab image
row_input = round(2160/binning); %corresponding y in camware, and x in matlab image
col_base = col_input;
row_base = row_input;

%column and row size of cropped images 
col_input_crop = rowRegion_binned(1,2) - rowRegion_binned(1,1);
row_input_crop = colRegion_binned(1,2) - colRegion_binned(1,1);
col_base_crop = rowRegion_binned(2,2) - rowRegion_binned(2,1);
row_base_crop = colRegion_binned(2,2) - colRegion_binned(2,1);


%% convert input and base points before image flipping
base_points_original = base_points;
if strcmp(fliplr_base, 'y')%x-axis
    base_points_original(:,1) = row_base - base_points(:,1) + 1;
end
if strcmp(flipud_base, 'y')%y-axis
    base_points_original(:,2) = col_base - base_points(:,2) + 1;
end

input_points_original = input_points;
if strcmp(fliplr_input, 'y')
    input_points_original(:,1) = row_input - input_points(:,1) + 1;
end
if strcmp(flipud_input, 'y')
    input_points_original(:,2) = col_input - input_points(:,2) + 1;
end

%to check
% subplot(121);imagesc(inputImage);axis image;hold on
% plot(input_points_original(:,1), input_points_original(:,2), 'r.');
% 
% subplot(122);imagesc(baseImage);axis image;hold on
% plot(base_points_original(:,1), base_points_original(:,2), 'r.');


%% convert coordinates of input and base points to (un)cropped coordinates
[c_x, c_y] = idxcrop(input_points_original(:,1), ...
    input_points_original(:,2), rowRegion(1,:), colRegion(1,:), binning);
input_points_c_original(:,1) = c_x;
input_points_c_original(:,2) = c_y;

[c_x, c_y] = idxcrop(base_points_original(:,1), ...
    base_points_original(:,2), rowRegion(2,:), colRegion(2,:), binning);
base_points_c_original(:,1) = c_x;
base_points_c_original(:,2) = c_y;

%% convert input and base points in (un)cropped coordinates after image flipping
%need to know size of uncropped input and base beforehand

base_points_converted = base_points_c_original;
if strcmp(fliplr_base, 'y')
    base_points_converted(:,1) = col_input_crop  - base_points_c_original(:,1) + 1;
end
if strcmp(flipud_base, 'y')
    base_points_converted(:,2) = row_input_crop  - base_points_c_original(:,2) + 1;
end

input_points_converted = input_points_c_original;
if strcmp(fliplr_input, 'y')
    input_points_converted(:,1) = col_base_crop - input_points_c_original(:,1) + 1;
end
if strcmp(flipud_input, 'y')
    input_points_converted(:,2) = row_base_crop - input_points_c_original(:,2) + 1;
end

% to check
% subplot(121);imagesc(fliplr(flipud(cam2)));axis image;hold on
% plot(base_points(:,1), base_points(:,2), 'r.');
% 
% subplot(122);
% imagesc(fliplr(flipud(cam2_max)));axis image;hold on
% plot(base_points_unc(:,1), base_points_unc(:,2), 'r.');
% 
% 
% Exps = readExpsDatabase('ExpsDatabase_tg.m', 47);
% p{1} = ProtocolLoadDS(Exps);
% 
% tools.LoadRotateInfo('\\zserver3.ioo.ucl.ac.uk\Data\Stacks', p{1})


