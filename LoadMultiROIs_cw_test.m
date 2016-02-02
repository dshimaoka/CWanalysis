function [ROI, roisuffix] = LoadMultiROIs_cw_test(loadDir, expid)
    
% [ROI, roisuffix] = LoadMultiROIs_test(loadDir, expid)
%    
% TO DO: 
% when magnification factor is different across exps (expid=108)
% input resize factor

DIRS.Stacks = '\\zserver3.ioo.ucl.ac.uk\Data\Stacks';
shapeSuffix = '2Dgauss_minor';

Exps = readExpsDatabase('ExpsDatabase_cw.m', expid);
p = ProtocolLoadDS(Exps);
   
switch expid
        case {88, 91} %animal 66
        roisuffix = {'_stationary','_stationary_lm','_stationary_al','_stationary_rl',...
            '_stationary_am','_barrel','_auditory','_stationary__active_mslimb'};
        expid_roi = [91 91 91 91 91 92 92 88];
        iStim_roi = [ 1  1  1  1  1  1  2  1];
        useAllStim = [1  1  1  1  1  0  0  1];
        dataSuffix_roi = {'ar gd','ar gd','ar gd','ar gd','ar gd','','','ar gd'};
        
        expGroup = {[91 92], 88};%this can be automatically determined from p-files
        
    case {40, 41, 44} %animal 56
        roisuffix = {'_stationary','_stationary_lm','_stationary_al','_stationary_rl',...
            '_stationary_am','_barrel','_auditory'};%mslimb is not clear
        expid_roi = [40 40 40 40 40 42 42];
        iStim_roi =  [1 1 1 1 1 1 2];
        useAllStim = [1 1 1 1 1 0 0];
        dataSuffix_roi = {'ar gd','ar gd','ar gd','ar gd','ar gd','',''};
        
        expGroup = {[40 41 42], 44};
        
    case {53, 65} %animal 58
        roisuffix = {'_stationary','_stationary_lm','_stationary_rl','_barrel',...
            '_stationary__active_mslimb'};
        expid_roi = [53 53 53 66 65];
        iStim_roi = [1 1 1 1 1];
        useAllStim = [1 1 1 0 1];
        dataSuffix_roi = {'ar gd','ar gd','ar gd','','ar gd'};
        
        expGroup = {[53 66], 65};
        
    case {27, 38} %animal 52
        roisuffix = {'_stationary','_stationary_lm','_stationary_rl','_stationary_am',...
            '_barrel','_auditory','_stationary__active_mslimb'};
        expid_roi = [38 38 38 38 39 39 27];
        iStim_roi = [1 1 1 1 1 2 1];
        useAllStim = [1 1 1 1 0 0 1];
        dataSuffix_roi = {'ar gd','ar gd','ar gd','ar gd','','', 'ar gd'};
        expGroup = {27,[38 39]};
        
    case {109, 108, 116, 117} %animal 79
        roisuffix = {'_stationary','_stationary_lm','_stationary_rl',...
            '_barrel','_auditory','_stationary__active_mslimb','_al','_pm'};
        expid_roi = [109 109 109 107 107 108 105 105];
        iStim_roi = [1 1 1 1 2 1 1 1];
        useAllStim = [1 1 1 0 0 1 1 1];
        dataSuffix_roi = {'ar gd','ar gd','ar gd','','','ar gd','',''};
        expGroup = {[108 109], [107 105], [116 117]};
        
    case {19, 32, 34} %animal 53
        roisuffix = {'_stationary','_stationary_lm','_stationary_al','_stationary_rl',...
            '_barrel','_auditory','_stationary__active_mslimb'};
        expid_roi = [32 32 32 32 33 33 34];
        iStim_roi =  [1 1 1 1 1 2 1];
        useAllStim = [1 1 1 1 0 0 1];
        dataSuffix_roi = {'ar gd','ar gd','ar gd','ar gd','','','ar gd'};
        expGroup = {[32 33 34], 19};
        
    case {70} % animal 62
        roisuffix = {'_lm','_al','_rl',...
            '_am','_barrel','_stationary__active_mslimb'};
        expid_roi = [70 70 70 70 69 70];
        iStim_roi =  [3  3  3  3  1 1];
        useAllStim = [0  0  0  0  0 1];
        dataSuffix_roi = {'detrend','detrend','detrend','detrend','detrend','ar gd'};
        expGroup = {[69 70]};
        
    case {83, 90} %animal 65
        roisuffix = {'_al','_lm','_rl','_am','_pm','_barrel','_auditory',...
            '_stationary__active_mslimb'};
        expid_roi = [89 89 89 89 89 113 113 83];
        iStim_roi =  [1  1  1  1  1  1   2   1];
        useAllStim = [1  1  1  1  1  0  0  1];
        dataSuffix_roi = {'','','','','','','','ar gd'};
        expGroup = {83, [89 90], 113};
        
    case 100 %animal 77
        roisuffix = {'_al','_lm','_rl','_am','_pm'};
        expid_roi = [101 101 101 101 101];
        iStim_roi = [1 1 1 1 1];
        useAllStim = [1 1 1 1 1];
        dataSuffix_roi = {'','','','','',''};
        expGroup = {[99 100],101};
        
% for cw 
%     case{24, 25, 26, 27} %animal 65
%         roisuffix = {'_v1','_rl','_lm','_am','_pm','_barrel','_auditory'};
%         expid_roi = [122 122 122 122 122 113 113];
%         iStim_roi =  [1  1   1   1   1   1   2 ];
%         useAllStim = [1  1   1   1   1   0   0];
%         dataSuffix_roi = {'','','','','','',''};
%         expGroup = {[24 25 26 27 122], 113};

end

nROI = length(roisuffix);

 %% load ROI
        ROI_c = [];
        for rr = 1:nROI
            Exps_roi = readExpsDatabase('ExpsDatabase_tg.m', expid_roi(rr));
            p_roi = ProtocolLoadDS(Exps_roi);
            MyStackDir_roi = tools.getDirectory( loadDir, p_roi, Exps_roi.ResizeFac, ...
                iStim_roi(rr), 1, Exps_roi.Cam, dataSuffix_roi{rr});
            
            if useAllStim(rr)
                ROIname = sprintf('%s/ROI%s_%s_%d_%d_allStim_%s', MyStackDir_roi, ...
                shapeSuffix, Exps_roi.animal, Exps_roi.iseries,...
                Exps_roi.iexp,  roisuffix{rr});
            else
                ROIname = sprintf('%s/ROI%s_%s_%d_%d_%d_%s', MyStackDir_roi, ...
                shapeSuffix, Exps_roi.animal, Exps_roi.iseries,...
                Exps_roi.iexp, iStim_roi(rr), roisuffix{rr});
            end
            
            %load cannot be used inside parfor...
            load(ROIname, 'ROI');
            
            %judge if exp for the trace and for ROI belong same exp group
            for ll = 1:length(expGroup)
                thisGroup = sum(ismember(expGroup{ll}, expid));
                if thisGroup
                    sameGroup = sum(ismember(expGroup{ll}, expid_roi(rr)));
                end
            end
                
            if ~sameGroup
                p2exps{1} = p_roi; %to be rotated
                
                Exps_hack = readExpsDatabase('ExpsDatabase_tg.m', 122);%24/11/15
                p_hack = ProtocolLoadDS(Exps_hack);%24/11/15

                p2exps{2} = p_hack;
                [t_concord, fliplr_input, fliplr_base, flipud_input, flipud_base, ...
                    inputImage, inputImage_registered, baseImage] = ...
                    tools.LoadRotateInfo_2exps( DIRS.Stacks, p2exps, [Exps_roi.ResizeFac Exps_hack.ResizeFac ]);
                
                %                 R = makeresampler({'cubic', 'nearest'}, 'symmetric');
                %                R = makeresampler('nearest', 'symmetric');
                R = makeresampler('nearest', 'bound');

                ROI_rotated = [];
                for rrr = 1:size(ROI,3)
                    ROI_rotated(:,:,rrr) = imtransform(ROI(:,:,rrr), ...
                        t_concord, R, ...
                        'XData',[1 size(baseImage,2)], 'YData',[1 size(baseImage,1)]);
                end
                ROI = ROI_rotated;
            end
            
            if rr == 1
                ROI_c = squeeze(ROI(:,:,1));
            else
                ROI_c = cat(3, ROI_c, squeeze(ROI(:,:,1)));
            end
        end
        ROI = ROI_c;
        