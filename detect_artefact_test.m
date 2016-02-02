ServerDir = '\\zserver3.ioo.ucl.ac.uk\Data\Stacks';
SetDefaultDirs;
DIRS.Stacks = '\\zserver3.ioo.ucl.ac.uk\Data\Stacks';
localDir = '\\zombie\Users\daisuke\Documents\MATLAB\stackset\Data';
suffix_in ='ar gd';
suffix_out = suffix_in;%'ar gd z vsM151015_SD_1_1';

iStim = 1;
evoke = 'stimCue';%'moveOn';%'onsetTone';% %'feedbackSt';%'onsetTone';%
window = [-0.5 0.7];% analyzing period [s]

th = 0.1;% threshold to detect artefact using sum signal at bregma in percentage

marginTime = 0.1; %[s]%for evokeTime.m

resptime = '';
nrep = 'nrep';

cwExp = [27];


%     [ROI,roisuffix] = LoadMultiROIs_cw_test(localDir,cwExp);
%     ROI = imresize(ROI, [120 92],'nearest');%complete hack. should be in LoadMultiROIs


Exps = readExpsDatabase('ExpsDatabase_cw.m', cwExp);

acqEndDelay = getExtraDelayInfo(Exps);

[p, block] = ProtocolLoadDS(Exps);
MyStackDir = tools.getDirectory( localDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_out);
mkdir(MyStackDir);


contrast_cache = [];
for itr = 1:block.numCompletedTrials
    contrast_cache = [contrast_cache block.trial(itr).condition.visCueContrast];
end
%contrastList = setdiff(unique(contrast_cache), 0);
contrastList = unique(contrast_cache);


seriesStr = num2str(Exps.iseries);
serieName = [seriesStr(1:4) '-' seriesStr(5:6) '-' seriesStr(7:8)];

tlname = sprintf('//zserver.ioo.ucl.ac.uk/Data/expInfo/%s/%s/%d/%s_%d_%s_Timeline.mat',...
    Exps.animal, ...
    serieName, ...
    Exps.iexp, ...
    serieName, Exps.iexp, Exps.animal);

if exist(tlname, 'file')
    display('Loading Timeline structure..');
    load(tlname);
end


if cwExp == 29
    ROI(1,:)=[48.2222222222223 92.4444444444445 9.88888888888887 11.5555555555555];
    ROI(2,:)=[86.4166666666667 46.5555555555556 4.88888888888889 14.8888888888889];
    ROI(3,:) = [67.1564611607584 55.8888888888889 8.22222222222223 12.6666666666667];%putative limb area
elseif cwExp == 27
    %(aligned to M151015_SD)
    ROI(1,:) = [42.25 85.8333333333333 11.8333333333333 16.4166666666667];
    ROI(2,:) = [77.75 25.1666666666667 5.99999999999999 8.08333333333333];
    ROI(3,:) = [60.75 43.5 9.74999999999999 9.33333333333333];
    %     ROI(1,:) = [43.6666666666667 88.1666666666666 13.0833333333333 14.3333333333333];
    %     ROI(2,:) = [78.8333333333334 27.5833333333333 4.3333333333333 5.58333333333333];
    %     ROI(3,:) = [50.25 49.0833333333333 8.91666666666666 8.08333333333331];
   % th = 2;% SD
   
elseif cwExp == 24
    ROI(3,:) = [56.0833333333334 46.4166666666667 5.99999999999999 8.91666666666666];
    ROI(2,:) = [80.25 22.6666666666667 5.99999999999999 8.91666666666666];
    ROI(1,:) = [44 95.5833333333334 5.99999999999999 6.83333333333331];
end
roisuffix = {'_v1','_bregma','limb'};

%11/3/15
repeatList = availableStackRepeats(DIRS.Stacks, p, Exps.ResizeFac, iStim, Exps.Cam.FileString, suffix_in);
repeatList = cell2mat(repeatList);


for cc = 1:2%1:3 %stimulus side
    for rr = 0%1:3 %response side
        for cidx = 1:length(contrastList)
            tevoke_total = [];
            
            stimContrast = contrastList(cidx);
            
            StackDir = tools.getDirectory( DIRS.Stacks, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_in);
            if cc == 0
                cue ='';
            elseif cc==1
                cue='l';
            elseif cc==2
                cue='r';
            elseif cc==3
                cue='o';
            end
            
            if rr==0
                resp = [];
            elseif rr==1
                resp='l';
            elseif rr==2
                resp='r';
            elseif rr==3
                resp='o';
            end
            
            disp(['Cue: ' cue ', Resp: ' resp ', Contrast: ' num2str(stimContrast')]);
            
            tgtTr = selecttrials(block.trial, cue, resp, resptime,nrep,stimContrast);
            tgtTr = intersect(tgtTr, repeatList); %11/3/15
            
            if ~isempty(tgtTr)
                Cam_evoke.FileString = '_cam2'; %hack
                suffix_evoke = ''; %hack
                STIMON = cell2mat(evokeTime(Exps, evoke, tgtTr, ...
                    marginTime, Cam_evoke, suffix_evoke, Timeline));
                
                tevoke_selected = cell(1);
                tgtTr_selected = [];
                for jjj = 1:length(tgtTr)
                    
                    iTr = tgtTr(jjj);
                    
                    tevoke_selected{jjj} = STIMON(jjj);
                    tgtTr_selected(jjj) = tgtTr(jjj);
                end
                
                %tgtTr_selected;
                %tevoke_selected;
                
                
                % ngTrials = [28 45 48 110 127 198 279];
                % [~,~,ngTrIdx]=intersect(ngTrials, tgtTr_selected);
                % tevoke_ng = tevoke_selected(ngTrIdx);
                %
                % okTrials = setxor(tgtTr_selected, ngTrials);
                % [~,~,okTrIdx]=intersect(okTrials, tgtTr_selected);
                % tevoke_ok = tevoke_selected(okTrIdx);
                
                
                % stackDir = tools.getDirectory( DIRS.Stacks, p, Exps.ResizeFac, 1, 1, Exps.Cam, dataSuffix);
                
                for fs = 1:2
                    if fs==1
                        Exps.Cam.FileString = '_ratio';
                        %% trials
                        [S_avg, S_SE, validEventTimes, validFrames, traces_r] =...
                            StackEvRepeats(DIRS.Stacks, Exps.Cam, p, iStim, tgtTr_selected, Exps.ResizeFac, ...
                            tevoke_selected, window, ROI, [], suffix_in);
                        
                        
                        traces_r = traces_r - repmat(median(traces_r,1),size(traces_r,1),1,1,1);
                        
                        
                    elseif fs==2
                        Exps.Cam.FileString = '_sum';
                        %% trials
                        [S_avg, S_SE, validEventTimes, validFrames, traces_s] =...
                            StackEvRepeats(DIRS.Stacks, Exps.Cam, p, iStim, tgtTr_selected, Exps.ResizeFac, ...
                            tevoke_selected, window, ROI, [], suffix_in);
                        
                        ThisFileName_v = sprintf('avgEvent_v_%s_%s_%s_%s', ...
                                evoke, cue, resp, num2str(100*stimContrast'));
                        tools.SaveMyStacks(S_avg, MyStackDir,[],ThisFileName_v);
                        
                        traces_s = traces_s - repmat(median(traces_s,1),size(traces_s,1),1,1,1);
                        if ~exist('traces_s_cat','var')
                            traces_s_cat = traces_s;
                        else
                            traces_s_cat = cat(4,traces_s_cat, traces_s);
                        end
                    end
                end
                
                noiseIdx = find(max(abs(squeeze(traces_s(:,:,2,:))) > th));
                [~,nonoiseIdx] = setxor(1:size(validEventTimes,1),noiseIdx);
                %noiseTrials =  tgtTr_selected(noiseIdx);
                
                mtraces_r = squeeze(mean(traces_r,4));
                mtraces_r_noise = squeeze(mean(traces_r(:,:,:,noiseIdx),4));
                mtraces_r_nonoise = squeeze(mean(traces_r(:,:,:,nonoiseIdx),4));
                
                setraces_r = 1/sqrt(size(traces_r,1))*squeeze(std(traces_r,[],4));
                setraces_r_noise = 1/sqrt(length(noiseIdx))*squeeze(std(traces_r(:,:,:,noiseIdx),[],4));
                setraces_r_nonoise = 1/sqrt(length(nonoiseIdx))*squeeze(std(traces_r(:,:,:,nonoiseIdx),[],4));
                
                %% figure 1 avg trace across trials
                for ridx=1:3
                    subplot(3,1,ridx);
                    shadedErrorBar(S_avg.TimeVec, mtraces_r(:,ridx),setraces_r(:,ridx),'k',1)
                    hold on
                    shadedErrorBar(S_avg.TimeVec, mtraces_r_noise(:,ridx),setraces_r_noise(:,ridx),'r',1)
                    shadedErrorBar(S_avg.TimeVec, mtraces_r_nonoise(:,ridx),setraces_r_nonoise(:,ridx),'b',1)
                    
                    %                     plot(S_avg.TimeVec, mtraces_r(:,ridx),'k', ...
                    %                         S_avg.TimeVec, mtraces_r_noise(:,ridx),'r', ...
                    %                         S_avg.TimeVec, mtraces_r_nonoise(:,ridx),'b');
                    xlim([S_avg.TimeVec(1) S_avg.TimeVec(end)]);
                    ylabel([roisuffix{ridx} '   dR/R']);
                    grid on;
                    
                    if ridx == 1
                        title(['NG tr:' num2str(length(noiseIdx)) ',  OK tr:' num2str(length(nonoiseIdx)) ])
                        legend('all','noise','no noise','location','southwest');
                    end
                end
                screen2png(['avgTrace_contrast_' cue '_' num2str(100*stimContrast')]);
                clf
                
                %% figure 2 individual trace
                for ridx = 1:3
                    subplot(3,2,2*ridx-1);imagesc(S_avg.TimeVec, 1:length(nonoiseIdx), squeeze(traces_r(:,:,ridx,nonoiseIdx))');grid on;
                    set(gca,'ytick',1:length(nonoiseIdx),'yticklabel',tgtTr_selected(nonoiseIdx));
                    title(['OK trials ratio at ' roisuffix(ridx)]);mcolorbar(gca,0.5);
                    subplot(3,2,2*ridx);imagesc(S_avg.TimeVec, 1:length(nonoiseIdx), squeeze(traces_s(:,:,ridx,nonoiseIdx))');grid on;
                    set(gca,'ytick','');
                    title(['sum at ' roisuffix(ridx)]);mcolorbar(gca,0.5);                    
                end
                screen2png(['everyTrace_contrast_' cue '_' num2str(100*stimContrast')]);
                clf

            end
        end
    end
end

%% summary plot across contrasts
%% figure 3: sum signal at bregma
maxVal = max(abs(traces_s_cat));
test_c = maxVal(:,:,2,:);
hist(test_c(:),0:0.01:0.4)
xlabel('max(abs(dF/F))')
ylabel('#trials sum at lambda')
screen2png('sum_at_bregma');

% S_ng = S_ng.SubtractBase([],0);
% traces_ng = S_ng.SpaceAverages([],[],[],ROI);
% S_ng_avg = S_ng.AverageConditions;
%
%
% %% ok trials
% [S_ok] =...
%     StackEvAll(DIRS.Stacks, Exps, 1, okTrials, 1, ...
%     tevoke_ok, window, [], suffix_in);
%
% S_ok = S_ok.SubtractBase([],0);
% traces_ok = S_ok.SpaceAverages([],[],[],ROI);
% S_ok_avg = S_ok.AverageConditions;
%
% %combine two conditions
% S_all_avg = S_ok_avg;
% S_all_avg = S_all_avg.AssignCondition(S_ng_avg.Values,2);
%
% %% individual trace


% subplot(232);imagesc(S_ok.TimeVec, 1:S_ok.nConds, squeeze(traces_ok(:,:,2))');grid on;
% set(gca,'ytick',1:S_ok.nConds,'yticklabel',okTrials);
% title('OK trials at bregma');mcolorbar(gca,0.5);
% subplot(235);imagesc(S_ok.TimeVec, 1:S_ng.nConds, squeeze(traces_ng(:,:,2))');grid on;
% set(gca,'ytick',1:S_ng.nConds,'yticklabel',ngTrials);
% title('NG trials at bregma');mcolorbar(gca,0.5);
%
% subplot(233);imagesc(S_ok.TimeVec, 1:S_ok.nConds, squeeze(traces_ok(:,:,3))');grid on;
% set(gca,'ytick',1:S_ok.nConds,'yticklabel',okTrials);
% title('OK trials at limb');mcolorbar(gca,0.5);
% subplot(236);imagesc(S_ok.TimeVec, 1:S_ng.nConds, squeeze(traces_ng(:,:,3))');grid on;
% set(gca,'ytick',1:S_ng.nConds,'yticklabel',ngTrials);
% title('NG trials at limb');mcolorbar(gca,0.5);
%
% %% individual trace
% subplot(211);
% shadedErrorBar(S_all_avg.TimeVec',mean(squeeze(traces_ok(:,:,1)),2), std(squeeze(traces_ok(:,:,1)),[],2),'k',1)
% hold on
% shadedErrorBar(S_all_avg.TimeVec',mean(squeeze(traces_ng(:,:,1)),2), std(squeeze(traces_ng(:,:,1)),[],2),'r',1)
% title('at v1');
%
% subplot(212);
% shadedErrorBar(S_all_avg.TimeVec',mean(squeeze(traces_ok(:,:,2)),2), std(squeeze(traces_ok(:,:,2)),[],2),'k',1)
% hold on
% shadedErrorBar(S_all_avg.TimeVec',mean(squeeze(traces_ng(:,:,2)),2), std(squeeze(traces_ng(:,:,2)),[],2),'r',1)
% title('at bregma');
%
% %% deviation in time
% std_ok = squeeze(std(traces_ok,[],1));
% std_ng = squeeze(std(traces_ng,[],1));
% plot(std_ok(:,1), std_ok(:,2),'ko');
% hold on
% plot(std_ng(:,1), std_ng(:,2),'ro');
% legend('ok','ng');
% xlabel('v1 std');ylabel('bregma std');
%
%
%
%
