
cd 'C:\Users\daisuke\Documents\MATLAB\CWanalysis'
clear

localDir = '\\zombie\Users\daisuke\Documents\MATLAB\stackset\Data';
saveDir = '\\zombie\Users\daisuke\Documents\MATLAB\CWanalysis\CWanalysis20151125_z_all';
mkdir(saveDir);
doTraces = true;
scrsz = get(0, 'screensize');
whitebg(0,'w');close;
set(0,'defaultfigurecolor','w');
   



iStim = 1;
SetDefaultDirs;
DIRS.Stacks = '\\zserver3.ioo.ucl.ac.uk\Data\Stacks';
suffix = 'ar gd z vsM151015_SD_1_1_all';


tgtexps = [24 25 26];%[1:4 6 8];

evoke = 'stimCue';%'onsetTone';%'feedbackSt';%



%trial condition
resptime = '';
nrep = 'nrep';
contrastList = [0 .12 .25 .5];
% allContrast = sort(100*[-contrastList contrastList]);
allContrast = (100*[-contrastList contrastList]);
ntcontrasts = 2*length(contrastList);

%%output
Savg = [];
avgtcourse = [];
setcourse = [];
nEvents = nan(2*length(contrastList),3,max(tgtexps));
avgtcontResp = [];
setcontResp = [];

%% parameters for visualization
baseTime = [-0.1 0];
respTime = [0.1 0.2];


for cc = 1:2 %stimulus side
    for rr = 0:3 %response side
        for cidx = 1:length(contrastList)
            clear Sall
            
            %output
            tcourseall = [];
            numValidEvents = 0;
            
            
            stimContrast = contrastList(cidx);
            
            
            if cc==1
                cue='l';
            elseif cc==2
                cue='r';
            end
            
            %             if cc == 1
            %                 contrastIdx =  1 + length(contrastList) - cidx;
            %             elseif cc == 2
            %                 contrastIdx = length(contrastList) + cidx;
            %             end
            
            contrastIdx = length(contrastList)*(cc-1) + cidx; %25/11/15
            
            if rr==0
                resp = [];
            elseif rr==1
                resp='l';
            elseif rr==2
                resp='r';
            elseif rr==3
                resp='o';
            end
            
            disp(['Cue: ' cue ', Resp: ' resp ', Contrast: ' num2str(stimContrast)]);
            
            eee = 1;
            for cwExp = tgtexps;
                
                Exps = readExpsDatabase('ExpsDatabase_cw.m', cwExp);
                [p, block] = ProtocolLoadDS(Exps);
                
                
                MyStackDir = tools.getDirectory( localDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix);
                
                %% load MEAN across events
                %S_avgEvent_h = S_nanMean.GetOneCondition(2);%S_allEvent_h.AverageConditions;
                %S_avgEvent_v = S_nanMean.GetOneCondition(1);%S_allEvent_v.AverageConditions;
                
                %ThisFileName_h = sprintf('avgEvent_h_%s_%s_%s_%d', evoke, cue, resp, 100*stimContrast);
                ThisFileName_v = sprintf('avgEvent_v_%s_%s_%s_%d', evoke, cue, resp, 100*stimContrast);
                
                if exist([fullfile(MyStackDir,ThisFileName_v) '.mat'], 'file')
                    
                    %tools.SaveMyStacks(S_avgEvent_h, StackDir, DIRS.Temp, ThisFileName_h);
                    S = tools.LoadMyStacks(MyStackDir, ThisFileName_v);
                    
                    
                    %Spop{contrastIdx,rr+1} = Spop{contrastIdx,rr+1}.AssignCondition(S.Values);
                    if ~exist('Sall','var')
                        Sall = S; %Sall is a stackset containing all sessions of particular stim/resp condition
                        aaa = 1;
                    else
                        aaa = aaa+1;
                        Sall = Sall.AssignCondition(S.Values,aaa);
                    end
                end
                
                %% traces of all trials
                %                     [timeCourse, evTimeVec, validEventTimes] =...
                %                         AllEventTcourse(DIRS.Stacks, Cam, p, iStim, tgtTr, ...
                %                         ResizeFac, tevoke_selected, window, ROI, zsc, suffix);
                
                ThisFileName_timeCourse = sprintf('singleTcourse_%s_%s_%s_%d', evoke, cue, resp, 100*stimContrast);
                
                if exist([fullfile(MyStackDir,ThisFileName_timeCourse) '.mat'], 'file')
                    load(fullfile(MyStackDir,ThisFileName_timeCourse), 'timeCourse',...
                        'evTimeVec','validEventTimes','ROI','roisuffix','tevoke_selected');
                    
                    %                 TIMECOURSE{contrastIdx,rr+1} = cat(4, TIMECOURSE{contrastIdx,rr+1}, timeCourse);
                    tcourseall = cat(4, tcourseall, timeCourse);
                    %                     numValidEvents(eee) = size(validEventTimes, 1);
                    %
                    %                     nEvents(contrastIdx, rr+1, cwExp) = numValidEvents(eee);
                    %
                    %
                    %                     eee = eee+1;
                end
                
                clear timeCourse validEventTimes
                
                
                ThisFileName_validEventTimes = sprintf('validEventTimes_%s_%s_%s_%s.mat', ...
                    evoke, cue, resp, num2str(100*stimContrast));
                if exist(fullfile(MyStackDir,ThisFileName_validEventTimes), 'file')
                    load(fullfile(MyStackDir,ThisFileName_validEventTimes), 'validEventTimes');
                    numValidEvents(eee) = size(validEventTimes, 1);
                    
                    nEvents(contrastIdx, rr+1, cwExp) = numValidEvents(eee);
                    eee = eee+1;
                end
            end
            
            totEvents = sum(numValidEvents);
            weightEvents = numValidEvents/totEvents;
            
            if totEvents>0
                %Savg is a stackset containing avg of sessions of particular stim/resp condition
                Sall.Description = ['avg across ' num2str(totEvents) 'events'];
                Savg{contrastIdx,rr+1} = Sall.AverageConditions(weightEvents);
                
                
                %% SUBTRACT BASELINE PER TRIAL 26/11/15
                evTimeVec = Sall.TimeVec;
                [~, baseIdx(1)] = min(abs(evTimeVec - baseTime(1)));
                [~, baseIdx(2)] = min(abs(evTimeVec - baseTime(2)));
                [~, respIdx(1)] = min(abs(evTimeVec - respTime(1)));
                [~, respIdx(2)] = min(abs(evTimeVec - respTime(2)));
                
                Savg{contrastIdx,rr+1} = Savg{contrastIdx,rr+1}.SubtractBase(baseTime(1), baseTime(2));
                
                %avg and SE time course across trials in one condition
                baseTavg_cache = repmat(mean(tcourseall(baseIdx(1):baseIdx(2),:,:,:)),...
                    size(tcourseall,1),1,1,1);
                
                tcourseall = tcourseall - baseTavg_cache;
                avgtcourse(:,:,:,contrastIdx, rr+1) = mean(tcourseall,4);
                setcourse(:,:,:,contrastIdx, rr+1) = 1/sqrt(totEvents) * std(tcourseall, 1, 4);
                
                
                %% contrast-response curve, time averaged
                
                %                 contResp = squeeze(mean(tcourseall(respIdx(1):respIdx(2),:,:,:)) - ...
                %                     mean(tcourseall(baseIdx(1):baseIdx(2),:,:,:)));
                
                % THIS NASTY SQUEEZE!!!
                contResp = mean(tcourseall(respIdx(1):respIdx(2),:,:,:));
                contResp = reshape(contResp, size(contResp,1),size(contResp,3),size(contResp,4));
                
                avgtcontResp(:,:,contrastIdx, rr+1) = mean(contResp,3);
                setcontResp(:,:,contrastIdx, rr+1) = 1/sqrt(totEvents) * std(contResp, 1, 3);
                
                
                %% save all single traces
                ThisFileName_timeCourse_all = ...
                    sprintf('tcourse_all%s_%s_%s_%s_%d', suffix, evoke, cue, resp, 100*stimContrast);
                save(fullfile(saveDir,ThisFileName_timeCourse_all), 'tcourseall',...
                    'evTimeVec','ROI','roisuffix');
            end
            %nEvents(contrastIdx, rr+1) = totEvents;
            
            clear Sall
        end
    end
end

for cwExp = tgtexps;
    rateL(:,cwExp) = nEvents(:,2,cwExp)./nEvents(:,1,cwExp);
    rateR(:,cwExp) = nEvents(:,3,cwExp)./nEvents(:,1,cwExp);
end
sumnEvents = nansum(nEvents, 3);
avgrateL = sumnEvents(:,2)./sumnEvents(:,1);
avgrateR = sumnEvents(:,3)./sumnEvents(:,1);
%avgrateO = sumnEvents(:,4)./sumnEvents(:,1);

% plot(allContrast, avgrateR,'r.-')
% hold on
% plot(allContrast, avgrateL,'b.-')
% plot(allContrast, avgrateO,'g.-')
% marginplot
% legend('choose right','choose left','timeout');
% ylabel('probability');
% xlabel('stimulus contrast (%)')
% savePaperFigure(gcf, ['bhv' num2str(tgtexps)]);

save(fullfile(saveDir,['avgEvent_pop' num2str(tgtexps) suffix]),'Savg','nEvents');

save(fullfile(saveDir,['avgEvent_pop' num2str(tgtexps) suffix]), 'avgtcourse','setcourse',...
        'evTimeVec','avgtcontResp','setcontResp','-append');




%% d-prime
dprimetcourse = [];
for cc = 1:2 %stimulus side
    for cidx = 1:length(contrastList)
        
        noRecords = false;
        
        if cc==1
            cue='l';
        elseif cc==2
            cue='r';
        end
        
        contrastIdx = length(contrastList)*(cc-1) + cidx; %25/11/15
        
        
        stimContrast = contrastList(cidx);
        
        for rr = 1:2
            switch rr
                case 1
                    resp = 'r';
                case 2
                    resp = 'l';
            end
            
            ThisFileName_timeCourse_all = sprintf('tcourse_all%s_%s_%s_%s_%d', suffix, evoke, cue, resp, 100*stimContrast);
            
            if exist(fullfile(saveDir, [ThisFileName_timeCourse_all,'.mat']), 'file')
                load(fullfile(saveDir,ThisFileName_timeCourse_all), 'tcourseall',...
                    'evTimeVec','ROI');
            else
                noRecords = true;
                break;
            end
            
            switch rr
                case 1
                    tcourseall_r = tcourseall;
                case 2
                    tcourseall_l = tcourseall;
            end
        end
        
        if ~noRecords
            tcourseall_rl = cat(4, tcourseall_r, tcourseall_l);
            sdtcourseall_rl = std(tcourseall_rl, 1, 4);
            
            %is this correct?
            dprimetcourse(:,:,:,contrastIdx) = (mean(tcourseall_r,4) - mean(tcourseall_l,4)) ...
                ./ sdtcourseall_rl;
        end
        clear tcourseall_r tcourseall_l
    end
end


%% visualization

tgtConstIdx = [6 7 8]; %1:8;

%% figure 2: contrast-time image
for roiIdx = 1:length(roisuffix);
    figure('position',scrsz,'visible','on');
    
    for rr = 0:2 %response direction
        if rr==0
            resp = [];
        elseif rr==1
            resp='l';
        elseif rr==2
            resp='r';
        end
        
        for vh = 1%:2 %ratio
            subplot(3,2,2*rr+vh);
            imagesc(evTimeVec, allContrast(tgtConstIdx), squeeze(avgtcourse(:,vh,roiIdx,tgtConstIdx,rr+1))');
            set(gca,'ytick',allContrast(tgtConstIdx),'yticklabel',allContrast(tgtConstIdx));
            caxis([-1 1]);
            
            if rr==2
                xlabel(['Time from' evoke]);
            end
            if vh == 1
                title(sprintf('ratio\nresp %s, roi%d',resp, roiIdx))
            elseif vh == 2 && rr ==0
                title(['sum'])
            end
            
            mcolorbar;
            marginplot
        end
    end
    ThisFileName = sprintf('neurometric%s%s_%s_roi%d',suffix,num2str(tgtexps), evoke, roiIdx);
    screen2png(fullfile(saveDir,ThisFileName));
    close
end

%% figure 3: contrast response curve at particular time
figure('position',scrsz,'visible','on');
for roiIdx = 1:length(roisuffix);
    
    for vh = 1%:2 %ratio
        %subplot(4,2,2*rr+vh);
        %         subplot(2,length(roisuffix),roiIdx + (length(roisuffix))*(vh-1));
        subplot(length(roisuffix), 1, roiIdx);
        for rr = 1:2 %response direction
            if rr == 0;
                colorVec = [0 0 0];
            elseif rr == 1;
                colorVec = [0 0 1];
            elseif rr == 2;
                colorVec = [1 0 0];
            elseif rr == 3;
                colorVec = [0 1 0];
            end
            shadedErrorBar(allContrast(tgtConstIdx), squeeze(avgtcontResp(vh,roiIdx,tgtConstIdx, rr+1)),...
                squeeze(setcontResp(vh,roiIdx,tgtConstIdx, rr+1)),...
                {'-','markerfacecolor',colorVec,'markeredgecolor','none',...
                'color', colorVec}, 1);
            ylim([-0.2 1])
            
            set(gca,'xtick',allContrast(tgtConstIdx),'xticklabel',allContrast(tgtConstIdx));
            
            hold on
        end
        marginplot
    end
end
ThisFileName = sprintf('contResp%s%s_%s_roiall', suffix,num2str(tgtexps), evoke);
savePaperFigure(gcf,fullfile(saveDir,ThisFileName));
close


%% figure 4: mean+se time trace at one contrast
for tgtConstIdx = 1:8;
    figure('position',scrsz,'visible','on');
    for roiIdx = 1:length(roisuffix);
        for vh = 1 % rati
            subplot(1,length(roisuffix), roiIdx);
            for rr = [1 2]
                
                if rr == 1;%blue choose left
                    colorVec = [0 0 1];
                elseif rr == 2;%red choose right
                    colorVec = [1 0 0];
                elseif rr == 3;
                    colorVec = [0 1 0];
                end
                shadedErrorBar(evTimeVec, avgtcourse(:,vh,roiIdx,tgtConstIdx,rr+1), ...
                    setcourse(:,vh,roiIdx,tgtConstIdx,rr+1),...
                    {'-','markerfacecolor',colorVec,'markeredgecolor','none',...
                    'color', colorVec}, 1);
                hold on
                %xlim([evTimeVec(1) evTimeVec(end)]);
                xlim([-0.1 evTimeVec(end)]);
                ylim([-0.5 1.5]);
                
            end
            marginplot
        end
    end
    ThisFileName = sprintf('avgEvent_const%s%d_%s_%s_roiall', ...
        suffix,allContrast(tgtConstIdx), num2str(tgtexps), evoke);
    savePaperFigure(gcf,fullfile(saveDir,ThisFileName));
    close
end



%% figure 5: d' time trace at one contrast, with
for constIdx = 1:8;    
    figure('position',scrsz,'visible','on');
    for roiIdx = 1:length(roisuffix);
        for vh = 1%:2 % ratio and sum
            %             subplot(2,length(roisuffix),roiIdx + (length(roisuffix))*(vh-1));
            subplot(1,length(roisuffix), roiIdx);
            
            plot(evTimeVec, squeeze(dprimetcourse(:,vh,roiIdx,constIdx)));
            %xlim([evTimeVec(1) evTimeVec(end)])
            xlim([-0.1 0.6])
            ylim([-0.6 0.6])
            %             if vh == 1
            %                 title(sprintf('ratio\nresp roi%d',roiIdx))
            %             elseif vh == 2
            %                 title(['sum'])
            %             end
            marginplot
        end
    end
    ThisFileName = sprintf('dprime_const%s%d_%s_%s_roiall', suffix,allContrast(constIdx), num2str(tgtexps), evoke);
    savePaperFigure(gcf,fullfile(saveDir,ThisFileName));
    close
end
