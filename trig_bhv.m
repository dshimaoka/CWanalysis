
rmpath('\\zombie\Users\daisuke\Documents\MATLAB\MATLAB_20140318_ZOOROPA')
scrsz = get(0, 'screensize');


iStim = 1;
evoke = 'moveOn';%'moveOn';%'onsetTone';% %'feedbackSt';%'onsetTone';%
if strcmp(evoke, 'moveOn')
    window = [-1 0.15];
elseif strcmp(evoke, 'stimCue')
    window = [-0.5 0.7];% analyzing period [s]
end

%trial condition
resptime = '';
nrep = 'nrep';

%create time traces of spedified ROI
doTraces = false;

%detection of wheel movement
useWheel = false; %whether to discard trials with wheel movements
ThSpeed = 1;%0.01; %[m/s] %0.2 is enough high to be detected in ratio movie
duration = 0.1;
winSpeed = 0.01;%0.1;%[s]

marginTime = 0.1; %[s]%for evokeTime.m

SetDefaultDirs;
DIRS.Stacks = '\\zserver3.ioo.ucl.ac.uk\Data\Stacks';
localDir = '\\zombie\Users\daisuke\Documents\MATLAB\stackset\Data';
shareDir = '\\zserver\Lab\Share\Shimaoka';
suffix_in = 'ar gd';
suffix_out = suffix_in;%'ar gd z vsM151015_SD_1_1';

if doTraces    
    %for expid=29
    ROI(1,:)=[48.2222222222223 92.4444444444445 9.88888888888887 11.5555555555555];
    ROI(2,:)=[86.4166666666667 46.5555555555556 4.88888888888889 14.8888888888889];
    roisuffix = {'_v1','_beta'};
end

    for cwExp = [24];
        
        %     [ROI,roisuffix] = LoadMultiROIs_cw_test(localDir,cwExp);
        %     ROI = imresize(ROI, [120 92],'nearest');%complete hack. should be in LoadMultiROIs
        
        
        Exps = readExpsDatabase('ExpsDatabase_cw.m', cwExp);

        for cams = 1:4
        switch cams
            case 1
                Exps.Cam.FileString = '_cam1';
            case 2
                Exps.Cam.FileString = '_cam2';
            case 3
                Exps.Cam.FileString = '_ratio';
            case 4
                Exps.Cam.FileString = '_sum';
        end

        acqEndDelay = getExtraDelayInfo(Exps);
        
        [p, block] = ProtocolLoadDS(Exps);
        MyStackDir = tools.getDirectory( localDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_out);
        mkdir(MyStackDir);
        
        zStackDir = tools.getDirectory( shareDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_out);
        mkdir(zStackDir);

        
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
        
        
        %11/3/15
        repeatList = availableStackRepeats(DIRS.Stacks, p, Exps.ResizeFac, iStim, Exps.Cam.FileString, suffix_in);
        repeatList = cell2mat(repeatList);
        
        clear S_allEvent_*;
        for cc = 0%1:2 %stimulus side
            for rr = 0%1:3 %response side
                for cidx = 1%:length(contrastList)
                    tevoke_total = [];
                    
                    stimContrast = [];%contrastList(cidx);
                    
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
                    
                    %% prepare event time (tevoke), to which data is aligned
                    if isempty(tgtTr)
                        tevoke_selected = [];
                        disp('No trials in this condition.')
                    else
                        disp([num2str(length(tgtTr)) ' trials detected.']);
                        
                        %acqStartTime = getacqStartTimeCW(Exps, tgtTr, Exps.Cam, suffix);
                        
                        if useWheel
                            Exps_evoke = Exps;
                            Exps_evoke.Cam.FileString = '_cam2';
                            [~,~, SPEED, TIME] = ...
                                tl.bodyMoveStatsTL_cw(Exps_evoke, tgtTr, 1, winSpeed);
                        end
                        
                        %STIMON = CWGetStimOnOff(Exps, tgtTr, '_cam2', 'ar');
                        Cam_evoke.FileString = '_cam2'; %hack
                        suffix_evoke = ''; %hack
                        STIMON = cell2mat(evokeTime(Exps, evoke, tgtTr, ...
                            marginTime, Cam_evoke, suffix_evoke, Timeline));
                        
                        tevoke_selected = cell(1);
                        tgtTr_selected = [];
                        kkk = 1;
                        for jjj = 1:length(tgtTr)
                            
                            iTr = tgtTr(jjj);
                            %disp(['trial: ' num2str(iTr)]);
                            
                            %% wheel movement
                            if useWheel
                                
                                motion = detectLocomotion(abs(SPEED{jjj}), TIME{jjj}, ...
                                    ThSpeed, duration);
                                
                                [~,winStIdx] = min(abs(TIME{jjj} - (STIMON(jjj)+window(1))));
                                %[~,winStIdx] = min(abs(TIME{jjj} - (STIMON(jjj)+0)));%test 27/11/15
                                [~,winEnIdx] = min(abs(TIME{jjj} - (STIMON(jjj)+window(2))));
                                
                                %whether motion is detected within the
                                %window around the specified event
                                motionInWindow = sum(motion(winStIdx:winEnIdx));
                                if ~motionInWindow
                                    tevoke_selected{kkk} = STIMON(jjj);
                                    tgtTr_selected(kkk) = iTr;
                                    kkk = kkk + 1;
                                end
                                
                            else
                                tevoke_selected{jjj} = STIMON(jjj);
                                tgtTr_selected(jjj) = tgtTr(jjj);
                            end
                        end
                        
                        disp([num2str(length(tgtTr_selected)) ' trials selected.']);
                        
                        %% average stackset across selected trials
                        if ~isempty(tevoke_selected{1})
                            
                            if ~doTraces
                                [S_avg, S_SE, validEventTimes, validFrames] = ...
                                    StackEvRepeats(DIRS.Stacks, Exps.Cam, p, iStim, tgtTr_selected, Exps.ResizeFac, ...
                                    tevoke_selected, window, [], [], suffix_in);
                            else
                                [S_avg, S_SE, validEventTimes, validFrames, timeCourse] =...
                                    StackEvRepeats(DIRS.Stacks, Exps.Cam, p, iStim, tgtTr_selected, Exps.ResizeFac, ...
                                    tevoke_selected, window, ROI, [], suffix_in);
                                
                                
                                Exps_evoke = Exps;
                                Exps_evoke.Cam.FileString = '_cam2';
                                [~,~, SPEED_c, TIME_c] = ...
                                    tl.bodyMoveStatsTL_cw(Exps_evoke, validEventTimes(:,1), 1, winSpeed);%slow
                                
                                for nt = 1:size(validEventTimes,1);
                                    thisTr = validEventTimes(nt,1);
                                    thisTime = validEventTimes(nt,2);
                                    
                                    ttrace = squeeze(timeCourse(:,:,:,nt));
                                    mttrace = median(ttrace,1);
                                    ttrace = ttrace - repmat(mttrace, size(ttrace,1),1);
                                    h1 = tools.showTtrace_Pspec(ttrace, S_avg.TimeVec, 'on');
                                    
                                    
                                    h2 = figure('visible','on');
                                    subplot(1,2,1);
                                    plot(TIME_c{nt} - thisTime, SPEED_c{nt});
                                    %xlim([TIME{thisTr}(winStIdx) TIME{thisTr}(winEnIdx)]);
                                    xlim([window(1) window(2)])
                                    title(['trial: ' num2str(thisTr)]);
                                    ylabel('wheel speed');
                                    h = mergefigs([h1; h2]);
                                    %clear h h1 h2
                                    close all
                                    screen2png(['trajectories_trial' num2str(thisTr)]);
                                end
                                
                            end
                            
                            evTimeVec = S_avg.TimeVec;
                            
                            ThisFileName_validEventTimes = sprintf('validEventTimes_%s_%s_%s_%s', ...
                                evoke, cue, resp, num2str(100*stimContrast));
                            if ~useWheel
                                ThisFileName_validEventTimes = [ThisFileName_validEventTimes '_all'];
                            end
                            
                            save(fullfile(MyStackDir,ThisFileName_validEventTimes), 'validEventTimes');
                            
                            %% save MEAN across events
                            S_avgEvent_v = S_avg.GetOneCondition(1);%S_allEvent_v.AverageConditions;
                            ThisFileName_v = sprintf('avgEvent_v_%s_%s_%s_%s', ...
                                evoke, cue, resp, num2str(100*stimContrast'));
                            if ~useWheel
                                ThisFileName_v = [ThisFileName_v '_all'];
                            end
                            tools.SaveMyStacks(S_avgEvent_v, MyStackDir, DIRS.Temp, ThisFileName_v);
                            
                            
                            %transient on 2/2/16
%                             save(fullfile(zStackDir,ThisFileName_validEventTimes), 'validEventTimes');
%                             tools.SaveMyStacks(S_avgEvent_v, zStackDir, DIRS.Temp, ThisFileName_v);
                            
                            %% save individual across events (warning: can be heavy!)
                            %                             [S_allEvent_v,validEventTimes, validFrames] =...
                            %                                 StackEvAll(DIRS.Stacks, Exps, 1, tgtTr_selected, 1, ...
                            %                                 tevoke_selected, window, [], suffix_in);
                            %                             ThisFileName_a = sprintf('allEvent_v_%s_%s_%s_%s', evoke, cue, resp, num2str(100*stimContrast'));
                            %                             if ~useWheel
                            %                                 ThisFileName_a = [ThisFileName_a '_all'];
                            %                             end
                            %                             tools.SaveMyStacks(S_allEvent_v, MyStackDir, DIRS.Temp, ThisFileName_a);
                            
                            %if doTraces
                            %% traces of all trials
                            %                             [timeCourse, evTimeVec, validEventTimes] =...
                            %                                 AllEventTcourse_test(DIRS.Stacks, Cam, p, iStim, tgtTr, ...
                            %                                 Exps.ResizeFac, tevoke_selected, window, ROI, [], suffix);
                            
                            
                            %                             ThisFileName_timeCourse = sprintf('singleTcourse_%s_%s_%s_%s', evoke, cue, resp, num2str(100*stimContrast'));
                            %                             save(fullfile(MyStackDir,ThisFileName_timeCourse), 'timeCourse',...
                            %                                 'evTimeVec','validEventTimes','ROI','roisuffix','tevoke_selected');
                            
                            clear timeCourse evTimeVec validEventTimes
                            %end
                            clear S_avg S_SE S_avgEvent_h S_avgEvent_v tevoke_total;
                            %end
                        end
                    end
                end
            end
        end
    end
end





