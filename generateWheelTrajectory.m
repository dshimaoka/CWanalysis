function psychoDat = generateWheelTrajectory(allBlocks, figDestination, plotFit)
% generatePsychometric plots the psychometic data from the
% blocks entered and optionally fits the data with phychoFit or GLM
% (for 2AFC and 2AUC respectively), as well as saving the figure to the
% path 'figDestination'.
%
% allBlocks         Block structure or cell array of structs
% figDestination    String; Path where fig is to be saved or 'none'
% plotFit           Logical; If true, a plot of the fitted data is made
%
% psychoDat         Struct containing expRef, contrasts, response made,
%                   feedback, repeat numbers and response times

% 2015-12-20 DS created from generatePsychometric


%% Initialize & collect data
if isa(allBlocks, 'struct')
    allBlocks = {allBlocks};
elseif ~isa(allBlocks, 'struct')&&~isa(allBlocks, 'cell')
    error('allBlocks must be a single block structure or a cell array of blocks');
end

numCompletedTrials = 0;
for b = 1:length(allBlocks)
    numCompletedTrials = numCompletedTrials + allBlocks{b}.numCompletedTrials;
end

if numCompletedTrials == 0;
    psychoDat = struct(...
        'expRef',{''},...
        'contrast',[],...
        'resp',[],...
        'feedback',[],...
        'repeatNum',[],...
        'rt',[],...
        'rewardSize',[]);
    return
end

twindow = [-0.5 1.0];%[s]
trajectoryTimes = linspace(twindow(1), twindow(2), 200);
trajectory_SC = zeros(length(trajectoryTimes),numCompletedTrials);
trajectory_IS = zeros(length(trajectoryTimes),numCompletedTrials);

contrast = zeros(2,numCompletedTrials);
resp = zeros(1,numCompletedTrials);
feedback = zeros(1,numCompletedTrials);
repeatNum = zeros(1,numCompletedTrials);
rt = zeros(1,numCompletedTrials);
rewardSize = zeros(2,numCompletedTrials);
expRef = cell(1:length(allBlocks));
tInd = 1;
for b = 1:length(allBlocks)
    fOneInd = 1;
    expRef{b} = allBlocks{b}.expRef;
    %     resp = [allBlocks{b}.trial.responseMadeID];
    for t = 1:allBlocks{b}.numCompletedTrials
        contrast(:,tInd) = allBlocks{b}.trial(t).condition.visCueContrast;
        resp(tInd) = allBlocks{b}.trial(t).responseMadeID;
        feedback(tInd) = allBlocks{b}.trial(t).feedbackType;
        repeatNum(tInd) = allBlocks{b}.trial(t).condition.repeatNum;
        if allBlocks{b}.trial(t).feedbackType==1
            if size(allBlocks{b}.rewardDeliveredSizes,2)>1
                rewardSize(1,tInd) = allBlocks{b}.rewardDeliveredSizes(fOneInd,1);
                rewardSize(2,tInd) = allBlocks{b}.rewardDeliveredSizes(fOneInd,2);
            else
                rewardSize(tInd,1) = allBlocks{b}.rewardDeliveredSizes(fOneInd);
            end
            fOneInd = fOneInd+1;
        end
        rt(tInd) = allBlocks{b}.trial(t).responseMadeTime - allBlocks{b}.trial(t).interactiveStartedTime;
        
        
        
        %% calculate wheel trajectory, with respect to cue onset
        [~,startIdx_SC] = min(abs(allBlocks{b}.inputSensorPositionTimes - ...
            (allBlocks{b}.trial(t).stimulusCueStartedTime + twindow(1))));
        [~,endIdx_SC] = min(abs(allBlocks{b}.inputSensorPositionTimes - ...
            (allBlocks{b}.trial(t).stimulusCueStartedTime + twindow(2))));
        [~,zeroIdx_SC] = min(abs(allBlocks{b}.inputSensorPositionTimes - ...
            (allBlocks{b}.trial(t).stimulusCueStartedTime)));
        
        trajectoryTimes_c = allBlocks{b}.inputSensorPositionTimes(startIdx_SC:endIdx_SC)- allBlocks{b}.trial(t).stimulusCueStartedTime;
        trajectory_c = allBlocks{b}.inputSensorPositions(startIdx_SC:endIdx_SC);% - allBlocks{b}.trial(t).interactiveZeroInputPos;
        
        trajectory_c = trajectory_c - mean(allBlocks{b}.inputSensorPositions(startIdx_SC:zeroIdx_SC));
        
        trajectory_SC(:,tInd) = interp1(trajectoryTimes_c, trajectory_c, trajectoryTimes);
        
        
        %% calculate wheel trajectory, with respect to interactive start
        [~,startIdx_IS] = min(abs(allBlocks{b}.inputSensorPositionTimes - ...
            (allBlocks{b}.trial(t).interactiveStartedTime + twindow(1))));
        [~,endIdx_IS] = min(abs(allBlocks{b}.inputSensorPositionTimes - ...
            (allBlocks{b}.trial(t).interactiveStartedTime + twindow(2))));
%         [~,zeroIdx] = min(abs(allBlocks{b}.inputSensorPositionTimes - ...
%             (allBlocks{b}.trial(t).interactiveStartedTime)));
        
        trajectoryTimes_c = allBlocks{b}.inputSensorPositionTimes(startIdx_IS:endIdx_IS)- allBlocks{b}.trial(t).interactiveStartedTime;
        trajectory_c = allBlocks{b}.inputSensorPositions(startIdx_IS:endIdx_IS);% - allBlocks{b}.trial(t).interactiveZeroInputPos;
        
        trajectory_c = trajectory_c - mean(allBlocks{b}.inputSensorPositions(startIdx_SC:zeroIdx_SC));%not use SC rather than IS
        
        trajectory_IS(:,tInd) = interp1(trajectoryTimes_c, trajectory_c, trajectoryTimes);
        
        tInd = tInd+1;
    end
end

%% Calculate psychometric data
% rwdTypes = unique(rewardSize(rewardSize>0)); %Array of unique reward types
% rwd = rewardSize';
% perf = (perf / histc(repeatNum,1))*100; %Performance = number of correct / frequency of repeatNum 1 * 100
perf = sum(feedback==1&repeatNum==1)/sum(repeatNum==1)*100;
respTypes = unique(resp(resp>0));
numRespTypes = numel(respTypes);

cDiff = diff(contrast);
cVals = unique(cDiff);

psychoM = zeros(numRespTypes,length(cVals));
psychoMCI = zeros(numRespTypes,length(cVals));
meanRTs = zeros(numRespTypes,length(cVals));
meanRTsCIlow = zeros(numRespTypes,length(cVals));
meanRTsCIup = zeros(numRespTypes,length(cVals));
numTrials = zeros(1,length(cVals));
numChooseR = zeros(numRespTypes, length(cVals));
for r = 1:numRespTypes
    for c = 1:length(cVals) %For the number of unique contrasts
        %         incl_lsrOn = repeatNum==1&cDiff==cVals(c)&rewardSize(2,:)>0; %Logical array of trials that aren't repeats of each contrast
        %         incl_lsrOff = repeatNum==1&cDiff==cVals(c)&rewardSize(1,:)==0; %Logical array of trials that aren't repeats of each contrast
        incl = repeatNum==1&cDiff==cVals(c); %Logical array of trials that aren't repeats of each contrast
        numTrials(c) = sum(incl); %Number of trails is equal to number of non-repeat trials
        numChooseR(r,c) = sum(resp==respTypes(r)&incl);
        
        psychoM(r, c) = numChooseR(r,c)/numTrials(c);
        psychoMCI(r, c) = 1.96*sqrt(psychoM(r, c)*(1-psychoM(r, c))/numTrials(c));
        
        %         meanRTs(r, c) = mean(rt(resp==respTypes(r)&incl));
        %         meanRTsCIlow(r, c) = std(rt(resp==respTypes(r)&incl))./sqrt(sum(resp==respTypes(r)&incl));
        %         meanRTsCIup(r, c) = std(rt(resp==respTypes(r)&incl))./sqrt(sum(resp==respTypes(r)&incl));
        
        q = quantile(rt(resp==respTypes(r)&incl), 3);
        meanRTs(r, c) = q(2);
        meanRTsCIlow(r, c) = q(2)-q(1);
        meanRTsCIup(r, c) = q(3)-q(2);
        
        
    end
end


%% plotting trajectory


%% Plotting curve
if ~isempty(figDestination)
    colors(1,:) = [1 0.8 0]; % Yellow (Stim on left, turns right - correct resp)
    colors(2,:) = [1 0 1]; % Magenta (Stim on right, turns left - correct resp)
    colors(3,:) = [0.2 0.2 0.2];
    %     labels = {'turn right','turn left'};
    
    expRefMod = allBlocks{1}.expRef; expRefMod(expRefMod=='_') = '-';
    scrsz = get(0, 'screensize');
    f(1) = figure('Name', expRefMod, 'NumberTitle', 'Off', 'position',scrsz);
    
    for r = 1:numRespTypes
        subplot(3,1,1);
        %     plot(cVals, psychoM(r, :), 'o', 'Color', colors(r,:));
        %     hold on;
        %     plot(reshape([cVals; cVals; nan(size(cVals))],1,3*length(cVals)), ...
        %         reshape([psychoM(r,:)-psychoMCI(r,:); psychoM(r,:)+psychoMCI(r,:); nan(size(psychoM(r,:)))], 1, 3*length(psychoM(r,:))), 'Color', colors(r,:));
        plotWithErr(cVals, psychoM(r,:), psychoMCI(r,:), colors(r,:));
        hold on;
        plot(cVals, psychoM(r,:), 'ko');
        
        %       Reation time plot
        subplot(3,1,2);
        plot(cVals, meanRTs(r,:), 'o', 'Color', colors(r,:));
        hold on;
        plot(reshape([cVals; cVals; nan(size(cVals))],1,3*length(cVals)), ...
            reshape([meanRTs(r,:)-meanRTsCIlow(r,:); meanRTs(r,:)+meanRTsCIup(r,:); nan(size(meanRTs(r,:)))], 1, 3*length(meanRTs(r,:))), 'Color', colors(r,:));
        xlim([cVals(1) cVals(end)]*1.1);
        xlabel('contrast');
        ylabel('reaction time (sec)');
        
    end
    
    subplot(3,1,1);
    %     legend({'Left' 'Right'}, 'Location', 'EastOutside');
    plot([cVals(1) cVals(end)], [0.5 0.5], 'k:');
    ylim([0 1]);
    xlim([cVals(1) cVals(end)]*1.1);
    xlabel('contrast');
    ylabel('proportion choose R');
    title([expRefMod ', numTrials = ' num2str(numCompletedTrials) ', perf = ' num2str(perf,3) '%']);
    
    subplot(3,1,3);
    plot(cVals, numTrials, 'ko');
    xlim([cVals(1) cVals(end)]*1.1);
    ylabel('number of trials')
    xlabel('contrast');
    yl = ylim(); ylim([0 yl(2)]);
end

if length(unique(resp))>2&&plotFit==true
    plotFit = 2;
elseif length(unique(resp))<3&&plotFit==true
    plotFit = 1;
end

%% Plotting the fit
switch plotFit
    case 1 % PsychoFit
        nn = numTrials;
        xx = cVals*100; %Unique contrasts in percent
        if size(psychoM,1)>1
            pp = psychoM(2,:); %Performance
        else
            pp = psychoM;
        end
        
        if any(xx==0)
            intercept = psychoM(xx==0);
        else
            intercept=0;
            %             intercept = interp1(pp,xx, 0.5, 'spline');
        end
        
        %Calculate params for line fit & pass to Psychofit
        parstart = [intercept mean(nn) min(pp)];
        parmin = [-40 10 min(pp)];
        parmax = [40 max(nn) max(pp)];
        pars = mle_fit_psycho([xx;nn;pp],'erf_psycho',parstart,parmin,parmax,25); % Matteo's function
        %     pars = mle_fit_psycho([cVals;nn;pp],'erf_psycho_2gammas',parstart,parmin,parmax,1);
        
        
        %Plot the data
        f(2) = figure('Name', [expRefMod '_fit'], 'NumberTitle', 'Off', 'Units', 'Normalized', 'Position', [0.3 0.6 0.50 0.30]);
        plot([0 0],[0 1],'k:');
        hold on
        plot(xx,pp,'ko','markerfacecolor','k');  % Plot points
        plot(-100:100, erf_psycho( pars, -100:100 ),'b','LineWidth',1.1); % Plot the fit
        plot([-100 100], [0.5 0.5], 'k:'); % Extend axis limit
        xlabel('contrast (%)', 'FontWeight', 'bold');
        ylabel('% Right', 'FontWeight', 'bold');
        hold off
        title(sprintf('b = %2.0f, t = %2.0f, l = %0.2f, n = %2.0f',pars, sum(nn)));
    case 2 % GLM
        if exist('GLM','file')==2||exist('GLM','builtin')==2
            data.contrast_cond = contrast';
            data.response = resp';
            data.repeatNum = repeatNum';
            if any(all(contrast))
                modelString = 'C^N';
            else
                modelString = 'C^N-subset';
            end
            g = GLM(data);
            g = g.setModel(modelString);
            g = g.fit;
            f(2) = figure('Name', [expRefMod '_' modelString '-fit'], 'NumberTitle', 'Off', 'Units', 'Normalized', 'Position', [0.3 0.6 0.50 0.30]);
            g.plotData;g.plotFit;
        else
            warning('Could not plot fit; GLM toolbox not found');
        end
end

%% plot wheel trajectory and probability
f(2) = figure('position',scrsz);
for c = 1:length(cVals) %For the number of unique contrasts
    for r = 1:numRespTypes
        incl = repeatNum==1&cDiff==cVals(c); %Logical array of trials that aren't repeats of each contrast
        
        trIdx = find(resp==respTypes(r)&incl);
        
        if isempty(trIdx) continue;
        else
            subplot(3,length(cVals),c);
            switch r
                case 1
                    lc = 'r';
                case 2
                    lc = 'b';
                case 3
                    lc = 'g';
            end
            
            
            plot(trajectoryTimes, trajectory_SC(:,trIdx)','color',lc);
            xlim([trajectoryTimes(1) trajectoryTimes(end)]);
            ylim([-500 500]);
            grid on;
            hold on
            if c == 1
                title([expRefMod ' contrast:' num2str(cVals(c))])
            else
                title(['contrast:' num2str(cVals(c))]);
            end
            if c == 1 && r == 1
                xlabel('time from stimCueOn');
            end
            
            subplot(4,length(cVals),c+r*length(cVals));
            probDensity = zeros(size(trajectory_SC,1),51);
            trajectoryPosition = linspace(-500,500,51);
            for tt=1:size(trajectory_SC,1)
                probDensity(tt,:) = hist(trajectory_SC(tt,trIdx),trajectoryPosition);
            end
            probDensity = probDensity/length(trIdx);
            imagesc(trajectoryTimes, trajectoryPosition, probDensity');
            grid on
            axis xy;
            caxis([0 1]);
            
            if c==1;
                title('Probability density');
            end
        end
        
    end
end
xlabel('time from cue onset [s]');
ylabel('wheel position');
%legend('left choice','right choice');


%% plot wheel trajectory and probability
f(3) = figure('position',scrsz);
for c = 1:length(cVals) %For the number of unique contrasts
    for r = 1:numRespTypes
        incl = repeatNum==1&cDiff==cVals(c); %Logical array of trials that aren't repeats of each contrast
        
        trIdx = find(resp==respTypes(r)&incl);
        
        if isempty(trIdx) continue;
        else
            subplot(3,length(cVals),c);
            switch r
                case 1
                    lc = 'r';
                case 2
                    lc = 'b';
                case 3
                    lc = 'g';
            end
            
            
            plot(trajectoryTimes, trajectory_IS(:,trIdx)','color',lc);
            xlim([trajectoryTimes(1) trajectoryTimes(end)]);
            ylim([-500 500]);
            grid on;
            hold on
            if c == 1
                title([expRefMod ' contrast:' num2str(cVals(c))])
            else
                title(['contrast:' num2str(cVals(c))]);
            end
            if c == 1 && r == 1
                xlabel('time from interactiveStart');
            end
            
            subplot(4,length(cVals),c+r*length(cVals));
            probDensity = zeros(size(trajectory_IS,1),51);
            trajectoryPosition = linspace(-500,500,51);
            for tt=1:size(trajectory_IS,1)
                probDensity(tt,:) = hist(trajectory_IS(tt,trIdx),trajectoryPosition);
            end
            probDensity = probDensity/length(trIdx);
            imagesc(trajectoryTimes, trajectoryPosition, probDensity');
            grid on
            axis xy;
            caxis([0 1]);
            
            if c==1;
                title('Probability density');
            end
        end
        
    end
end
xlabel('time from interactive onset [s]');
ylabel('wheel position');

%% Save Plots
figName = fullfile(figDestination, [allBlocks{1}.expRef '_psychometric']);
type = 'png';
if ~any(strcmp(figDestination, {'none' ''}));
    if ~exist(figDestination, 'dir')
        mkdir(figDestination);
    end
    
    switch type
        case 'fig'
            savefig(f,[figName '.fig']);
        case {'jpg' 'png' 'jpeg'}
            if numel(f)>1
                for i = 1:length(f)
                    saveSameSize(f(i),'file',[figName '(' num2str(i) ')'],'format',type)
                end
            else
                saveSameSize(f,'file',figName,'format',type)
            end
        otherwise
            error('Format not supported.  Figure not saved.');
    end
end

%% Save data
psychoDat.expRef = expRef;
psychoDat.contrast = contrast;
psychoDat.resp = resp;
psychoDat.feedback = feedback;
psychoDat.repeatNum = repeatNum;
psychoDat.rt = rt;
psychoDat.rewardSize = rewardSize;
