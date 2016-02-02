
cwExp = 29;
useWheel = true;
suffix_out = 'ar gd z';


Exps = readExpsDatabase('ExpsDatabase_cw.m', cwExp);
    
    
[p, block] = ProtocolLoadDS(Exps);

% contrast_cache = [];
% for itr = 1:block.numCompletedTrials
%     contrast_cache = [contrast_cache block.trial(itr).condition.visCueContrast];
% end
% contrastList = unique(contrast_cache);

%for expid 24
% cueList = {'l','l','l','r','r','r'};
% respList = {'r','r','r','r','r','r'};%{'l','l','l','l','l','l'};
% contrastList = [0.5 0.25 0.12 0.12 0.25 0.5];

%for expid 29
cueList = {'l','l','l','l','r','r','r','r'};
respList = {'','','','','','','',''};%{'r','r','r','r','r','r'};%{'l','l','l','l','l','l'};
contrastList = [0.5 0.25 0.12 0.06 0.06 0.12 0.25 0.5];


MyStackDir = tools.getDirectory( localDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_out);



iii = 1;
Sall = StackSet;
for nc = 1:length(cueList)
    
    cue = cueList{nc};
    resp = respList{nc};
    stimContrast = contrastList(nc);
    
    disp(['Cue: ' cue ', Resp: ' resp ', Contrast: ' num2str(stimContrast)]);
    
    
    ThisFileName_v = sprintf('avgEvent_v_%s_%s_%s_%s', evoke, cue, resp, num2str(100*stimContrast));
    if ~useWheel
        ThisFileName_v = [ThisFileName_v '_all'];
    end
    if exist(fullfile(MyStackDir, [ThisFileName_v '.mat']))
        S = tools.LoadMyStacks(MyStackDir, ThisFileName_v);
        
        S = S.SubtractBase(-0.1,0);
        [~,~,~,nTrials] = strread(S.Description, '%s %s %s %d');
        
        if isempty(Sall.Values)
            Sall = S;
            Sall.Description = [];
        else
            Sall = Sall.AssignCondition(S.Values,iii);
        end
        Sall.Description{iii} = ['Cue: ' cue ', Resp: ' resp ...
            ', Contrast: ' num2str(stimContrast) ', nTrials:' num2str(nTrials)];
        iii=iii+1;
    end
end

%% save stacks
saveDir = '\\zserver\Lab\Share\Shimaoka';
saveStackDir = tools.getDirectory( saveDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_out);
ThisFileName = sprintf('avgEvent_v_stimCue__%s__all', respList{1});
tools.SaveMyStacks(Sall, saveStackDir, [], ThisFileName);

