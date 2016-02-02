
shareDir = '\\zserver\Lab\Share\Shimaoka';

Exps = readExpsDatabase('ExpsDatabase_cw.m', 24);
[p, block] = ProtocolLoadDS(Exps);
iStim = 1;
suffix_out = 'ar gd';

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
    
    MyStackDir = tools.getDirectory( shareDir, p, Exps.ResizeFac, iStim, 1, Exps.Cam, suffix_out);
    
    % MyStackDir = '\\zserver\Lab\Share\Shimaoka\65_ratio\20151012\005\Resize 25 Stim 001 Repeat All ar gd z';
    
    %load avg movie triggered by movement onset
    ThisFileName_v = 'avgEvent_v_moveOn____all';
    Smov = tools.LoadMyStacks(MyStackDir, ThisFileName_v);
    
    if cams==1
        Smov_a = Smov;
    else
        Smov_a = Smov_a.AssignCondition(Smov.Values, cams);
    end
end

% %load avg movie triggered by stimCue. choose right
% ThisFileName_v = 'avgEvent_v_stimCue__r__all';
% Sr = tools.LoadMyStacks(MyStackDir, ThisFileName_v);
%
% %load avg movie triggered by stimCue. choose left
% ThisFileName_v = 'avgEvent_v_stimCue__l__all';
% Sl = tools.LoadMyStacks(MyStackDir, ThisFileName_v);


