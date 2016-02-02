%1st induction of (0,0)
%load('\\zserver\Data\expInfo\Zorro62\2015-07-20\2\2015-07-20_2_Zorro62_Block.mat');

%last time animal did a lot, but biased. Cone is in place
%load('\\zserver\Data\expInfo\Zorro62\2015-07-17\3\2015-07-17_3_Zorro62_Block.mat');

resptime = '';
nrep = 'nrep';

contrastList_c = [];
for bb = 1:length(block.trial)
    contrastList_c = [contrastList_c [block.trial(bb).condition.visCueContrast]];
end

contrastList = unique(contrastList_c)';

numTrials = [];
fractionTrials = [];
for cc = 1:3 %stimulus side
    for rr = 1:3 %response side
        for cidx = 1:length(contrastList)
            stimContrast = contrastList(cidx);
            
            
            if cc==1
                cue='l';
            elseif cc==2
                cue='r';
            elseif cc==3
                cue='o';
            end
            
            if rr==3
                resp = 'o';
            elseif rr==1
                resp='l';
            elseif rr==2
                resp='r';
            end
            
            tgtTr = selecttrials(block.trial, cue, resp,resptime,nrep,stimContrast);
            numTrials(cc,rr,cidx) = length(tgtTr);
            tgtTr = selecttrials(block.trial, cue, [],resptime,nrep,stimContrast);
            fractionTrials(cc,rr,cidx) = numTrials(cc,rr,cidx)/length(tgtTr);
        end
    end
end

numTrials_conc = cat(2, fliplr(squeeze(numTrials(1,:,2:length(contrastList)))), squeeze(numTrials(3,:,1))', squeeze(numTrials(2,:,2:length(contrastList))));
fractionTrials_conc = cat(2, fliplr(squeeze(fractionTrials(1,:,2:length(contrastList)))), squeeze(fractionTrials(3,:,1))', squeeze(fractionTrials(2,:,2:length(contrastList))));
contrastList_conc = cat(2, -fliplr(contrastList(2:length(contrastList))), 0, contrastList(2:length(contrastList)));

subplot(121);
plot(contrastList_conc, numTrials_conc);

subplot(122);
plot(contrastList_conc, fractionTrials_conc);
legend('choose R', 'choose L' ,'no Go','location','northoutside');

dt = between(datetime(block.startDateTimeStr), datetime(block.endDateTimeStr));
title([num2str(block.numCompletedTrials) 'trials in ' char(dt)]);

