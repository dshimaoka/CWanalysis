
tgttrials = [];
for cues = 1:2
    for resps = 1:2
        for contrasts = 1:4
            
            if cues == 1;
                cue = 'l';
            elseif cues == 2
                cue = 'r';
            end
            
            if resps == 1;
                resp = 'l';
            elseif resps == 2
                resp = 'r';
            end
            
            if contrasts == 1;
                contrast = 1;
            elseif contrasts == 2
                contrast = 0.5;
            elseif contrasts == 3
                contrast = 0.2;
            elseif contrasts == 4
                contrast = 0.1;
            end
            
            tgttrials(cues,resps,contrasts) = ...
            length(selecttrials(block.trial, cue, resp, '', 'nrep', contrast));
            
        end
    end
end

for contrasts = 1:4
    hitrate_cueleft(contrasts) = tgttrials(1,2,contrasts) / sum(tgttrials(1,:,contrasts));
    hitrate_cueright(contrasts) = tgttrials(2,1,contrasts) / sum(tgttrials(2,:,contrasts));
end
    
    