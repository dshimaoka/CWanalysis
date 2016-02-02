function [allBlocks, psychoDat] = multiPsycho_withTrajectory(subject, dates, dateRange, ...
    separateBlocks, plotFit, makeFig)
% multiPsycho gathers the blocks and (optionally) plots the psychometic data
% from the sessions of the subject and dates provided.
%
% subject           string; Name of mouse, e.g. 'M150326_TEST'
% dates             string / cell array of strings / vector of datenums; 
%                   date of experiment(s), e.g. {'2015-03-20' '2015-03-26'}
%                   default calls last expriment date; 'last'
% (dateRange)       logical (optional); if true, all sessions between the 
%                   dates are gathered
% (separateBlocks)  logical (optional); if true each training session on
%                   the same day is analysed separately, if false all
%                   training sessions are collapsed (default)
% (plotFit)         logical (optional); if true, plot with PsychoFit (2AFC) 
%                   or GLM (2AUC, requires GLM toolbox)
% (makeFig)         logical (optional); if true psychometric curves are
%                   saved to \\zserver\Data\behavior\ChoiceWorld\ (default)

addpath('\\zserver\Code\PsychoGen');
addpath('\\zserver\Code\Rigging\main');
addpath('\\zserver\Code\Rigging\cortexlab');%dat.mpepMessageParse

%% Initialize inputs
if nargin < 1
    error('Error in multiPsycho: Must specify at least the subject');
end
if nargin < 2
   dates = 'last';
end
if nargin < 3
   dateRange = false;
end
if nargin < 4
   separateBlocks = false;
end
if nargin < 5
   plotFit = false;
end
if nargin < 6
   makeFig = true;
end    
% define path to save figure
if makeFig 
    figDestination = ['\\zserver\Data\behavior\ChoiceWorld\' subject '\'];
else
    figDestination = 'none';
end
% get list of references and dates for subject
[expRef, expDate] = dat.listExps(subject);
% convert the date(s) to vector of datenums
if ~isa(dates,'numeric')
    if strcmpi(dates,'last')
        dates = max(expDate);
    else
        dates = datenum(dates,'yyyy-mm-dd');
    end
elseif isa(dates,'cell')
    dates = cell2mat(dates);
end
% get date nums between specified range
if numel(dates)==2&&dateRange==1
    dates = sort(dates);
    dates = dates(1):dates(2);
end
if size(dates,2)>size(dates,1)
    dates = dates';
end
% find paths to existing block files
idx = cell2mat(arrayfun(@(x)find(expDate==x),dates,'UniformOutput',0));
filelist = mapToCell(@(r) dat.expFilePath(r, 'block', 'master'),expRef(idx));
existfiles = cell2mat(mapToCell(@(l) file.exists(l),filelist));
% load block files
allBlocks = cellfun(@load,filelist(existfiles));

%% Pass to generatePsychometric
expDate = expDate(idx);
[~,~,ic] = unique(expDate(existfiles));
psychoDat = struct(...
    'expRef',{''},...
    'contrast',[],...
    'resp',[],...
    'feedback',[],...
    'repeatNum',[],...
    'rt',[],...
    'rewardSize',[]);
i=1;
for j = 1:max(ic)
    if separateBlocks
         for k = 1:length(allBlocks(ic==j))
            psychoDat(i) = generateWheelTrajectory({allBlocks(k).block}, figDestination, plotFit);
            i=i+1;
         end
    else
        psychoDat(j) = generateWheelTrajectory({allBlocks(ic==j).block}, figDestination, plotFit);
        i=i+1;
    end
end

%% Option to save psychometric data
% if ~isempty([psychoDat(:).resp])
%     str = input('Do you want to save psychometric data? Y/N [N]: ','s');
%     if isempty(str)||strcmpi(str,'n')
%         clear psychoDat;
%     else
%         [filename, path] = uiputfile([subject '_psychoDat.m'],'Save psychometric data');
%         save([path filename],'psychoDat');
%     end
% end
