rawDir = '\\zserver2.ioo.ucl.ac.uk\Data\GCAMP';
%Exps = readExpsDatabase('ExpsDatabase_cw.m', 14);
Exps.animal = 'Affogato52';
Exps.iseries = '20150402';
Exps.iexp = 2;
p = ProtocolLoadDS(Exps,1);
Exps.Cam.FileString = '_cam2';

for iStim = 1
    directory = fullfile(rawDir, [Exps.animal Exps.Cam.FileString], num2str(Exps.iseries), num2str(Exps.iexp));
    directory_cache = [directory '_cache'];
    
    if ~exist(directory_cache, 'dir')
        mkdir(directory_cache)
    end
    
    names = dir(directory);
    names = {names(~[names.isdir]).name};
    
    nFiles = numel(names);
    
    for nn = 1:nFiles
        
        
        [~, name_nosuffix, ext] = fileparts(names{nn});
        
        name_separate = textscan(name_nosuffix, '%s %s %d','delimiter','_');
        
        %file number created by camware
        fileNumber = name_separate{3};
        
        newName = [Exps.animal Exps.Cam.FileString '_' num2str(fileNumber+1) ext];
        
        Source = fullfile(directory, names{nn});
        Dest = fullfile(directory_cache, newName);
       
        [Status(nn), Msg{nn}] = FileRename(Source, Dest);
    end
    
end

if ~any(Status-1)
    rmdir(directory);
    movefile(directory_cache, directory);
end
 
    

