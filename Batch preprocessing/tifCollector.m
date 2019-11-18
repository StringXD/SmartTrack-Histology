function[] = tifCollector()
% Put the preprocessed tiff image for each mouse together into a new folder
% and sort
% folder name should be like: 'Npxl_MouseNum_SlideNum'
scans = uigetdir2;
for i = 1:length(scans)
    tempath = [scans{i} '\'];
    cd(tempath);
    slicefolders = dir('Npxl*');
    mouseIDstr = cell(length(slicefolders),1);
    for j = 1:length(slicefolders)
        C = strsplit(slicefolders(j).name,'_');
        mouseIDstr{j} = C{2};
    end
    mouseID = unique(mouseIDstr);
    for id = 1:length(mouseID)
        mkdir([mouseID{id} '_indexed']);
        tifNum = 0;
        currMouseDirs = dir(['Npxl_' mouseID{id} '*']);
        currMouseDirNames = {currMouseDirs.name};
        [currMouseDirNames,~] = sort_nat(currMouseDirNames,'ascend');
        for sliceID = 1:length(currMouseDirs)
            tempSliceFolder = [scans{i} '\' currMouseDirNames{sliceID} '\'];
            cd(tempSliceFolder);
            images = dir('*.tif');
            imgNames = {images.name};
            [imgNames,~] = sort_nat(imgNames,'ascend');
            for img = 1:length(imgNames)
                tifNum = tifNum + 1;
                copyfile(imgNames{img},[scans{i} '\' mouseID{id} '_indexed\' num2str(tifNum) '.tif']);
            end
            cd ..
        end
    end
end
                
                
            
    end
    
    
        
    