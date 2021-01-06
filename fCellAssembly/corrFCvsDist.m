% correlation between FC density and distance
%% remap the su paths
load 114_sorted_file_path.mat
homedir='/home/zx/neupix/wyt';
path_list=h5read('transient_6.hdf5','/path');

metaFolders = cell(length(path_list),1);
sorted_meta_fpath = cell(length(sorted_fpath),1);

for i = 1:length(path_list)
    folder=regexp(path_list{i},'(\w|\\|-)*','match','once');
    [~,~,~,metaFolders{i},~]=jointFolder(folder,cell(0),homedir);
end

for j = 1:length(sorted_fpath)
    folder=regexp(sorted_fpath{j},'(\w|\\|-)*','match','once');
    [~,~,~,sorted_meta_fpath{j},~]=jointFolder(folder,cell(0),homedir);
end

save('remappedPaths.mat','metaFolders','sorted_meta_fpath');

remapIdx = cellfun(@(x) find(startsWith(string(sorted_meta_fpath),x)),metaFolders,'UniformOutput',false);
abnormalpaths = unique(metaFolders(cellfun(@(x) length(x) ~= 1,remapIdx)));

%% translate cluster_id to uniqueId and link to coordinates
load sucoords.mat
cid_list=h5read('transient_6.hdf5','/cluster_id');

pathOK = cellfun(@(x) length(x) == 1,remapIdx);
probe2ndLabel = contains(path_list,'imec1') & ~contains(metaFolders,'singleProbe');
uid_list = cell2mat(remapIdx(pathOK))*100000 + 10000*probe2ndLabel(pathOK) + double(cid_list(pathOK));
coord = coord(pathOK,:);

%% construct inter-neuron distance dataset
distPairSum = [];
for sessIdx = 1:length(sorted_fpath)
    currCidSel = uid_list>sessIdx*100000 & uid_list<(sessIdx+1)*100000;
    currCids = uid_list(currCidSel);
    currCoords = coord(currCidSel,:);
    % use combination trick to faster computing
    currCidIndices = 1:length(currCids);
    C = nchoosek(currCidIndices,2);
    coordDiffs = num2cell(currCoords(C(:,1),:) - currCoords(C(:,2),:),2);
    dists = cellfun(@norm,coordDiffs);
    distPairSum = [distPairSum;currCids(C(:,1)),currCids(C(:,2)),dists];
end
distPairSum = [sort(distPairSum(:,1:2),2),distPairSum(:,3)];
save('distPairSum.mat','distPairSum','-v7.3');

%% construct inter-neuron diatance dataset (1 probe)
distPairSum = [];
validOriPaths = path_list(pathOK);
uniPaths = unique(validOriPaths);
for trackIdx = 1:length(uniPaths)
    currCidSel = contains(validOriPaths,uniPaths{trackIdx});
    currCids = uid_list(currCidSel);
    currCoords = coord(currCidSel,:);
    % use combination trick to faster computing
    currCidIndices = 1:length(currCids);
    C = nchoosek(currCidIndices,2);
    coordDiffs = num2cell(currCoords(C(:,1),:) - currCoords(C(:,2),:),2);
    dists = cellfun(@norm,coordDiffs);
    distPairSum = [distPairSum;currCids(C(:,1)),currCids(C(:,2)),dists];
end

%% corr FC with dist
for bin=1:6
    load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    pair_regs{bin}=pair_reg;
    reg_chains_S1{bin}=reg_chain_S1;
    reg_chains_S2{bin}=reg_chain_S2;
    conn_chains_S1{bin}=conn_chain_S1;
    conn_chains_S2{bin}=conn_chain_S2;
    pair_chains{bin}=pair_chain;
    clear pair_reg reg_chain_S1 reg_chain_S2 pair_chain conn_chain_S1 conn_chain_S2
end
edges = 0:0.5:3.5;
thres = 100;
mns = cell(length(edges)-1,1);
mmn = zeros(length(edges)-1,1);
stdn = zeros(length(edges)-1,1);
for distIdx = 1:length(edges)-1
    conncount = nan(length(sorted_fpath),6,2);
    pairscount = nan(length(sorted_fpath),1);
    for ses = 1:length(sorted_fpath)
        locSel = distPairSum(:,1)>ses*100000 & distPairSum(:,1)<(ses+1)*100000 & distPairSum(:,3)>edges(distIdx)*100 & distPairSum(:,3)<=edges(distIdx+1)*100;
        locPairs = distPairSum(locSel,1:2);
        pairsel = pair_chains{1}(:,1)>ses*100000 & pair_chains{1}(:,1)<(ses+1)*100000;
        pairDiff = diff(pair_regs{1},1,2) & pairsel;
        currPairs = sort(pair_chains{1}(pairDiff,:),2);
        validPairs = intersect(locPairs,currPairs,'rows');
        pairscount(ses) = length(validPairs);
        for bin=1:6
            selS1 = conn_chains_S1{bin}(:,1)>ses*100000 & conn_chains_S1{bin}(:,1)<(ses+1)*100000;
            selS2 = conn_chains_S2{bin}(:,1)>ses*100000 & conn_chains_S2{bin}(:,1)<(ses+1)*100000;
            connDiffS1 = diff(reg_chains_S1{bin},1,2) & selS1;
            connDiffS2 = diff(reg_chains_S2{bin},1,2) & selS2;
            currConnS1 = sort(conn_chains_S1{bin}(connDiffS1,:),2);
            currConnS2 = sort(conn_chains_S2{bin}(connDiffS2,:),2);
            validConnS1 = intersect(locPairs,currConnS1,'rows');
            validConnS2 = intersect(locPairs,currConnS2,'rows');
            conncount(ses,bin,1) = size(validConnS1,1);
            conncount(ses,bin,2) = size(validConnS2,1);
        end
    end
    mn = [];
    for ses = 1:length(sorted_fpath)  
        sescountcba = [];
        if pairscount(ses) > thres
            for bin = 1:6
                sescountcba=[sescountcba;squeeze(conncount(ses,bin,:)/pairscount(ses))];
            end
            mn=[mn;ses,mean(sescountcba)];
        end
    end
    mns{distIdx} = mn;
    mmn(distIdx) = mean(mn(:,2));
    stdn(distIdx) = std(mn(:,2));
end
%% plot
errorbar(0.5:0.5:3.5,mmn,stdn/sqrt(length(sorted_fpath)-1));
xlabel('Distance(mm)');
ylabel('Functional coupling density');

%% ANOVA
for distIdx = 1:length(edges)-1
    mns{distIdx} = [mns{distIdx},distIdx*ones(length(mns{distIdx}),1)];
end
mnss = cell2mat(mns)';
groups = string(num2cell(mnss(3,:)*0.5));
[p,tbl,stats] = anova1(mnss(2,:),groups);    













%%
function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list,homedir)
metaFolder=replace(folder,'\','/');
metaFolder=fullfile(fullfile(homedir,'DataSum'),metaFolder);
if isfolder(metaFolder)
    spkFolder=replace(metaFolder,'imec1','imec0');
    file=dir(fullfile(spkFolder,'spike_info.mat'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 2-tracks');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        %             pause;
        return
    end
    folderType=2;
    C = strsplit(metaFolder,'/');
    metaFolder = strjoin(C(1:end-1),'/');
else
    metaFolder=replace(metaFolder,'DataSum','DataSum/singleProbe');
    spkFolder=metaFolder;
    file=dir(fullfile(spkFolder,'spike_times.npy'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 1-track');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        return
    end
    folderType=1;
end
end