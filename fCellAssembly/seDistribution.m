% finer analysis for specify start/end cell
%standard : at least 1 showup in one trial in average
% load rings.mat
% load ringSampleTrials.mat

% ringsm = cell(1,3);
% metadata = cell(1,3);
% samples = cell(1,3);


% for midx = 1:3;
% r1=[];
% for bb=1:6
%     r1=[r1;cell2mat(rings(midx,:,bb,1)')];
% end
% r2=[];
% for bb=1:6
%     r2=[r2;cell2mat(rings(midx,:,bb,2)')];
% end
% r1(:,midx+3)=1;r2(:,midx+3)=2;ringsm{midx}=[r1(:,1:(midx+2));r2(:,1:(midx+2))];samps=[r1(:,midx+3);r2(:,midx+3)];
% samples{midx} = samps;
% metadata{midx} = zeros(length(ringsm{midx}),midx+3);
% for ssidx=1:length(ringsm{midx})
%     sessIdx=idivide(int32(ringsm{midx}(ssidx,1)),int32(100000));
%     transIds=int32(rem(ringsm{midx}(ssidx,:),100000))';
%     metadata{midx}(ssidx,1) = sessIdx;
%     metadata{midx}(ssidx,2:end) = transIds';
% end
% end

% %% 1 cell for all rings
% seSumLoose = cell(length(startEndStat),1);
% seSumStrict = cell(length(startEndStat),1);
% for ridx = 1:length(startEndStat)
%     seSumLoose{ridx} = cell2mat(startEndStat{ridx,1});
%     seSumStrict{ridx} = cell2mat(startEndStat{ridx,2});
% end
% transData = metadata(:,2:end);
% uniqueSuList = unique(transData);
% initRate = zeros(length(uniqueSuList),2); % col1 for loose; col2 for strict
% endRate = zeros(length(uniqueSuList),2);
% for i = 1:length(uniqueSuList)
%     [row,col] = find(transData == uniqueSuList(i));
%     tempLS = [];
%     tempSS = [];
%     tempLE = [];
%     tempSE = [];
%     for j = 1:length(row)
%         if ~isempty(seSumLoose{row(j)})
%             tempLS = [tempLS;seSumLoose{row(j)}(:,1) == col(j)];
%             tempLE = [tempLE;seSumLoose{row(j)}(:,2) == col(j)];
%         end
%         if ~isempty(seSumStrict{row(j)})
%             tempSS = [tempSS;seSumStrict{row(j)}(:,1) == col(j)];
%             tempSE = [tempSE;seSumStrict{row(j)}(:,2) == col(j)];
%         end
%     end
%     if length(tempLS) < 5
%         initRate(i,1) = -1;
%     else
%         initRate(i,1) = sum(tempLS)/length(tempLS);
%     end
%     if length(tempSS) < 5
%         initRate(i,2) = -1;
%     else
%         initRate(i,2) = sum(tempSS)/length(tempSS);
%     end
%     if length(tempLE) < 5
%         endRate(i,1) = -1;
%     else
%         endRate(i,1) = sum(tempLE)/length(tempLE);
%     end
%     if length(tempSE) < 5
%         endRate(i,2) = -1;
%     else
%         endRate(i,2) = sum(tempSE)/length(tempSE);
%     end
% end
% edges = 0:0.02:1;
% validInitRate = initRate(initRate(:,2)~=-1,2);
% [counts,centers] = hist(validInitRate,edges);
% bar(centers,counts/length(validInitRate));
% ylabel('Probability');
% xlabel('Fraction of times the neuron initiate the ring');
% hist(initRate(:,2),50)
% set(gca,'XLim',[0,1]);
% hist(endRate(:,2),50)
% set(gca,'XLim',[0,1]);
% %% 1 cell for 1 ring
%% plot
load over1HzRingData.mat
initRateSum = [];
for ssidx = 1:length(over1HzRingData)
    currData = over1HzRingData{ssidx};
    if isempty(currData)
        continue;
    end
    ringSuData = currData(:,7:11);
    uniqueSus = unique(ringSuData);
    uniqueSus(uniqueSus == -1) = [];
    for i = uniqueSus(:)'
        [rows,cols] = find(ringSuData == i);
        initOutcomes = currData(rows,1);
        initCount = sum(initOutcomes == cols);
        initRate = initCount/length(rows);
        initRateSum = [initRateSum;initRate];
    end
end


hist(initRateSum,50)
ylabel('Probability');
xlabel('Initiation Rate');

% shuffled data
shuffledRateSum = [];
for iter = 1:1000
    shuffledData = over1HzRingData;
    initRateSum = [];
    for ssidx = 1:length(shuffledData)
        currData = shuffledData{ssidx};
        if isempty(currData)
            continue;
        end
        ringSuData = currData(:,7:11);
        uniqueSus = unique(ringSuData);
        uniqueSus(uniqueSus == -1) = [];
        rsize = sum(ringSuData ~= -1,2);
        psedoInitLabel = arrayfun(@(x) randi(x),rsize);
        for i = uniqueSus(:)'
            [rows,cols] = find(ringSuData == i);
            initCount = sum(psedoInitLabel(rows) == cols);
            initRate = initCount/length(rows);
            initRateSum = [initRateSum;initRate];
        end
    end
    shuffledRateSum = [shuffledRateSum,initRateSum];
end
meanShuffledData = mean(shuffledRateSum,2);
hist(meanShuffledData)

%% plot loops with different lengths separately
load over1HzRingData.mat

initRateSum3 = [];
initRateSum4 = [];
initRateSum5 = [];

for ssidx = 1:length(over1HzRingData)
    currData = over1HzRingData{ssidx};
    if isempty(currData)
        continue;
    end
    ringSuData = currData(:,7:11);
    uniqueSus = unique(ringSuData);
    uniqueSus(uniqueSus == -1) = [];
    rsize = sum(ringSuData ~= -1,2);
    for i = uniqueSus(:)'
        [rows,cols] = find(ringSuData == i);
        initOutcomes = currData(rows,1);
        initCount3 = sum(initOutcomes == cols & rsize(rows) == 3);
        initRate3 = initCount3/sum(rsize(rows) == 3);
        initCount4 = sum(initOutcomes == cols & rsize(rows) == 4);
        initRate4 = initCount4/sum(rsize(rows) == 4);
        initCount5 = sum(initOutcomes == cols & rsize(rows) == 5);
        initRate5 = initCount5/sum(rsize(rows) == 5);
        initRateSum3 = [initRateSum3;initRate3];
        initRateSum4 = [initRateSum4;initRate4];
        initRateSum5 = [initRateSum5;initRate5];
    end
end
subplot(1,3,1);
hist(initRateSum3)
subplot(1,3,2)
hist(initRateSum4)
subplot(1,3,3)
hist(initRateSum5)

% shuffled data
shuffledRateSum3 = [];
shuffledRateSum4 = [];
shuffledRateSum5 = [];
for iter = 1:1000
    shuffledData = over1HzRingData;
    initRateSum3 = [];
    initRateSum4 = [];
    initRateSum5 = [];
    for ssidx = 1:length(shuffledData)
        currData = shuffledData{ssidx};
        if isempty(currData)
            continue;
        end
        ringSuData = currData(:,7:11);
        uniqueSus = unique(ringSuData);
        uniqueSus(uniqueSus == -1) = [];
        rsize = sum(ringSuData ~= -1,2);
        psedoInitLabel = arrayfun(@(x) randi(x),rsize);
        for i = uniqueSus(:)'
            [rows,cols] = find(ringSuData == i);
            initCount3 = sum(psedoInitLabel(rows) == cols & rsize(rows) == 3);
            initRate3 = initCount3/sum(rsize(rows) == 3);
            initCount4 = sum(psedoInitLabel(rows) == cols & rsize(rows) == 4);
            initRate4 = initCount4/sum(rsize(rows) == 4);
            initCount5 = sum(psedoInitLabel(rows) == cols & rsize(rows) == 5);
            initRate5 = initCount5/sum(rsize(rows) == 5);
            initRateSum3 = [initRateSum3;initRate3];
            initRateSum4 = [initRateSum4;initRate4];
            initRateSum5 = [initRateSum5;initRate5];
        end
    end
    shuffledRateSum3 = [shuffledRateSum3,initRateSum3];
    shuffledRateSum4 = [shuffledRateSum4,initRateSum4];
    shuffledRateSum5 = [shuffledRateSum5,initRateSum5];
end
save('shuffledInitRateSums.mat','shuffledRateSum3','shuffledRateSum4','shuffledRateSum5');
meanShuffledData3 = mean(shuffledRateSum3,2);
meanShuffledData4 = mean(shuffledRateSum4,2);
meanShuffledData5 = mean(shuffledRateSum5,2);
meanShuffledData3(isnan(meanShuffledData3)) = [];
meanShuffledData4(isnan(meanShuffledData4)) = [];
meanShuffledData5(isnan(meanShuffledData5)) = [];
initRateSum3(isnan(initRateSum3)) = [];
initRateSum4(isnan(initRateSum4)) = [];
initRateSum5(isnan(initRateSum5)) = [];


shuffledRateSum3(isnan(shuffledRateSum3(:,1)),:) = [];
shuffledRateSum4(isnan(shuffledRateSum4(:,1)),:) = [];
shuffledRateSum5(isnan(shuffledRateSum5(:,1)),:) = [];

bootci3 = bootci(1000,@mean,shuffledRateSum3');
bootci4 = bootci(1000,@mean,shuffledRateSum4');
bootci5 = bootci(1000,@mean,shuffledRateSum5');

[h1,p1] = kstest2(initRateSum3,meanShuffledData3);
[h2,p2] = kstest2(initRateSum4,meanShuffledData4);
[h3,p3] = kstest2(initRateSum5,meanShuffledData5);

edges3 = 0.3:0.0006:0.36;
edges4 = 0.2:0.001:0.3;
edges5 = 0.15:0.001:0.25;

cdfu3 = histcounts(bootci3(2,:),edges3,'Normalization','cdf');
cdfb3 = histcounts(bootci3(1,:),edges3,'Normalization','cdf');
cdfu4 = histcounts(bootci4(2,:),edges4,'Normalization','cdf');
cdfb4 = histcounts(bootci4(1,:),edges4,'Normalization','cdf');
cdfu5 = histcounts(bootci5(2,:),edges5,'Normalization','cdf');
cdfb5 = histcounts(bootci5(1,:),edges5,'Normalization','cdf');

X3 = (edges3(1:end-1)+edges3(2:end))/2;
X4 = (edges4(1:end-1)+edges4(2:end))/2;
X5 = (edges5(1:end-1)+edges5(2:end))/2;

subplot(1,3,1);
fill([X3,fliplr(X3)],[cdfb3,fliplr(cdfu3)],[0.8,0.8,0.8],'edgecolor','none','FaceAlpha',0.5);
hold on;
cdfplot(initRateSum3);
cdfplot(meanShuffledData3);
subplot(1,3,2);
fill([X4,fliplr(X4)],[cdfb4,fliplr(cdfu4)],[0.8,0.8,0.8],'edgecolor','none','FaceAlpha',0.5);
hold on;
cdfplot(initRateSum4);
cdfplot(meanShuffledData4);
subplot(1,3,3);
fill([X5,fliplr(X5)],[cdfb5,fliplr(cdfu5)],[0.8,0.8,0.8],'edgecolor','none','FaceAlpha',0.5);
hold on;
cdfplot(initRateSum5);
cdfplot(meanShuffledData5);
legend

%%
%idces from hebb_pattern_showcase.m
if ~exist('rings','var')
    load over1HzRings.mat
    load rings.mat
    load 114_sorted_file_path.mat
    %freg=load('sel_conn_chain_duo_6s_1_2.mat','pair_reg','pair_chain');
    %load reg_keep.mat
    cd('/home/xd/pkg')
    addpath(fullfile('npy-matlab-master','npy-matlab'));
    addpath('fieldtrip-20200521');
    ft_defaults
    delay=6;
    if isunix
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
end


over1HzRingData = cell(length(ringConnListFull),1);
for ssidx = 1:length(ringConnListFull)
    currRings = ringConnListFull{ssidx};
    if isempty(currRings)
        continue;
    end
    over1HzRingData{ssidx} = [];
    folder=regexp(sorted_fpath{ssidx},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    for ridx = 1:size(currRings,1)
        currRing = currRings(ridx,:);
        msize = length(find(currRing ~= -1));
        transIds = currRing(1:msize)';
        [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel);
        cfg.trials = find(spktrial.trialinfo(:,8)==delay);
        for tidx=cfg.trials(:)'
            ts_id=[];
            for seqid=1:msize
                tsel=spktrial.trial{seqid}==tidx;
                ts=spktrial.time{seqid}(tsel);
                ts=ts(ts>1 & ts<7);
                ts(2,:)=seqid;
                ts_id=[ts_id;ts'];
            end
            if length(ts_id) < msize + 1
                continue;
            end
            [~,s]=sort(ts_id(:,1));
            ts_id=ts_id(s,:);
            tempdatasum=relax_tag(ts_id,msize);
            if isempty(tempdatasum)
                continue;
            end
            ringNum = size(tempdatasum,1);
            tempInfoMat = -1*ones(ringNum,11);
            tempInfoMat(:,1:4) = tempdatasum;
            tempInfoMat(:,5) = tidx;
            tempInfoMat(:,6) = msize;
            tempInfoMat(:,7:6+msize) = ones(ringNum,1)*transIds(:)';
            over1HzRingData{ssidx} = [over1HzRingData{ssidx};tempInfoMat];
        end
    end
end
save('over1HzRingData.mat','over1HzRingData','-v7.3');

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


function [avail,out]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model)
sps=30000;
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
if isempty(trials)
    avail=false;
    out=[];
    return
end

%     trials=double(trials);
%     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];
if strcmp(model, 'full')
    cluster_ids=[sustIds;transIds;nonselIds];
elseif startsWith(model,'selec')
    cluster_ids=[sustIds;transIds];
elseif startsWith(model,'nonsel')
    cluster_ids=nonselIds;
else
    keyboard
end

%  single-unit candidate

if folderType==1
    spkTS=readNPY(fullfile(spkFolder,'spike_times.npy'));
    spkId=readNPY(fullfile(spkFolder,'spike_clusters.npy'));
elseif folderType==2
    fstr=load(fullfile(spkFolder,'spike_info.mat'));
    spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
    spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);
end
FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
end
%  continuous format F T struct file
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end



function out=clearBadPerf(facSeq, mode)
if strcmp(mode, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
else
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0);
            facSeq(i-39:i,10)=good;
            if nnz(good)>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(all(facSeq(:,9:10),2),:);
    else
        out=[];
    end
end
end


function out=relax_tag(in,msize)
out = [];
skiptag=0;
for i=1:length(in)
    if i<skiptag
        continue
    end
    curr_su=in(i,2);
    targets=[(curr_su+1):(curr_su+msize-1),curr_su];
    targets(targets>msize)=targets(targets>msize)-msize;
    tsseq=[i,in(i,1:2)];
    for t=targets
        rows=tsseq(end,1)+(1:msize*10);
        rows(rows>length(in))=[];
        if isempty(rows)
            break
        end
        didx=find( ...
            in(rows,2)==t ... %post unit
            & in(rows,1)<tsseq(end,2)+0.01 ...
            & in(rows,1)>tsseq(end,2)+0.0005 ... %matching time window, assuming 1kHz
            ,1);
        if isempty(didx)
            break
        else
            tsseq=[tsseq;tsseq(end,1)+didx,in(tsseq(end,1)+didx,1:2)];
        end
    end
    if length(tsseq)<msize+1
        continue
    else
        out = [out;in(i,2),in(i,1),tsseq(end,3),tsseq(end,2)];
        skiptag=tsseq(2,1);
    end
end
end