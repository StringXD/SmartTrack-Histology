% Are there isolated rings?
load ringsDataSum.mat
ringIDList = cell(length(ringsDset),2);
for ssidx = 23
    currData = ringsDset{ssidx};
    if isempty(currData)
        continue;
    end
    ringIDList{ssidx,1} = unique(currData(:,7:11),'rows');
    tempRings = ringIDList{ssidx,1};
    if size(tempRings,1) < 10
        continue;
    end
    for ridx = 1:size(tempRings,1)
        shareSuRingStat = xor(any(tempRings == tempRings(ridx,:) & tempRings ~= -1,2),all(tempRings == tempRings(ridx,:),2));
        if ridx == 1
            ringIDList{ssidx,2} = shareSuRingStat;
        else
            ringIDList{ssidx,2} = ringIDList{ssidx,2} | shareSuRingStat;
        end
    end
end

%% showcase
% > 1Hz rings meaning > 106
for ssidx = 23
    ringSS23Data = ringsDset{23};
    uniRingsAll = unique(ringSS23Data(:,7:end),'rows');
    ringFreq = zeros(size(uniRingsAll,1),1);
    for ridx = 1:size(uniRingsAll,1)
        thisRing = uniRingsAll(ridx,:);
        if sum(ismember(ringSS23Data(:,7:end),thisRing,'rows')) >= 106
            ringFreq(ridx) = 1;
        end
    end
    tid = 17;
    currTrialData = ringsDset{23}(ismember(ringSS23Data(:,7:end),uniRingsAll(logical(ringFreq),:),'rows') & ringSS23Data(:,5) == tid,:);
    %uniRings = unique(currTrialData(:,7:end),'rows');
    %uniStat = zeros(size(uniRings,1),1);
    % share in entire delay period
%     for ridx = 1:size(uniRings,1)
%         transIds = uniRings(ridx,uniRings(ridx,:)~=-1);
%         for nidx = 1:length(transIds)
%             if nidx == 1
%                 shareSuRingStat = any(uniRings == transIds(nidx),2);
%             else
%                 shareSuRingStat = shareSuRingStat | any(uniRings == transIds(nidx),2);
%             end
%         end    
%         shareSuRingCount = sum(shareSuRingStat)-1;
%         if shareSuRingCount > 0
%             uniStat(ridx) = 1;
%         end
%     end
    shareSuRingDataBinned = zeros(6,2);
    midTimes = (currTrialData(:,2) + currTrialData(:,4))/2;
    for bin = 1:6
        currBinData = currTrialData(midTimes > bin & midTimes < bin+1,:);
        % share in each 1s bin
        uniRingBinned = unique(currBinData(:,7:11),'rows');
        urbStat = zeros(size(uniRingBinned,1),1);
        for tempridx = 1:size(uniRingBinned,1)
            transIds = uniRingBinned(tempridx,uniRingBinned(tempridx,:)~=-1);
            for nidx = 1:length(transIds)
                if nidx == 1
                    shareSuRingBinnedStat = any(uniRingBinned == transIds(nidx),2);
                else
                    shareSuRingBinnedStat = shareSuRingBinnedStat | any(uniRingBinned == transIds(nidx),2);
                end
            end
            shareSuRingBinnedCount = sum(shareSuRingBinnedStat) - 1;
            if shareSuRingBinnedCount > 0
                urbStat(tempridx) = 1;
            end
        end
        shareSuRingDataBinned(bin,1) = sum(ismember(currBinData(:,7:end),uniRingBinned(logical(urbStat),:),'rows'));
        shareSuRingDataBinned(bin,2) = size(currBinData,1);
    end
end
scatter(1:6,shareSuRingDataBinned(:,1)./shareSuRingDataBinned(:,2),'filled');
hold on;
plot(1:6,shareSuRingDataBinned(:,1)./shareSuRingDataBinned(:,2));
set(gca,'YLim',[0 1.1]);     
xlabel('Time bins(s)');
ylabel('Fraction of inter-connected rings');

%% stepped binsize
% Compute the overall fraction for each bin and average over time

load ringsDataSum.mat
load tnumList.mat
allSSData = cell(114,1);
for ssidx = 1:114
    ringSSData = ringsDset{ssidx};
    if size(ringSSData,1) < 1000
        continue;
    end
    uniRingsAll = unique(ringSSData(:,7:end),'rows');
    ringFreq = zeros(size(uniRingsAll,1),1);
    for ridx = 1:size(uniRingsAll,1)
        thisRing = uniRingsAll(ridx,:);
        if sum(ismember(ringSSData(:,7:end),thisRing,'rows')) >= tnumList{ssidx}
            ringFreq(ridx) = 1;
        end
    end
    tidList = unique(ringSSData(:,5));
    binData = cell(100,length(tidList));
    for tidx = 1:length(tidList)
        tid = tidList(tidx);
        currTrialData = ringSSData(ismember(ringSSData(:,7:end),uniRingsAll(logical(ringFreq),:),'rows') & ringSSData(:,5) == tid,:);
        binSize = 10:10:1000;
        for binSizeIdx = 1:100
            binEdges = 1:binSize(binSizeIdx)/1000:7;
            shareSuRingDataBinned = zeros(length(binEdges)-1,2);
            midTimes = (currTrialData(:,2) + currTrialData(:,4))/2;
            for bin = 1:length(binEdges)-1
                currBinData = currTrialData(midTimes > binEdges(bin) & midTimes < binEdges(bin+1),:);
                uniRingBinned = unique(currBinData(:,7:11),'rows');
                urbStat = zeros(size(uniRingBinned,1),1);
                for tempridx = 1:size(uniRingBinned,1)
                    transIds = uniRingBinned(tempridx,uniRingBinned(tempridx,:)~=-1);
                    for nidx = 1:length(transIds)
                        if nidx == 1
                            shareSuRingBinnedStat = any(uniRingBinned == transIds(nidx),2);
                        else
                            shareSuRingBinnedStat = shareSuRingBinnedStat | any(uniRingBinned == transIds(nidx),2);
                        end
                    end
                    shareSuRingBinnedCount = sum(shareSuRingBinnedStat) - 1;
                    if shareSuRingBinnedCount > 0
                        urbStat(tempridx) = 1;
                    end
                end
                shareSuRingDataBinned(bin,1) = sum(ismember(currBinData(:,7:end),uniRingBinned(logical(urbStat),:),'rows'));
                shareSuRingDataBinned(bin,2) = size(currBinData,1);
            end
            binData{binSizeIdx,tidx} = shareSuRingDataBinned;
        end
    end
    allSSData{ssidx} = binData;
%     ratio = zeros(100,size(binData,2));
%     for ii = 1:100
%         for tidx = 1:size(binData,2)
%             ratio(ii,tidx) = nanmean(binData{ii,tidx}(:,1)./binData{ii,tidx}(:,2));
%         end
%     end
%     plot(nanmean(ratio,2));
%     pause;
%     close all;
end
% overall plot
allBinnedData = [];
for ssidx = 1:114
    currSum = allSSData{ssidx};
    if isempty(currSum)
        continue;
    end
    ratio = zeros(100,size(currSum,2));
    for ii = 1:100
        for tidx = 1:size(currSum,2)
            ratio(ii,tidx) = nanmean(currSum{ii,tidx}(:,1)./currSum{ii,tidx}(:,2));
        end
    end
    allBinnedData = [allBinnedData;nanmean(ratio,2)'];
end
% remove abnormal
allBinnedData(23,:) = [];
meanCurve = mean(allBinnedData);
stdCurve = std(allBinnedData);
X = 10:10:1000;
fill([X,fliplr(X)],[meanCurve-stdCurve/sqrt(size(allBinnedData,1)-1),fliplr(meanCurve+stdCurve/sqrt(size(allBinnedData,1)-1))],[0.8,0.8,0.8],'edgecolor','none','FaceAlpha',0.5);
hold on;
plot(X,meanCurve,'Color',[0.75,0,0],'LineWidth',1.5)
xlabel('Bin Size(ms)');
ylabel('Inter-loop Coupling Rate');






