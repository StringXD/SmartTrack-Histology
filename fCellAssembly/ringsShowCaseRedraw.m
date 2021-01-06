% redraw ring connection showcase
load rings.mat
load s1s2TrialStatLen3.mat
%load ringSeqList23.mat
load ringsDataSum.mat

% find s1 and s2 trial id
midx = 1;
r1=[];
for bb=1:6
    r1=[r1;cell2mat(rings(midx,:,bb,1)')];
end
r2=[];
for bb=1:6
    r2=[r2;cell2mat(rings(midx,:,bb,2)')];
end
r1(:,midx+3)=1;r2(:,midx+3)=2;ringsm=[r1(:,1:(midx+2));r2(:,1:(midx+2))];samps=[r1(:,midx+3);r2(:,midx+3)];
metadata = zeros(length(ringsm),midx+3);
for ssidx=1:length(ringsm)
    sessIdx=idivide(int32(ringsm(ssidx,1)),int32(100000));
    transIds=int32(rem(ringsm(ssidx,:),100000))';
    metadata(ssidx,1) = sessIdx;
    metadata(ssidx,2:end) = transIds';
end
trials32 = s1s2Trials(metadata(:,1)==23,:);
s1trials = trials32{1,1};
s2trials = trials32{1,2};

%trials = unique(ringSeqList23(:,5));
trials = unique(ringsDset{23}(:,5));
%colors = [0.93,0.69,0.13;0.47,0.67,0.19;0,0.45,0.74;1,0.07,0.65];
% same color
colors = [0.75,0,0;0.75,0,0;0.75,0,0;0.75,0,0];
cidx = 0;
%% part1. Plot
for tidx = 17
    subplot(2,1,1);
    tid = trials(tidx);
    %currTrialData = ringSeqList23(ringSeqList23(:,5) == tid,:);
    % screen for > 1Hz rings
    ringSS23Data = ringsDset{23};
    uniRingsAll = unique(ringSS23Data(:,7:end),'rows');
    ringFreq = zeros(size(uniRingsAll,1),1);
    for ridx = 1:size(uniRingsAll,1)
        thisRing = uniRingsAll(ridx,:);
        if sum(ismember(ringSS23Data(:,7:end),thisRing,'rows')) >= length(trials)
            ringFreq(ridx) = 1;
        end
    end
    currTrialData = ringsDset{23}(ismember(ringSS23Data(:,7:end),uniRingsAll(logical(ringFreq),:),'rows') & ringSS23Data(:,5) == tid,:);
    currSuList = currTrialData(:,7:end);
    rseData = [currTrialData,zeros(size(currTrialData,1),1)];
    for lineIdx = 1:size(currTrialData,1)
        if rseData(lineIdx,end) ~= 0
            continue;
        end
        rows = find(rseData(:,2) == rseData(lineIdx,4) |  rseData(:,4) == rseData(lineIdx,2));
        if ~isempty(rows)
            for crow = rows(:)'
                if rseData(crow,6+rseData(crow,1)) == rseData(lineIdx,6+rseData(lineIdx,3)) ||  rseData(crow,6+rseData(crow,3)) == rseData(lineIdx,6+rseData(lineIdx,1))
                    if rseData(crow,end) == 0
                        cidx = cidx + 1;
                        if cidx > 4 
                            cidx = cidx - 4;
                        end
                        rseData(crow,end) = cidx; 
                    end
                    rseData(lineIdx,end) = rseData(crow,end);
                end
            end
        end
    end
    uniqueRings = unique(currSuList,'rows');
    
    
    
    
    colors = [0.5,0.5,0.5;colors];
    
    rseData = sortrows(rseData,7:11);
    plotrowId = 1;
    for lineIdx = 1:size(rseData,1)
        if lineIdx ~= 1 && ~all(rseData(lineIdx,7:11) == rseData(lineIdx-1,7:11))
            plotrowId = plotrowId + 1;
        end
        line([rseData(lineIdx,2),rseData(lineIdx,4)],[plotrowId,plotrowId],'Color',colors(rseData(lineIdx,end)+1,:),'LineWidth',2);
        hold on;
    end
%     % 
%     rseSum = []; % store already connected rings
%     connPairId = [];
%     for ridx = 1:size(uniqueRings,1)
%         currRingData = currTrialData(ismember(currSuList,uniqueRings(ridx,:),'rows'),:);
%         if isempty(currRingData)
%             continue;
%         end
%         for rpt = 1:size(currRingData,1)
%             if ~isempty(rseSum)
%                 if ~isempty(find(rseSum(:,2) == currRingData(rpt,4),1))
%                     [row,~] = find(rseSum(:,7:11) == currRingData(rpt,6+currRingData(rpt,3)));
%                     for rowidx = 1:length(row)
%                         if rseSum(row(rowidx),2) == currRingData(rpt,4)
%                             if sum(rseSum(row(rowidx),end-2:end)==[0.5,0.5,0.5]) == 0
%                                 cidx = cidx + 1;
%                                 if cidx > 4
%                                     cidx = cidx - 4;
%                                 end
%                                 rseSum(row(rowidx),end-2:end) = colors(cidx,:);
%                             end
%                         end
%                         if rowidx == 1
%                             rseSum = [rseSum;[currRingData(rpt,:),rseSum(row(rowidx),end-2:end)]];
%                         end
%                     end
%                     %h2=line([currRingData(rpt,4),currRingData(rpt,4)],[rseSum(find(rseSum(:,2) == currRingData(rpt,4),1),1),ridx],'Color','red','LineStyle','--');
%                     %connPairId = [connPairId;[rseSum(find(rseSum(:,2) == currRingData(rpt,4),1),1),ridx]];
%                     %xline(currRingData(rpt,4),'--r');
%                     %hold on;
%                 elseif ~isempty(find(rseSum(:,4) == currRingData(rpt,2),1))
%                     [row,~] = find(rseSum(:,7:11) == currRingData(rpt,6+currRingData(rpt,1)));
%                     for rowidx = 1:length(row)
%                         if rseSum(row(rowidx),4) == currRingData(rpt,2)
%                             if sum(rseSum(row(rowidx),end-2:end)==[0.5,0.5,0.5]) == 0
%                                 cidx = cidx + 1;
%                                 if cidx > 4
%                                     cidx = cidx - 4;
%                                 end
%                                 rseSum(row(rowidx),end-2:end) = colors(cidx,:);
%                             end
%                         end
%                         if rowidx == 1
%                             rseSum = [rseSum;[currRingData(rpt,:),rseSum(row(rowidx),end-2:end)]];
%                         end
%                     end
%                     %h2=line([currRingData(rpt,2),currRingData(rpt,2)],[rseSum(find(rseSum(:,3) == currRingData(rpt,2),1),1),ridx],'Color','red','LineStyle','--');
%                     %connPairId = [connPairId;[rseSum(find(rseSum(:,3) == currRingData(rpt,2),1),1),ridx]];
%                     %xline(currRingData(rpt,2),'--r');
%                     %hold on;
%                 else
%                     rseSum = [rseSum;[currRingData,ones(size(currRingData,1),1)*[0.5,0.5,0.5]]]; 
%                     %h1=line([currRingData(rpt,2),currRingData(rpt,4)],[ridx,ridx],'Color',[0.5,0.5,0.5],'LineWidth',2);
%                     %hold on;
%                 end
%             else
%                rseSum = [rseSum;[currRingData,ones(size(currRingData,1),1)*[0.5,0.5,0.5]]]; 
%             end
%             %rseSum = [rseSum;ridx,currRingData(rpt,2),currRingData(rpt,4)];
%         end
%     end
    %plot
    %legend([h1,h2],{'Ring Activity','Consecutive Ring-ring Conn.'});
    ax = gca;
    ax.XLim = [1 7];
    xlabel('Time (s)');
    ylabel('Unique Rings')
    title(sprintf('Sess#%d, Trial#%d',23,tid));
    subplot(2,1,2);
    edges = 1:0.05:7;
    counts = zeros(1,length(edges)-1);
    pcounts = zeros(1,length(edges)-1);
    midpoint = zeros(1,length(edges)-1);
    for i = 1:length(counts)
        counts(i) = sum(~(rseData(:,2)> edges(i+1) | rseData(:,4)<edges(i)));
        pcounts(i) = sum(~(rseData(:,2)> edges(i+1) | rseData(:,4)<edges(i)) & rseData(:,end) ~= 0);
        midpoint(i) = mean(edges(i:i+1));
    end
    % frequency normalized unique ring number
    freq = counts/(0.05*plotrowId);
    pfreq = pcounts/(0.05*plotrowId);
    b = bar(midpoint,[pfreq;freq-pfreq]','stacked');
    b.FaceColor = [0.8,0.8,0.8];
    b.CData(:,1) = [0.75,0,0];
    b.EdgeColor = 'none';
    ylabel('Frequency(Hz)');
    xlabel('Time (s)');        
end

%% Part2. sort connected rings closer
cRings = unique(connPairId);
treeList = [];
candidateList = [];
connCounts = zeros(length(cRings),1);
for i = 1:length(cRings)
    connCounts(i) = length(find(connPairId==cRings(i)));
end
while ~isempty(cRings)
    if isempty(treeList)
        winnerId = find(connCounts == max(connCounts),1);
        treeList = [treeList;cRings(winnerId)];
        connCounts(cRings == treeList) = [];
        cRings(cRings == treeList) = [];
        candidateList = cRings;
    else
        scores = zeros(length(candidateList),1);
        for cidx = 1:length(candidateList)
            [~,ia,~] = intersect(connPairId(:,1),treeList);
            [~,ib,~] = intersect(connPairId(:,2),treeList);
            [ic,~] = find(connPairId == candidateList(cidx));
            scores(cidx) = length(intersect(union(ia,ib),ic));
        end
        if length(find(scores==max(scores))) > 1
            temp = candidateList(scores==max(scores));
            tempScore = connCounts(scores==max(scores));
            treeList = [treeList;temp(find(tempScore==max(tempScore),1))];
        else
            treeList = [treeList;candidateList(scores==max(scores))];
        end
        connCounts(candidateList == treeList(end)) = [];
        cRings(cRings == treeList(end)) = [];
        candidateList(candidateList == treeList(end)) = [];
    end
end
ringListNew = (1:size(uniqueRings,1))';
[~,idd,~] = intersect(ringListNew,treeList);
ringListNew(idd) = [];
ringListNew = [treeList;ringListNew];
%% Part3.
for tidx = 17
    subplot(2,1,1);
    tid = trials(tidx);
    currTrialData = ringSeqList23(ringSeqList23(:,5) == tid,:);
    currSuList = currTrialData(:,7:end);
    uniqueRings = unique(currSuList,'rows');
    rseSum = [];
    connPairId = [];
    ii = 1:size(uniqueRings,1);
    for ridx = ringListNew'
        currRingData = currTrialData(ismember(currSuList,uniqueRings(ridx,:),'rows'),:);
        if isempty(currRingData)
            continue;
        end
        for rpt = 1:size(currRingData,1)
            h1=line([currRingData(rpt,2),currRingData(rpt,4)],[ridx,ridx],'LineWidth',2);
            hold on;
            if ~isempty(rseSum)
                if ~isempty(find(rseSum(:,2) == currRingData(rpt,4),1))
                    h2=line([currRingData(rpt,4),currRingData(rpt,4)],[ii(ringListNew==rseSum(find(rseSum(:,2) == currRingData(rpt,4),1),1)),ii(ringListNew==ridx)],'Color','red','LineStyle','--');
                    connPairId = [connPairId;[rseSum(find(rseSum(:,2) == currRingData(rpt,4),1),1),ridx]];
                    %xline(currRingData(rpt,4),'--r');
                    hold on;
                end
                if ~isempty(find(rseSum(:,3) == currRingData(rpt,2),1))
                    h2=line([currRingData(rpt,2),currRingData(rpt,2)],[ii(ringListNew==rseSum(find(rseSum(:,3) == currRingData(rpt,2),1),1)),ii(ringListNew==ridx)],'Color','red','LineStyle','--');
                    connPairId = [connPairId;[rseSum(find(rseSum(:,3) == currRingData(rpt,2),1),1),ridx]];
                    %xline(currRingData(rpt,2),'--r');
                    hold on;
                end
            end
            rseSum = [rseSum;ridx,currRingData(rpt,2),currRingData(rpt,4)];
        end
    end
    legend([h1,h2],{'Ring Activity','Consecutive Ring-ring Conn.'});
    ax = gca;
    ax.XLim = [1 7];
    xlabel('Time (s)');
    ylabel('Unique Rings')
    title(sprintf('Sess#%d, Trial#%d',23,tid));
    subplot(2,1,2);
    edges = 1:0.05:7;
    counts = zeros(1,length(edges)-1);
    midpoint = zeros(1,length(edges)-1);
    for i = 1:length(counts)
        counts(i) = sum(~(rseSum(:,2)> edges(i+1) | rseSum(:,3)<edges(i)));
        midpoint(i) = mean(edges(i:i+1));
    end
    bar(midpoint,counts,'FaceColor',[0.8,0.8,0.8],'EdgeColor','none');
    ylabel('Counts');
    xlabel('Time (s)');        
end