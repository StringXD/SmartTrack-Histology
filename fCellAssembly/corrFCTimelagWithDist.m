% correlation between FC timelag and inter-neuron distance
%% construct timelag dataset
timelagWin = -9:2:9;
thresh=norminv(0.995); %0.05 bonferroni corrected
for bin = 1:6
    load(sprintf('full_duo_XCORR_sums_delay_6_%d_%d_2msbin.mat',bin,bin+1));
    [~,I] = sort(cell2mat(sums(:,1)));
    sums = sums(I,:);
    if bin == 1
        timelagData = cell(size(sums,1),6,2);
        sidList = cell2mat(sums(:,1));
    end
    for sidx = 1:size(sums,1)
        fprintf('%d of %d\n',sidx, size(sums,1));
        xc_s1=sums{sidx,5};
        xshuf_s1=sums{sidx,6};
        xc_s2=sums{sidx,7};
        xshuf_s2=sums{sidx,8};
        if isempty(xc_s1) || isempty(xc_s2)
            continue
        end
        timelagMat1 = [];
        timelagMat2 = [];
        for si=1:(size(xc_s1.xcorr,1)-1)
            su1id=str2double(xc_s1.label{si,1});
            for sj=(si+1):size(xc_s1.xcorr,2)
                su2id=str2double(xc_s1.label{sj,1});
                peaks1=NaN;
                peaks2=NaN;
                totalCount_S1=nansum(squeeze(xc_s1.xcorr(si,sj,:)));
                totalCount_S2=nansum(squeeze(xc_s2.xcorr(si,sj,:)));
                if numel(xc_s1.cfg.trials)<20 ||  numel(xc_s2.cfg.trials)<20 || (totalCount_S1<250 && totalCount_S2<250)
                    continue;
                end
                %sample 1 stats
                hists1=squeeze(xc_s1.xcorr(si,sj,:));
                shufs1=squeeze(xshuf_s1.shiftpredictor(si,sj,:));
                diffs1=hists1-smooth(shufs1);
                diffs1(50:51)=0;
                stds1=std(shufs1);
                scores1=diffs1(46:55)./stds1;
                if any(scores1>thresh)
                    [peaks1,timelagInd1] = nanmax(scores1);
                    % A peak at a negative lag (I.E. AI>0) for stat.xcorr(chan1,chan2,:) means that chan1 is leading chan2.
                    timelag1 = timelagWin(timelagInd1);
                    timelagMat1 = [timelagMat1;[su1id,su2id,timelag1]];
                end
                %sample 2 stats
                hists2=squeeze(xc_s2.xcorr(si,sj,:));
                shufs2=squeeze(xshuf_s2.shiftpredictor(si,sj,:));
                diffs2=hists2-smooth(shufs2);
                diffs2(50:51)=0;
                stds2=std(shufs2);
                scores2=diffs2(46:55)./stds2;
                if any(scores2>thresh)
                    [peaks2,timelagInd2] = nanmax(scores2);
                    % A peak at a negative lag (I.E. AI>0) for stat.xcorr(chan1,chan2,:) means that chan1 is leading chan2.
                    timelag2 = timelagWin(timelagInd2);
                    timelagMat2 = [timelagMat2;[su1id,su2id,timelag2]];
                end
            end
        end
        timelagData{sidx,bin,1} = timelagMat1;
        timelagData{sidx,bin,2} = timelagMat2;
    end
end        
save('timelagData.mat','timelagData','sidList','-v7.3');                                       
%% construct inter-nueron distance dataset
cid_list=h5read('transient_6.hdf5','/cluster_id');
path_list=h5read('transient_6.hdf5','/path');

load('sidList.mat');
load('sucoords.mat');
distDset = cell(length(sidList),2);
for i = 1:length(sidList)
    sidx = sidList(i);
    currPath = deblank(path_list{sidx});
    imecTxt = char(currPath);
    imec = str2double(imecTxt(end-8));
    if imec == 0
        alternatePath = strrep(currPath,'imec0','imec1');
        currCids = [cid_list(contains(path_list,currPath));10000+cid_list(contains(path_list,alternatePath))];
    else
        alternatePath = strrep(currPath,'imec1','imec0');
        currCids = [10000+cid_list(contains(path_list,currPath));cid_list(contains(path_list,alternatePath))];
    end
    currCoords = [coord(contains(path_list,currPath),:);coord(contains(path_list,alternatePath),:)];
    distMat = [];
    for si = 1:length(currCids)-1
        su1id = double(currCids(si));
        for sj = (si+1):length(currCids)
            su2id = double(currCids(sj));
            dist = norm(currCoords(si,:) - currCoords(sj,:));
            distMat = [distMat;[su1id,su2id,dist]];
        end
    end
    distDset{i,1} = sidx;
    distDset{i,2} = distMat;
end
save('distMat.mat','distDset','-v7.3');

%% correlation between timelag and distance and count
tList = -9:2:9;
dThres = 0:0.25:4.5;
for j = 1:6
    countMat1 = zeros(18,10);
    for i = 1:length(distDset)
        distdata = distDset{i,2};
        s1data = timelagData{i,j,1};
        if isempty(s1data) 
            continue;
        end
        [~,ia,ib] = intersect(s1data(:,1:2),distdata(:,1:2),'rows');
        s1dataRF = s1data(ia,3);
        dist1 = double(distdata(ib,3))/100;
        for tbin = 1:10
            for distbin = 1:18
                countMat1(distbin,tbin) = countMat1(distbin,tbin) + sum(s1dataRF == tList(tbin) & dist1 > dThres(distbin) & dist1 <= dThres(distbin+1));
            end
        end
    end
end
countMat1 = countMat1 / 6;
% put negative value counts to positive
countMat1(:,6:10) = countMat1(:,6:10) + countMat1(:,5:-1:1);

bar3(countMat1(end:-1:1,6:10));

set(gca,'zscale','linear')
xlabel('Time lag(ms)');
ylabel('Distance(mm)');
zlabel('Connection Number');

%% peek cross bin and sample time lag variation Full
for i = 1:105
    tempCell = cell(1,12);
    for j = 1:6
        tempCell{j} = timelagData{i,j,1};
    end
    for j = 7:12
        tempCell{j} = timelagData{i,j-6,2};
    end
    C = tempCell{7}(:,1:2);
    for selturn = 8:12
        C = union(C,tempCell{selturn}(:,1:2),'rows');
    end
    %out = zeros(length(C),12);
    out = zeros(length(C),6);
    for idx = 7:12
        [~,ia,ib] = intersect(tempCell{idx}(:,1:2),C,'rows');
        out(ib,idx-6) = tempCell{idx}(ia,3);
    end
    for tbin = 7:12
        if tbin == 7
            out(out(:,1) < 0,:) = -out(out(:,1) < 0,:);
        else
            nnzeros = zeros(length(out),1);
            for tp = 1:length(out)
                nnzeros(tp) = nnz(out(tp,1:tbin-7));
            end
            out(out(:,tbin-6) < 0 & nnzeros == 0,:) = -out(out(:,tbin-6) < 0 & nnzeros == 0,:);
        end
    end            
    
    out = sortrows(out);
    h = imagesc(out);
    set(h,'alphadata',out~=0);
    c = colorbar;
    c.Label.String = 'time lag(ms)';
    ylabel('Neuron pair');
    xlabel('Bin');
    pause;
end
