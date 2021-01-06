% Fraction of spikes contributing to FC
%% classification
recompute = 0;
if recompute == 1
    fstr=cell(1,6);
    for bin = 1:6
        fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    end
    connSum = cell(114,2);
    for I = 1:114
        currSessConnS1 = [];
        currSessConnS2 = [];
        for bin = 1:6
            lbound=100000*I;
            ubound=100000*(I+1);
            sel11=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
            sel12=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound & diff(fstr{bin}.reg_chain_S2,1,2);
            if nnz(sel11)>0 
                % suPre, suPost, connBin, selPre, selPost
                currSessConnS1 = [currSessConnS1;fstr{bin}.conn_chain_S1(sel11,:)-lbound,bin*ones(sum(sel11),1),max(fstr{bin}.pref_chain_S1(sel11,1:6),[],2),max(fstr{bin}.pref_chain_S1(sel11,7:12),[],2)];
            end
            if nnz(sel12)>0
                currSessConnS2 = [currSessConnS2;fstr{bin}.conn_chain_S2(sel12,:)-lbound,bin*ones(sum(sel12),1),max(fstr{bin}.pref_chain_S2(sel12,1:6),[],2),max(fstr{bin}.pref_chain_S2(sel12,7:12),[],2)];
            end
        end
        connSum{I,1} = currSessConnS1;
        connSum{I,2} = currSessConnS2;
    end
    save('connSum.mat','connSum','-v7.3');
end
%% prepare spike train data
if ~exist('connSum','var')
    load('connSum.mat');
end
load 114_sorted_file_path.mat
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

sumSum = cell(length(connSum),2); % sample1 sample2
for ssidx = 1:length(connSum)
    disp(ssidx)
    if isempty(connSum{ssidx,1})
        continue;
    end
    folder=regexp(sorted_fpath{ssidx},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    for sample = 1:2
        currData = connSum{ssidx,sample};
        suList = unique(currData(:,1:2));
        tempDataSum = zeros(length(suList),61); 
        % cols: pre&congru(6)+post&congru(6)+pre&incongru(6)+post&incongru(6)+ignoreSel(12)+ignoreIO(12)+totalfunc(6)+total(6)+sel 
        transIds = suList(:);
        [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel);
        if sample==1
            cfg.trials = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
        else
            cfg.trials = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);
        end
        for suIdx = 1:length(suList)
            [rows,loc] = find(currData(:,1:2) == suList(suIdx));
            workBins = unique(currData(:,3));
            preCase = currData(rows(loc==1),:);
            postCase = currData(rows(loc==2),:);
            for bin = workBins(:)'
                if ~isempty(preCase)
                    candidate = [suList(suIdx),preCase(1,4)];
                    post = preCase(preCase(:,3) == bin,[2,5]);
                else
                    candidate = [suList(suIdx),postCase(1,5)];
                    post = [];
                end
                if ~isempty(postCase)
                    pre = postCase(postCase(:,3) == bin,[1,4]);
                else
                    pre = [];
                end
                % main calculation here
                [cumul,total] = spikeFeeder(suList,spktrial,cfg.trials,bin,candidate,pre,post,0.002,0.01);
                tempDataSum(suIdx,1:54) = tempDataSum(suIdx,1:54) + cumul;
            end
            tempDataSum(suIdx,55:60) = total;
            tempDataSum(suIdx,end) = candidate(2);
        end
        sumSum{ssidx,sample} = tempDataSum;
    end
end
save('spike2FCRaw.mat','sumSum','-v7.3');


                
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
                    
                    
function [cumul,total] = spikeFeeder(suList,spktrial,trials,bin,candidate,pre,post,minlag,maxlag)
cumul = zeros(1,54);
total = zeros(1,6);
if ~isempty(post)
    sel_post_congru = find(ismember(suList,post(post(:,2) == candidate(2),1)));
    sel_post_incongru = find(~ismember(suList,post(post(:,2) == candidate(2),1)));
else
    sel_post_congru = [];
    sel_post_incongru = [];
end
if ~isempty(pre)
    sel_pre_congru = find(ismember(suList,pre(pre(:,2) == candidate(2),1)));
    sel_pre_incongru = find(~ismember(suList,pre(pre(:,2) == candidate(2),1)));
else
    sel_pre_congru = [];
    sel_pre_incongru = [];
end
for tidx = trials(:)'
    tsel = spktrial.trial{suList==candidate(1)}==tidx;
    ts = spktrial.time{suList==candidate(1)}(tsel);
    for tbin = 1:6
        tsTotalBin = ts(ts>tbin & ts<tbin+1);
        total(tbin) = total(tbin) + length(tsTotalBin);
    end
    tsBin = sort(ts(ts>bin & ts<bin+1));
    if isempty(tsBin)
        continue;
    end
    % put spikes together
    % pre-congru
    if ~isempty(sel_post_congru)
        spikeFeedprec = sumSpikesXcorr(spktrial,sel_post_congru,tidx,tsBin,minlag,maxlag,1);
    else
        spikeFeedprec = zeros(length(tsBin),1);
    end
    % pre-incongru
    if ~isempty(sel_post_incongru)
        spikeFeedprei = sumSpikesXcorr(spktrial,sel_post_incongru,tidx,tsBin,minlag,maxlag,1);
    else
        spikeFeedprei = zeros(length(tsBin),1);
    end
    % post-congru
    if ~isempty(sel_pre_congru)
        spikeFeedpostc = sumSpikesXcorr(spktrial,sel_pre_congru,tidx,tsBin,minlag,maxlag,2);
    else
        spikeFeedpostc = zeros(length(tsBin),1);
    end
    % post-incongru
    if ~isempty(sel_pre_incongru)
        spikeFeedposti = sumSpikesXcorr(spktrial,sel_pre_incongru,tidx,tsBin,minlag,maxlag,2);
    else
        spikeFeedposti = zeros(length(tsBin),1);
    end
    % counting
    spikeFeedpre = logical(spikeFeedprec) | logical(spikeFeedprei);
    spikeFeedpost = logical(spikeFeedpostc) | logical(spikeFeedposti);
    spikeFeedc = logical(spikeFeedprec) | logical(spikeFeedpostc);
    spikeFeedi = logical(spikeFeedprei) | logical(spikeFeedposti);
    spikeFeeda = spikeFeedpre | spikeFeedpost;
    cumul(bin) = cumul(bin) + sum(spikeFeedprec);
    cumul(bin+6) = cumul(bin+6) + sum(spikeFeedpostc);
    cumul(bin+12) = cumul(bin+12) + sum(spikeFeedprei);
    cumul(bin+18) = cumul(bin+18) + sum(spikeFeedposti);
    cumul(bin+24) = cumul(bin+24) + sum(spikeFeedpre);
    cumul(bin+30) = cumul(bin+30) + sum(spikeFeedpost);
    cumul(bin+36) = cumul(bin+36) + sum(spikeFeedc);
    cumul(bin+42) = cumul(bin+42) + sum(spikeFeedi);
    cumul(bin+48) = cumul(bin+48) + sum(spikeFeeda);
end
end


function spikeFeed = sumSpikesXcorr(spktrial,nidices,tidx,ts,minlag,maxlag,type)
% type 1: pre type 2: post
sumts = cellfun(@(x,y) y(x==tidx),spktrial.trial(nidices),spktrial.time(nidices),'UniformOutput',false);
if type == 1
    sumts = cellfun(@(x) x(x>ts(1) & x<ts(end)+maxlag),sumts,'UniformOutput',false);
else
    sumts = cellfun(@(x) x(x>ts(1)-maxlag & x<ts(end)),sumts,'UniformOutput',false);
end
sumts = sort(cell2mat(sumts)); 
% xcorr
% compute all distances at once using a multiplication trick
if length(ts)*length(sumts) < 2*10^7 % allow matrix to grow to about 150 MB, should always work
    if type == 1
        D = log(exp(-ts(:))*exp(sumts(:)'));
    else
        D = log(exp(ts(:))*exp(-sumts(:)'));
    end
    spikeFeed = any(D < maxlag & D > minlag,2);
else
    spikeFeed = zeros(length(ts),1);
    if type == 1
        for spkIdx = 1:length(ts)
            spikeFeed(spkIdx) = any(sumts-ts(spkIdx) < maxlag & sumts-ts(spkIdx) > minlag);
        end
    else
        for spkIdx = 1:length(ts)
            spikeFeed(spkIdx) = any(-sumts+ts(spkIdx) < maxlag & -sumts+ts(spkIdx) > minlag);
        end
    end
end
end