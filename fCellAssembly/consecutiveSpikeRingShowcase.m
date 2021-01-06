% peek inter-connected rings

% consecutive spike rings showcase

%idces from hebb_pattern_showcase.m
if ~exist('rings','var')
    load rings.mat
    load 114_sorted_file_path.mat
    load ringSeqListNew.mat % this version for screened data
    % reg
    cid = h5read('transient_6.hdf5','/cluster_id');
    regs = h5read('transient_6.hdf5','/reg');
    paths = h5read('transient_6.hdf5','/path');
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

for sidx = 1:length(ringSeqList)
    if isempty(ringSeqList{sidx})
        continue;
    end
    folder=regexp(sorted_fpath{sidx},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    datasum = ringSeqList{sidx};
    trials = unique(datasum(:,5));

    for tidx = 1:length(trials)
         tid = trials(tidx);
         currTrialData = datasum(datasum(:,5) == tid,:);
         currSuList = currTrialData(:,7:end);
         % list unique rings
        currSuList(isnan(currSuList)) = -1;
        uniqueRings = unique(currSuList,'rows');
        for ridx = 1:size(uniqueRings,1)
            transIds = uniqueRings(ridx,uniqueRings(ridx,:) ~= -1)';
            %if sum(ismember(currSuList,uniqueRings(ridx,:),'rows')) > 5 
            if sum(ismember(currSuList,uniqueRings(ridx,:),'rows'))>18 && numel(transIds)==3
                % if meet the condition, go on and plot
                %transIds = uniqueRings(ridx,uniqueRings(ridx,:) ~= -1)';
                % localize
                msize = numel(transIds);
                dset = cell(msize,2);
                probe2path = strrep(sorted_fpath{sidx},'imec0','imec1');
                for suid = 1:msize
                    try
                        if transIds(suid) < 10000
                            TF = contains(string(paths),sorted_fpath{sidx});
                            currCids = cid(TF);
                            currRegs = regs(TF);
                            reg = currRegs(currCids == transIds(suid));
                        else
                            TF = contains(string(paths),probe2path);
                            currCids = cid(TF);
                            currRegs = regs(TF);
                            reg = currRegs(currCids == transIds(suid)-10000);
                        end
                        dset{suid,1} = transIds(suid);
                        dset(suid,2) = reg;
                    catch
                        disp('check!')
                    end
                end
                [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel,tid);
                ts_id=[];
                spkrow= find(ismember(currSuList,uniqueRings(ridx,:),'rows'));
                for seqid = 1:msize
                    tsel=spktrial.trial{seqid}==1;
                    ts=spktrial.time{seqid}(tsel);
                    ts=ts(ts>1 & ts<7);
                    ts(2,:)=seqid;
                    ts(3,:) = 0;
                    for k = 1:length(spkrow)
                        ts(3,ts(1,:)> currTrialData(spkrow(k),2)-0.0005 & ts(1,:)< currTrialData(spkrow(k),4)+0.0005) = 1;
                    end
                    ts_id=[ts_id;ts'];
                end
                %[~,s]=sort(ts_id(:,1));
                %t=ts_id(s,:) 
                %pause;       
                % plot
                fh=figure('Color','w','Position',[100,100,2400,400]);
                hold on
                color={'r','b','c','m','g'};
                ph=matlab.graphics.chart.primitive.Line.empty(0,5);
                for i=1:msize
                    selp=ts_id(:,2)==i & ts_id(:,3)==1;
                    seln=ts_id(:,2)==i & ts_id(:,3)==0;
                    plot(repmat(ts_id(seln,1)',2,1),repmat([i-0.6;i+0.6],1,nnz(seln)),'-','Color',[0.8,0.8,0.8]);
                    hold on;
                    plot(repmat(ts_id(selp,1)',2,1),repmat([i-0.6;i+0.6],1,nnz(selp)),'-','Color',color{i});
                end
                xlabel('Time (s)')
                ylabel('Unique SU')
                set(gca,'XTick',1:0.5:7)
                if msize==3
                    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d',sidx,tid,transIds(1),transIds(2),transIds(3)));
                    saveas(gcf,sprintf('ShowcaseSess%dTrial%dRing%d.fig',sidx,tid,ridx));
                    save(sprintf('RegSess%dTrial%dRing%d.mat',sidx,tid,ridx),'dset');
                elseif msize==4
                    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d, %d',sidx,tid,transIds(1),transIds(2),transIds(3),transIds(4)));
                    saveas(gcf,sprintf('ShowcaseSess%dTrial%dRing%d.fig',sidx,tid,ridx));
                    %save(sprintf('RegSess%dTrial%dRing%d.mat',sidx,tid,ridx),'dset');
                elseif msize==5
                    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d, %d, %d',sidx,tid,transIds(1),transIds(2),transIds(3),transIds(4),transIds(5)));
                    saveas(gcf,sprintf('ShowcaseSess%dTrial%dRing%d.fig',sidx,tid,ridx));
                    %save(sprintf('RegSess%dTrial%dRing%d.mat',sidx,tid,ridx),'dset');
                end
                %saveas(gcf,sprintf('ShowcaseSess%dTrial%dRing%d.fig',sidx,tid,ridx));
                close(gcf);
                %save(sprintf('RegSess%dTrial%dRing%d.mat',sidx,tid,ridx),'reg');
            end
        end
    end
end
            
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


function [avail,out]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model,trialInfo)
sps=30000;
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
trials = trials(trialInfo,:);
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

