% ring list(incongruent)
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
%bin=-2;
%fbase=load(sprintf('0906_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));


% load('reg_keep.mat','reg_set')
%if (exist('delay_data','var') && delay_data) || (exist('delay_inact_data','var') && delay_inact_data)
rings_all=cell(3,114,6,2);
rings_inact_all=cell(3,114,6,2);
delay_data = true;
delay_inact_data = true;
%rings_inact=cell(3,114,6,2);
for I=[23,46,92]
    for midx=1:3
        msize=midx+2;
        disp(I)
        lbound=100000*I;
        ubound=100000*(I+1);
        for bin=1:6
            sel11=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
            if nnz(sel11)>0
                if  delay_data
                    onering1=count_motif_all(fstr{bin}.conn_chain_S1(sel11,:),fstr{bin}.reg_chain_S1(sel11,:),fstr{bin}.pref_chain_S1(sel11,:),bin,msize);
                    onering1=unique(flexsort(onering1),'rows');
                    sel12=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound & diff(fstr{bin}.reg_chain_S2,1,2);
                    onering2=count_motif_all(fstr{bin}.conn_chain_S2(sel12,:),fstr{bin}.reg_chain_S2(sel12,:),fstr{bin}.pref_chain_S2(sel12,:),bin,msize);
                    onering2=unique(flexsort(onering2),'rows');
                    rings_all(midx,I,bin,:)={onering1,onering2};
                end
                if  delay_inact_data
                    onering1=count_motif_all_inact(fstr{bin}.conn_chain_S1(sel11,:),fstr{bin}.reg_chain_S1(sel11,:),fstr{bin}.pref_chain_S1(sel11,:),msize);
                    onering1=unique(flexsort(onering1),'rows');
                    sel12=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound & diff(fstr{bin}.reg_chain_S2,1,2);
                    onering2=count_motif_all_inact(fstr{bin}.conn_chain_S2(sel12,:),fstr{bin}.reg_chain_S2(sel12,:),fstr{bin}.pref_chain_S2(sel12,:),msize);
                    onering2=unique(flexsort(onering2),'rows');
                    rings_inact_all(midx,I,bin,:)={onering1,onering2};
                end
            end
        end

    end
end


rings_all_fast = rings_all(:,[23,46,92]',:,:);
rings_inact_all_fast = rings_inact_all(:,[23,46,92]',:,:);



% load congruent rings
load rings.mat

rings_congru_fast = rings(:,[23,46,92]',:,:);
rings_inact_congru_fast = rings_inact(:,[23,46,92]',:,:);
rings_incongru = cell(3,3,6,2);
rings_incongru_inact = cell(3,3,6,2);
for midx = 1:3
    for ssidx = 1:3
        for bin = 1:6
            for sample = 1:2
                if isempty(rings_all_fast{midx,ssidx,bin,sample})
                    rings_incongru{midx,ssidx,bin,sample} = [];
                elseif isempty(rings_congru_fast{midx,ssidx,bin,sample})
                    rings_incongru{midx,ssidx,bin,sample} = rings_all_fast{midx,ssidx,bin,sample};
                else
                    rings_incongru{midx,ssidx,bin,sample} = setdiff(rings_all_fast{midx,ssidx,bin,sample},rings_congru_fast{midx,ssidx,bin,sample},'rows');
                end
                if isempty(rings_inact_all_fast{midx,ssidx,bin,sample})
                    rings_incongru_inact{midx,ssidx,bin,sample} = [];
                elseif isempty(rings_inact_congru_fast{midx,ssidx,bin,sample})
                     rings_incongru_inact{midx,ssidx,bin,sample} = rings_inact_all_fast{midx,ssidx,bin,sample};
                else
                     rings_incongru_inact{midx,ssidx,bin,sample} = setdiff(rings_inact_all_fast{midx,ssidx,bin,sample},rings_inact_congru_fast{midx,ssidx,bin,sample},'rows');
                end
            end
        end
    end
end



%rings_congru_remap = cellfun(@(x) unique(x,'rows'),rings_congru_remap,'UniformOutput',false);
%rings_congru_inact_remap = cellfun(@(x) unique(x,'rows'),rings_congru_inact_remap,'UniformOutput',false);
%rings_incongru = cellfun(@(x,y) setdiff(x,y,'rows'),rings_all_fast,rings_congru_fast,'UniformOutput',false);
%rings_incongru_inact = cellfun(@(x,y) setdiff(x,y,'rows'),rings_inact_all_fast,rings_inact_congru_fast,'UniformOutput',false);

%rings_missed = cellfun(@(x,y) setdiff(x,y,'rows'),rings_congru_remap,rings_all_remap,'UniformOutput',false); % should be empty

save('incongru_rings_bs.mat','rings_incongru','rings_incongru_inact','-v7.3');






function out=count_motif_all(in,reg,pref,bin,msize)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)>0 & pref(:,bin+6)>0 & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end
end

function out=count_motif_all_inact(in,reg,pref,msize)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0 & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end
end










function out=count_motif_congru(in,reg,pref,bin,msize,sample)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end
end

function out=count_motif_congru_inact(in,reg,pref,msize,sample)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==sample & max(pref(:,7:12),[],2)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end
end



function out=flexsort(in)
out=in;
for i=1:size(out,1)
    [~,I]=min(out(i,:));
    while I>1
        out(i,:)=circshift(out(i,:),1);
        [~,I]=min(out(i,:));
    end
end
end


function [outconn,outreg,outpref]=shuffle_conn_chain(in,inpair,inpair_reg,inpair_pref)
outconn=nan(size(in));
outreg=nan(size(in));
outpref=nan(size(in,1),12);
combined=false(size(in,1),1);
for i=1:114
    lbound=i*100000;
    ubound=(i+1)*100000;
    sel=in(:,1)>=lbound & in(:,1)<ubound;
    combined=combined | sel;
    if nnz(sel)>0
        selpair=find(inpair(:,1)>=lbound & inpair(:,1)<ubound);
        shufsel=randperm(nnz(selpair));
        %        keyboard
        shufdata=inpair(selpair(shufsel(1:nnz(sel))),:);
        flipsel=randi(2,size(shufdata,1),1)>1;
        shufdata(flipsel,:)=shufdata(flipsel,[2,1]);
        outconn(sel,:)=shufdata;
        if exist('inpair_reg','var')
            outreg(sel,:)=inpair_reg(selpair(shufsel(1:nnz(sel))),:);
        end
        
        if exist('inpair_pref','var')
            outpref(sel,:)=inpair_pref(selpair(shufsel(1:nnz(sel))),:);
        end
    end
end
% sum(combined)
end





function sesscnt=countConnedSession()
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
sesscnt=zeros(114,1);

if (exist('delay_data','var') && delay_data) || (exist('delay_inact_data','var') && delay_inact_data)
    for I=1:114
        for midx=1:3
            msize=midx+2;
            disp(I)
            lbound=100000*I;
            ubound=100000*(I+1);
            for bin=1:6
                sel11=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
%                 keyboard()
                if nnz(sel11)>0
                    sesscnt(I)=1;
                end
            end
            
        end
    end
end
end