% ring list(incongruent)
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
%bin=-2;
%fbase=load(sprintf('0906_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));


% load('reg_keep.mat','reg_set')
%if (exist('delay_data','var') && delay_data) || (exist('delay_inact_data','var') && delay_inact_data)
rings_all=cell(3,114,6);
rings_inact_all=cell(3,114,6);
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
                    rings_all{midx,I,bin}=[onering1;onering2];
                end
                if  delay_inact_data
                    onering1=count_motif_all_inact(fstr{bin}.conn_chain_S1(sel11,:),fstr{bin}.reg_chain_S1(sel11,:),fstr{bin}.pref_chain_S1(sel11,:),msize);
                    onering1=unique(flexsort(onering1),'rows');
                    sel12=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound & diff(fstr{bin}.reg_chain_S2,1,2);
                    onering2=count_motif_all_inact(fstr{bin}.conn_chain_S2(sel12,:),fstr{bin}.reg_chain_S2(sel12,:),fstr{bin}.pref_chain_S2(sel12,:),msize);
                    onering2=unique(flexsort(onering2),'rows');
                    rings_inact_all{midx,I,bin}=[onering1;onering2];
                end
            end
        end

    end
end


rings_all_fast = rings_all(:,[23,46,92]',:);
rings_inact_all_fast = rings_inact_all(:,[23,46,92]',:);
rings_all_remap = cell(3,3);
rings_inact_all_remap = cell(3,3);
for midx = 1:3
    for ssidx = 1:3
        r = [];
        ri = [];
        for bin = 1:6
            r = [r;rings_all_fast{midx,ssidx,bin}];
            ri = [ri;rings_inact_all_fast{midx,ssidx,bin}];
        end
        rings_all_remap{midx,ssidx} = mod(unique(r,'rows'),100000);
        rings_inact_all_remap{midx,ssidx} = mod(unique(ri,'rows'),100000);
    end
end


% load congruent rings
load rings.mat

ringsm = cell(1,3);
metadata = cell(1,3);
for midx=1:3
r1=[];
for bb=1:6
    r1=[r1;cell2mat(rings(midx,:,bb,1)')];
end
r2=[];
for bb=1:6
    r2=[r2;cell2mat(rings(midx,:,bb,2)')];
end
r1(:,midx+3)=1;r2(:,midx+3)=2;
ringsm{midx}=[r1(:,1:(midx+2));r2(:,1:(midx+2))];
metadata{midx} = zeros(length(ringsm{midx}),midx+3);
for ssidx=1:length(ringsm{midx})
    sessIdx=idivide(int32(ringsm{midx}(ssidx,1)),int32(100000));
    transIds=int32(rem(ringsm{midx}(ssidx,:),100000))';
    metadata{midx}(ssidx,1) = sessIdx;
    metadata{midx}(ssidx,2:end) = transIds';
end
end


ringInactsm = cell(1,3);
metadataInacgt = cell(1,3);
for midx=1:3
r1=[];
for bb=1:6
    r1=[r1;cell2mat(rings_inact(midx,:,bb,1)')];
end
r2=[];
for bb=1:6
    r2=[r2;cell2mat(rings_inact(midx,:,bb,2)')];
end
r1(:,midx+3)=1;r2(:,midx+3)=2;
ringInactsm{midx}=[r1(:,1:(midx+2));r2(:,1:(midx+2))];
metadataInact{midx} = zeros(length(ringInactsm{midx}),midx+3);
for ssidx=1:length(ringInactsm{midx})
    sessIdx=idivide(int32(ringInactsm{midx}(ssidx,1)),int32(100000));
    transIds=int32(rem(ringInactsm{midx}(ssidx,:),100000))';
    metadataInact{midx}(ssidx,1) = sessIdx;
    metadataInact{midx}(ssidx,2:end) = transIds';
end
end

pickedSessions = [23,46,92];
rings_congru_remap = cell(3,3);
rings_congru_inact_remap = cell(3,3);
for midx = 1:3
    for ssidx = 1:3
        rings_congru_remap{midx,ssidx} = metadata{midx}(metadata{midx}(:,1)==pickedSessions(ssidx),2:end);
        rings_congru_inact_remap{midx,ssidx} = metadataInact{midx}(metadataInact{midx}(:,1)==pickedSessions(ssidx),2:end);
    end
end

rings_congru_remap = cellfun(@(x) unique(x,'rows'),rings_congru_remap,'UniformOutput',false);
rings_congru_inact_remap = cellfun(@(x) unique(x,'rows'),rings_congru_inact_remap,'UniformOutput',false);
rings_incongru_remap = cellfun(@(x,y) setdiff(x,y,'rows'),rings_all_remap,rings_congru_remap,'UniformOutput',false);
rings_incongru_inact_remap = cellfun(@(x,y) setdiff(x,y,'rows'),rings_inact_all_remap,rings_congru_inact_remap,'UniformOutput',false);

%rings_missed = cellfun(@(x,y) setdiff(x,y,'rows'),rings_congru_remap,rings_all_remap,'UniformOutput',false); % should be empty

save('incongru_rings.mat','rings_incongru_remap','rings_incongru_inact_remap','-v7.3');



%if exist('delay_data','var') && delay_data
%    save('rings.mat','rings','-append');
%end
%if exist('delay_inact_data','var') && delay_inact_data
%    save('rings.mat','rings_inact','-append');
%end
%end


% if (exist('delay_shuf','var') && delay_shuf) || (exist('delay_shuf_inact','var') && delay_shuf_inact)
%     shufrpt=1000;
%     rings_shuf=cell(shufrpt,3,114,6,2);
%     rings_shuf_inact=cell(shufrpt,3,114,6,2);
%     for rpt=1:shufrpt
%         disp(rpt)
%         for bin=1:6
%             [shufchainS1,shufregS1,shufprefS1]=shuffle_conn_chain(fstr{bin}.conn_chain_S1,fstr{bin}.pair_chain,fstr{bin}.pair_reg,fstr{bin}.pref_pair);
%             [shufchainS2,shufregS2,shufprefS2]=shuffle_conn_chain(fstr{bin}.conn_chain_S2,fstr{bin}.pair_chain,fstr{bin}.pair_reg,fstr{bin}.pref_pair);
%             parfor I=1:114
%                 for midx=1:3
%                     msize=midx+2;
%                     lbound=100000*I;
%                     ubound=100000*(I+1);
%                     sel21=shufchainS1(:,1)>=lbound & shufchainS1(:,1)<ubound & diff(shufregS1,1,2);
%                     if nnz(sel21)>0
%                         if delay_shuf
%                             onering1=count_motif_congru(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),bin,msize,1);
%                             onering1=unique(flexsort(onering1),'rows');
%                             sel22=shufchainS2(:,1)>=lbound & shufchainS2(:,1)<ubound & diff(shufregS2,1,2);
%                             onering2=count_motif_congru(shufchainS2(sel22,:),shufregS2(sel22,:),shufprefS2(sel22,:),bin,msize,2);
%                             onering2=unique(flexsort(onering1),'rows');
%                             rings_shuf(rpt,midx,I,bin,:)={onering1,onering2};
%                         end
%                         if  delay_shuf_inact
%                             disp('Go')
%                             onering1=count_motif_congru_inact(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),msize,1);
%                             onering1=unique(flexsort(onering1),'rows');
%                             sel22=shufchainS2(:,1)>=lbound & shufchainS2(:,1)<ubound & diff(shufregS2,1,2);
%                             onering2=count_motif_congru_inact(shufchainS2(sel22,:),shufregS2(sel22,:),shufprefS2(sel22,:),msize,2);
%                             onering2=unique(flexsort(onering1),'rows');
%                             rings_shuf_inact(rpt,midx,I,bin,:)={onering1,onering2};
%                         end                        
%                     end
%                 end
%             end
%         end
%     end
%     if exist('delay_shuf','var') && delay_shuf
%         save('rings.mat','rings_shuf','-append');
%     end
%     if exist('delay_shuf_inact','var') && delay_shuf_inact
%         save('rings.mat','rings_shuf_inact','-append');
%     end
% end
% 
% %% baseline
% if exist('base_data','var') && base_data
%     base_rings=cell(3,114,2);
%     parfor I=1:114
%         for midx=1:3
%             msize=midx+2;
%             disp(I)
%             lbound=100000*I;
%             ubound=100000*(I+1);
%             sel11=fbase.conn_chain_S1(:,1)>=lbound & fbase.conn_chain_S1(:,1)<ubound & diff(fbase.reg_chain_S1,1,2);
%             if nnz(sel11)>0
%                 onering1=count_motif_congru_inact(fbase.conn_chain_S1(sel11,:),fbase.reg_chain_S1(sel11,:),fbase.pref_chain_S1(sel11,:),msize,1);
%                 onering1=unique(flexsort(onering1),'rows');
%                 sel12=fbase.conn_chain_S2(:,1)>=lbound & fbase.conn_chain_S2(:,1)<ubound & diff(fbase.reg_chain_S2,1,2);
%                 onering2=count_motif_congru_inact(fbase.conn_chain_S2(sel12,:),fbase.reg_chain_S2(sel12,:),fbase.pref_chain_S2(sel12,:),msize,2);
%                 onering2=unique(flexsort(onering2),'rows');
%                 base_rings(midx,I,:)={onering1,onering2};
%             end
%         end
%     end
%     save('rings.mat','base_rings','-append');
% end
% if exist('base_shuf','var') && base_shuf
%     shufrpt=1000;
%     base_rings_shuf=cell(shufrpt,3,114,2);
%     for rpt=1:shufrpt
%         disp(rpt)
%         [shufchainS1,shufregS1,shufprefS1]=shuffle_conn_chain(fbase.conn_chain_S1,fbase.pair_chain,fbase.pair_reg,fbase.pref_pair);
%         [shufchainS2,shufregS2,shufprefS2]=shuffle_conn_chain(fbase.conn_chain_S2,fbase.pair_chain,fbase.pair_reg,fbase.pref_pair);
%         parfor I=1:114
%             for midx=1:3
%                 msize=midx+2;
%                 lbound=100000*I;
%                 ubound=100000*(I+1);
%                 sel21=shufchainS1(:,1)>=lbound & shufchainS1(:,1)<ubound & diff(shufregS1,1,2);
%                 if nnz(sel21)>0
%                     onering1=count_motif_congru_inact(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),msize,1);
%                     onering1=unique(flexsort(onering1),'rows');
%                     sel22=shufchainS2(:,1)>=lbound & shufchainS2(:,1)<ubound & diff(shufregS2,1,2);
%                     onering2=count_motif_congru_inact(shufchainS2(sel22,:),shufregS2(sel22,:),shufprefS2(sel22,:),msize,2);
%                     onering2=unique(flexsort(onering1),'rows');
%                     base_rings_shuf(rpt,midx,I,:)={onering1,onering2};
%                 end
%             end
%         end
%     end
%     save('rings.mat','base_rings_shuf','-append');
% end


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