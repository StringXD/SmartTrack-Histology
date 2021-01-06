% Selectivity index of all selective neurons
cid_list=h5read('transient_6.hdf5','/cluster_id');
sus_trans=h5read('transient_6.hdf5','/sus_trans');
reg_list=h5read('transient_6.hdf5','/reg');


reg_logic=~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'Unlabeled') &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'aco')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'act')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'alv')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'amc')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'ar')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'bsc')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'ccb')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'ccg')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'ccs')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'cing')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'dhc')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'ec')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'ee')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'em')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'fa')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'fi')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'fiber')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'fp')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'fr')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'int')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'opt')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'or')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'root')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'scwm')&~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'sm')...,
    &~strcmp(regexp(reg_list,'(\w|\\|-)*','match','once'),'st');

sus_logic = sus_trans(:,1);
trans_logic = sus_trans(:,2) | sus_trans(:,4);


load FR_modulated_20210104.mat

% firing rate
frdataFAsus = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==4&y(:,9)==1&y(:,10)==1,:),FR_modulated(reg_logic&sus_logic,1),FR_modulated(reg_logic&sus_logic,2),'UniformOutput',false);
frdataFBsus = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==8&y(:,9)==1&y(:,10)==1,:),FR_modulated(reg_logic&sus_logic,1),FR_modulated(reg_logic&sus_logic,2),'UniformOutput',false);
frdataFAtrans = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==4&y(:,9)==1&y(:,10)==1,:),FR_modulated(reg_logic&trans_logic,1),FR_modulated(reg_logic&trans_logic,2),'UniformOutput',false);
frdataFBtrans = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==8&y(:,9)==1&y(:,10)==1,:),FR_modulated(reg_logic&trans_logic,1),FR_modulated(reg_logic&trans_logic,2),'UniformOutput',false);
meanFrDataFAsus = cellfun(@mean,frdataFAsus,'UniformOutput',false);
meanFrDataFBsus = cellfun(@mean,frdataFBsus,'UniformOutput',false);
meanFrDataFAtrans = cellfun(@mean,frdataFAtrans,'UniformOutput',false);
meanFrDataFBtrans = cellfun(@mean,frdataFBtrans,'UniformOutput',false);
meanFrDataFAsus = cell2mat(meanFrDataFAsus);
meanFrDataFBsus = cell2mat(meanFrDataFBsus);
meanFrDataFAtrans = cell2mat(meanFrDataFAtrans);
meanFrDataFBtrans = cell2mat(meanFrDataFBtrans);
% selectivity index
unitSelIdx6sSus = (meanFrDataFAsus - meanFrDataFBsus)./(meanFrDataFAsus + meanFrDataFBsus);
unitSelIdx6sSus(isnan(unitSelIdx6sSus)) = 0;
unitSelIdx6sTrans = (meanFrDataFAtrans - meanFrDataFBtrans)./(meanFrDataFAtrans + meanFrDataFBtrans);
unitSelIdx6sTrans(isnan(unitSelIdx6sTrans)) = 0;
% smooth
unitSelIdx6sSus = num2cell(unitSelIdx6sSus,2);
sunitSelIdx6sSus = cellfun(@(x) smooth(x)',unitSelIdx6sSus,'UniformOutput',false);
unitSelIdx6sTrans = num2cell(unitSelIdx6sTrans,2);
sunitSelIdx6sTrans = cellfun(@(x) smooth(x)',unitSelIdx6sTrans,'UniformOutput',false);
sunitSelIdx6sSus = cell2mat(sunitSelIdx6sSus);
sunitSelIdx6sTrans = cell2mat(sunitSelIdx6sTrans);
% find significant bins
isbinsig6sSus = zeros(size(sunitSelIdx6sSus,1),60);
isbinsig6sTrans = zeros(size(sunitSelIdx6sTrans,1),60);
for i = 1:size(isbinsig6sSus,1)
    for j = 1:60
        [~,isbinsig6sSus(i,j)] = ranksum(frdataFAsus{i}(:,j),frdataFBsus{i}(:,j));
    end
end
for i = 1:size(isbinsig6sTrans,1)
    for j = 1:60
        [~,isbinsig6sTrans(i,j)] = ranksum(frdataFAtrans{i}(:,j),frdataFBtrans{i}(:,j));
    end
end
% mask
sigimg6sSus = sunitSelIdx6sSus(:,1:60).*isbinsig6sSus;
sigimg6sTrans = sunitSelIdx6sTrans(:,1:60).*isbinsig6sTrans;
% Classify by significant bin num, sort by peak bin
% jump to heatMapSelGM.m to draw styled heatmap
sigbinNum6s = sum(isbinsig6sTrans(:,17:40),2);

sigimg6sTransPart1Pos = sigimg6sTrans(sigbinNum6s <= 4 & mean(sigimg6sTrans,2) >= 0,:);
sigimg6sTransPart2Pos = sigimg6sTrans(sigbinNum6s > 4 & sigbinNum6s <= 8 & mean(sigimg6sTrans,2) >= 0,:);
sigimg6sTransPart3Pos = sigimg6sTrans(sigbinNum6s > 8 & sigbinNum6s <= 12 & mean(sigimg6sTrans,2) >= 0,:);
sigimg6sTransPart4Pos = sigimg6sTrans(sigbinNum6s > 12 & sigbinNum6s <= 16 & mean(sigimg6sTrans,2) >= 0,:);
sigimg6sTransPart5Pos = sigimg6sTrans(sigbinNum6s > 16 & sigbinNum6s <= 20 & mean(sigimg6sTrans,2) >= 0,:);
sigimg6sTransPart6Pos = sigimg6sTrans(sigbinNum6s > 20 & mean(sigimg6sTrans,2) >= 0,:);

sigimg6sTransPart1Neg = sigimg6sTrans(sigbinNum6s <= 4 & mean(sigimg6sTrans,2) < 0,:);
sigimg6sTransPart2Neg = sigimg6sTrans(sigbinNum6s > 4 & sigbinNum6s <= 8 & mean(sigimg6sTrans,2) < 0,:);
sigimg6sTransPart3Neg = sigimg6sTrans(sigbinNum6s > 8 & sigbinNum6s <= 12 & mean(sigimg6sTrans,2) < 0,:);
sigimg6sTransPart4Neg = sigimg6sTrans(sigbinNum6s > 12 & sigbinNum6s <= 16 & mean(sigimg6sTrans,2) < 0,:);
sigimg6sTransPart5Neg = sigimg6sTrans(sigbinNum6s > 16 & sigbinNum6s <= 20 & mean(sigimg6sTrans,2) < 0,:);
sigimg6sTransPart6Neg = sigimg6sTrans(sigbinNum6s > 20 & mean(sigimg6sTrans,2) < 0,:);

sigimg6sSusPos = sigimg6sSus(mean(sigimg6sSus,2) >= 0,:);
sigimg6sSusNeg = sigimg6sSus(mean(sigimg6sSus,2) < 0,:);

ordered6sTransPos = [sortByPeakTimeAsc(sigimg6sTransPart1Pos,sigimg6sTransPart1Pos(:,13:44),2,'ascend');sortByPeakTimeAsc(sigimg6sTransPart2Pos,sigimg6sTransPart2Pos(:,13:44),2,'ascend');...
    sortByPeakTimeAsc(sigimg6sTransPart3Pos,sigimg6sTransPart3Pos(:,13:44),2,'ascend');sortByPeakTimeAsc(sigimg6sTransPart4Pos,sigimg6sTransPart4Pos(:,13:44),2,'ascend');...
    sortByPeakTimeAsc(sigimg6sTransPart5Pos,sigimg6sTransPart5Pos(:,13:44),2,'ascend')];
ordered6sSusPos = sortByPeakTimeAsc(sigimg6sSus(mean(sigimg6sSus,2) >= 0,:),sigimg6sSus(mean(sigimg6sSus,2) >= 0,:),2,'ascend');
ordered6sTransNeg = [sortByPeakTimeAsc(sigimg6sTransPart1Neg,sigimg6sTransPart1Neg(:,13:44),2,'descend');sortByPeakTimeAsc(sigimg6sTransPart2Neg,sigimg6sTransPart2Neg(:,13:44),2,'descend');...
    sortByPeakTimeAsc(sigimg6sTransPart3Neg,sigimg6sTransPart3Neg(:,13:44),2,'descend');sortByPeakTimeAsc(sigimg6sTransPart4Neg,sigimg6sTransPart4Neg(:,13:44),2,'descend');...
    sortByPeakTimeAsc(sigimg6sTransPart5Neg,sigimg6sTransPart5Neg(:,13:44),2,'descend');sigimg6sTransPart6Neg];
ordered6sSusNeg = sortByPeakTimeAsc(sigimg6sSus(mean(sigimg6sSus,2) < 0,:),sigimg6sSus(mean(sigimg6sSus,2) < 0,:),2,'descend');
% Calculate the length of matrices so as to draw a dashed line
dash6s = size([ordered6sTransPos;ordered6sTransNeg],1)+0.5;

map6s = [ordered6sTransPos;ordered6sTransNeg;ordered6sSusPos;ordered6sSusNeg];
map = [[[0:0.1:0.9;0:0.1:0.9]' ones(10,1)];[1,1,1];[ones(10,1) [0.9:-0.1:0;0.9:-0.1:0]']];
colormap(map);


imagesc(map6s);
c = colorbar;
set(c,'Ticks',[-1,0,1]);
c.Label.String = 'Selectivity Index';
caxis([-1 1]);
set(gca,'XTickLabel',{'0','5','10'},'XTick',[12.5,32.5,52.5]);
set(gca,'YDir','normal'); 
set(gca,'YTickLabel',{'0','2000','4000','6000','8000'},'YTick',[0.5,2000.5,4000.5,6000.5,8000.5]);
hold on;
ax = gca;
plot(ax.XLim,[dash6s dash6s],'k--','LineWidth',1);
hold on;
plot([12.5,12.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([16.5,16.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([40.5,40.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([44.5,44.5],ax.YLim,'k--','LineWidth',2);
ylabel('Neuron ID');
xlabel('Time (s)');
title('6s Sustain/Transient Neurons');