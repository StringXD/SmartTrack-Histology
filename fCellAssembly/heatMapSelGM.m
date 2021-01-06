% heatMap of Selective Neuron Firing Rate
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

sel_logic = sus_trans(:,1) | sus_trans(:,2) | sus_trans(:,4) ;


load FR_modulated_20210104.mat


frdataFA = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==4&y(:,9)==1&y(:,10)==1,:),FR_modulated(reg_logic&sel_logic,1),FR_modulated(reg_logic&sel_logic,2),'UniformOutput',false);
frdataFB = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==8&y(:,9)==1&y(:,10)==1,:),FR_modulated(reg_logic&sel_logic,1),FR_modulated(reg_logic&sel_logic,2),'UniformOutput',false);
baselineFA = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==4&y(:,9)==1&y(:,10)==1,1:10),FR_modulated(reg_logic&sel_logic,1),FR_modulated(reg_logic&sel_logic,2),'UniformOutput',false);
baselineFB = cellfun(@(x,y) x(y(:,8)==6&y(:,5)==8&y(:,9)==1&y(:,10)==1,1:10),FR_modulated(reg_logic&sel_logic,1),FR_modulated(reg_logic&sel_logic,2),'UniformOutput',false);
baselineFAmean = cellfun(@(x) mean(x,'all'),baselineFA,'UniformOutput',false);
baselineFBmean = cellfun(@(x) mean(x,'all'),baselineFB,'UniformOutput',false);
baselineFAstd = cellfun(@(x) std(x,[],'all'),baselineFA,'UniformOutput',false);
baselineFBstd = cellfun(@(x) std(x,[],'all'),baselineFB,'UniformOutput',false);
ZscoreA = cellfun(@(x,y,z) mean((x-y)/z),frdataFA,baselineFAmean,baselineFAstd,'UniformOutput',false);
ZscoreB = cellfun(@(x,y,z) mean((x-y)/z),frdataFB,baselineFBmean,baselineFBstd,'UniformOutput',false);
% smooth
sZscoreA = cellfun(@(x) smooth(x)',ZscoreA,'UniformOutput',false);
sZscoreA = cell2mat(sZscoreA);
sZscoreB = cellfun(@(x) smooth(x)',ZscoreB,'UniformOutput',false);
sZscoreB = cell2mat(sZscoreB);

% Sort by mean S1-S2 firing rate
meanDiff6s = mean(sZscoreA(:,17:40) - sZscoreB(:,17:40),2);
[~,order] = sort(meanDiff6s,'descend');

subplot(1,2,1);
colormap(jet);
imagesc(sZscoreA(order,1:60));
c = colorbar;
caxis([-2 2]);
set(c,'Ticks',[-2,0,2]);
c.Label.String = 'Z-scored firing rate';
set(gca,'XTickLabel',{'0','5','10'},'XTick',[12.5,32.5,52.5]);
set(gca,'YDir','normal'); 
set(gca,'YTickLabel',{'0','2000','4000','6000','8000'},'YTick',[0.5,2000.5,4000.5,6000.5,8000.5]);
hold on;
ax = gca;
plot([12.5,12.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([16.5,16.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([40.5,40.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([44.5,44.5],ax.YLim,'k--','LineWidth',2);
ylabel('Neuron ID');
title('6s Delay Follow Sample 1');
subplot(1,2,2);
colormap(jet);
imagesc(sZscoreB(order,1:60));
c = colorbar;
caxis([-2 2]);
set(c,'Ticks',[-2,0,2]);
c.Label.String = 'Z-scored firing rate';
set(gca,'XTickLabel',{'0','5','10'},'XTick',[12.5,32.5,52.5]);
set(gca,'YDir','normal'); 
set(gca,'YTickLabel',{'0','2000','4000','6000','8000'},'YTick',[0.5,2000.5,4000.5,6000.5,8000.5]);
hold on;
ax = gca;
plot([12.5,12.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([16.5,16.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([40.5,40.5],ax.YLim,'k--','LineWidth',2);
hold on;
plot([44.5,44.5],ax.YLim,'k--','LineWidth',2);
title('6s Delay Follow Sample 2');
