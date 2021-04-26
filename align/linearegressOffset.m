% linear regression and then compare shrinkage factor and decide probe offset
EphysLabelFolder = cd;
EphysLabelLists = dir('*.csv');
savepath = 'I:\NP histology\probe revise\figure';
load('I:\NP histology\shrinkF.mat');
load('I:\NP histology\histoLengthTip.mat');
load('I:\NP histology\insertLengthTip.mat');
revisedLength = zeros(size(insertLengthMat));
revisedLength(:,1) = insertLengthMat(:,1);
insertLengthMat(insertLengthMat > 3840) = 3840;
for tableIdx = 1:length(EphysLabelLists)
    C = strsplit(EphysLabelLists(tableIdx).name,'.');
    meta = strsplit(C{1},'-');
    miceId = meta{1};
    probeId = meta{2};
    insertLength = insertLengthMat(shrinkageMat(:,1) == str2double(miceId),str2double(probeId)+1);
    shrinkage_factor = shrinkageMat(shrinkageMat(:,1) == str2double(miceId),str2double(probeId)+1);
    probe_histo_length_tip = histoLengthMat(shrinkageMat(:,1) == str2double(miceId),str2double(probeId)+1);
    histo_active_length = insertLength*shrinkage_factor;
    depthListPx = readtable(EphysLabelLists(tableIdx).name);
    histoRanges = [range(depthListPx.zha_y_histo_probe) range(depthListPx.xd_y_histo_probe)];
    ephysRanges = [range(depthListPx.zha_y_ephys),range(depthListPx.xd_y_ephys)];
    if range(depthListPx.zha_y_histo_probe) > 0 && range(depthListPx.xd_y_histo_probe) > 0
        px2um_histo_psudo = insertLength / mean([range(depthListPx.zha_y_histo_probe),range(depthListPx.xd_y_histo_probe)]);
        px2um_histo = histo_active_length / mean([range(depthListPx.zha_y_histo_probe),range(depthListPx.xd_y_histo_probe)]);
        px2um_ephys = 3840 / mean([range(depthListPx.zha_y_ephys),range(depthListPx.xd_y_ephys)]);
    else
        px2um_histo_psudo = insertLength / histoRanges(histoRanges > 0);
        px2um_histo = histo_active_length / histoRanges(histoRanges > 0);
        px2um_ephys = 3840 / ephysRanges(ephysRanges > 0);
    end
    try
        depth_histo_zha_px = sort(depthListPx.zha_y_histo_probe(depthListPx.zha_y_histo_probe ~= 0));
        depth_histo_zha = px2um_histo_psudo * (depth_histo_zha_px(2:end-1) - depth_histo_zha_px(1));
    catch
        depth_histo_zha = [];
    end
    try
        depth_histo_xd_px = sort(depthListPx.xd_y_histo_probe(depthListPx.xd_y_histo_probe ~= 0));
        depth_histo_xd = px2um_histo_psudo * (depth_histo_xd_px(2:end-1) - depth_histo_xd_px(1));
    catch
        depth_histo_xd = [];
    end
    try
        depth_ephys_zha_px = sort(depthListPx.zha_y_ephys(depthListPx.zha_y_ephys ~= 0));
        depth_ephys_zha = px2um_ephys * (depth_ephys_zha_px(2:end-1) - depth_ephys_zha_px(1));
    catch
        depth_ephys_zha = [];
    end
    try
        depth_ephys_xd_px = sort(depthListPx.xd_y_ephys(depthListPx.xd_y_ephys ~= 0));
        depth_ephys_xd = px2um_ephys * (depth_ephys_xd_px(2:end-1) - depth_ephys_xd_px(1));
    catch
        depth_ephys_xd = [];
    end
        plot(depth_histo_zha,depth_ephys_zha,'r*',depth_histo_xd,depth_ephys_xd,'b*');
        hold on;
        X = [depth_histo_zha;depth_histo_xd];
        X = [ones(size(X)) X];
        Y = [depth_ephys_zha;depth_ephys_xd];
        b = regress(Y,X);
        y = b(2)*X(:,2) + b(1);
        plot(X(:,2),y,'linewidth',2);
        hold on;
        line([0 3840],[0,3840*shrinkage_factor],'Color',[0.2,0.2,0.2],'LineStyle','--');
        ax = gca;
        set(gca,'XLim',[0 3840],'YLim',[0 3840]);
        saveas(gcf,[savepath '\' C{1} '.png']);
        pause(0.05);
        close;
        estimated_max_depth = probe_histo_length_tip + 3840 - histo_active_length - b(1);
        revisedLength(revisedLength(:,1) == str2double(miceId),str2double(probeId)+1) = estimated_max_depth;
    save('I:\NP histology\reviseLengthTip.mat','revisedLength');
end
    
    
    
    
    
    