%
load spike2FCActiveRaw.mat

feedRateS1 = [];
feedRateS2 = [];
feedRateCAS1 = [];
feedRateCAS2 = [];
for ssidx = 1:length(sumSum)
    if isempty(sumSum{ssidx,1})
        continue;
    end
    s1data = sumSum{ssidx,1};
    s2data = sumSum{ssidx,2};
    feedRateCAS1 = [feedRateCAS1;nanmean((s1data(:,1:6)+s1data(:,7:12))./s1data(:,55:60))];
    feedRateCAS2 = [feedRateCAS2;nanmean((s2data(:,1:6)+s2data(:,7:12))./s2data(:,55:60))];
    feedRateS1 = [feedRateS1;nanmean(s1data(:,49:54)./s1data(:,55:60))];
    feedRateS2 = [feedRateS2;nanmean(s2data(:,49:54)./s2data(:,55:60))];
end
feedRateCA = [feedRateCAS1;feedRateCAS2];
feedRate = [feedRateS1;feedRateS2];
meanFeedRateCA = mean(feedRateCA);
meanFeedRate = mean(feedRate);
semFeedRateCA = std(feedRateCA)/sqrt(size(feedRateCA,1)-1);
semFeedRate = std(feedRate)/sqrt(size(feedRate,1)-1);

%h = plot(1:6,meanFeedRateCA,1:6,meanFeedRate);
hold on;
errorbar(1:6,meanFeedRateCA,semFeedRateCA);
errorbar(1:6,meanFeedRate,semFeedRate);

ylabel('Functional coupling spike ratio');
xlabel('Time(s)')
legend

% congru vs incongru
feedRateCS1 = [];
feedRateIS1 = [];
feedRateCS2 = [];
feedRateIS2 = [];
for ssidx = 1:length(sumSum)
    if isempty(sumSum{ssidx,1})
        continue;
    end
    s1data = sumSum{ssidx,1};
    s2data = sumSum{ssidx,2};
    feedRateCS1 = [feedRateCS1;nanmean(s1data(:,37:42)./s1data(:,55:60),'all')];
    feedRateCS2 = [feedRateCS2;nanmean(s2data(:,37:42)./s2data(:,55:60),'all')];
    feedRateIS1 = [feedRateIS1;nanmean(s1data(:,43:48)./s1data(:,55:60),'all')];
    feedRateIS2 = [feedRateIS2;nanmean(s2data(:,43:48)./s2data(:,55:60),'all')];
end
bar([mean(feedRateCS1),mean(feedRateCS2),mean(feedRateIS1),mean(feedRateIS2)]);