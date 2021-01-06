%% ring rounds exponential decay
results = cell(1,3);
for midx = 1:3
    load(sprintf('turnStatLen%dFull.mat',midx+2));
    results{midx} = [];
    for j = 1:length(turnStat)
        results{midx} = [results{midx},cell2mat(turnStat{j,1}')];
    end
end
rStat3 = results{1}/3;
rStat4 = results{2}/4;
rStat5 = results{3}/5;
N3 = histcounts(rStat3,5/6:1/3:31/6);
N4 = histcounts(rStat4,7/8:1/4:41/8);
N5 = histcounts(rStat5,9/10:1/5:51/10);
% count -> frequency
load tnumList.mat
totalTrialNum = sum(cell2mat(tnumList));
plot(1:1/3:5,N3/(totalTrialNum*6));
hold on;
plot(1:1/4:5,N4/(totalTrialNum*6));
plot(1:1/5:5,N5/(totalTrialNum*6));
ax = gca;
ax.XTick = 1:5;
ax.XLim = [1 5];
set(gca,'yscale','log')
ylabel('Frequency(Hz)');
xlabel('Rounds of Ring Activity');
legend