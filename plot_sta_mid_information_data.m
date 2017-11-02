function [x, y] = plot_sta_mid_information_data(info)
% plot_sta_mid_information -compare information values for sta and mids,
% and also determine degree of nonlinearity and synergy
%
%
% info - struct array holding information data
%
% The plots made are:
%
% (1) I(1) versus I(1,2)
% (2) I(1) versus I(2)
% (3) I(1,2) versus I(1)+I(2)
% (4) I(1) versus I(sta)
%
% 
%
% You can make pretty pre-publication eps plots via the following commands:
%
% exportfig(gcf,'mid_info_histograms.eps', 'color', 'gray', 'height', 4, 'width', 6, 'fontmode', 'fixed', 'fontsize', 8);
% exportfig(gcf,'mid_info_scatter.eps', 'color', 'gray', 'height', 4, 'width', 6, 'fontmode', 'fixed', 'fontsize', 8);

% caa 6/1/06

set(0,'defaultAxesFontName', 'Palatino')
set(0,'defaultAxesFontSize', 10)

info_sta = [];
info1 = [];
info2 = [];
info_both = [];

for i = 1:length(info)

   info_sta = [info_sta mean(info(i).sta.information)];
   info1 = [info1 mean(info(i).mid1.information)];
   info2 = [info2 mean(info(i).mid2.information)];
   info_both = [info_both mean(info(i).mid12.information)];

end % (for i)


% Error check to get rid of bad INF values:

index = ~isinf(info_sta) & ~isinf(info1) & ~isinf(info2) & ~isinf(info_both);

info_sta = info_sta(index);
info1 = info1(index);
info2 = info2(index);
info_both = info_both(index);


%exportfig(gcf,'mid_info_histograms.eps', 'color', 'gray', 'height', 4, 'width', 6, 'fontmode', 'fixed', 'fontsize', 8);
%exportfig(gcf,'mid_info_scatter.eps', 'color', 'rgb', 'height', 4, 'width', 6, 'fontmode', 'fixed', 'fontsize', 8);

% Compare the information values on linear axes:
% ----------------------------------------------

f1 = figure;

subplot(2,2,1)
hold on
xmax = max(info_both);
plot(info1, info_both, 'o')
fplot('x', [0 max(info_both)], 'k-')
fplot('2*x', [0 max(info_both)], 'k-')
xlim([0 1.1*xmax]);
ylim([0 1.1*xmax]);
axis square
xlabel('info(1 dir)')
ylabel('info(1 dir and 2 dir)')
title('difference between 1dir-model and 2dir-model')

subplot(2,2,2)
hold on
xmax = max([max(info1) max(info2)]);
plot(info1, info2,'o')
fplot('x', [0 xmax], 'k-')
xlim([0 1.1*xmax]);
ylim([0 1.1*xmax]);
axis square

xlabel('info(1 dir)')
ylabel('info(2 dir)')
title('relative importance of 2dir')

ind_synergy = find(info_both > info1+info2);
ind_strong_synergy = find(info_both > (info1+info2)*1.5);
length(ind_strong_synergy);
ind_not_synergy = find(info_both < info1+info2);
Nsynergy = length(ind_synergy);
info_both(ind_strong_synergy)./(info1(ind_strong_synergy)+info2(ind_strong_synergy));


subplot(2,2,3)
hold off
plot(info1(ind_synergy)+info2(ind_synergy), info_both(ind_synergy),'ro')
hold on
plot(info1(ind_not_synergy)+info2(ind_not_synergy),info_both(ind_not_synergy),'o')
fplot('x',[0 max([max(info1),max(info2),max(info_both)])])
axis square
xlabel('info(1dir)+info(2dir)')
ylabel('info(1dir + 2dir)')
title('synergy between directions - selectivity?')


subplot(2,2,4)
hold on;
plot(info1, info_sta, 'o')
fplot('x',[0 max(max(info1),max(info_sta))])
xlim([0 max(max(info1),max(info_sta))])
ylim([0 max(max(info1),max(info_sta))])
axis square
title('difference between  sta and 1dir')
xlabel('info(sta)')
ylabel('info(1dir)')

close(f1);



% Compare the information values on log axes:
% ----------------------------------------------


% Compare the first filter's information to the total information
figure;
hold on
xmax = max([max(info1) max(info_both)]);
xmin = min([min(info1) min(info_both)]);
plot(info1, info_both, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');

fplot('x', [0.5*xmin 1.75*xmax], 'k-')
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);

set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
axis square
xlabel('I(v_1)')
ylabel('I(v_1,v_2)')
title('Degree of Feature Selectivity')


% Compare the second filter's information to the total information
figure;
hold on
xmax = max([max(info2) max(info_both)]);
xmin = min([min(info2) min(info_both)]);
plot(info2, info_both, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');

fplot('x', [0.5*xmin 1.75*xmax], 'k-')
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);

set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
axis square
xlabel('I(v_2)')
ylabel('I(v_1,v_2)')
title('Contribution of Second Filter')


figure;
hold on
xmax = max([max(info1) max(info2)]);
xmin = min([min(info1) min(info2)]);
plot(info1, info2, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
fplot('x', [0.5*xmin 1.75*xmax], 'k-')
% xlim([0.01 1.75*xmax]);
% ylim([0.01 1.75*xmax]);
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);
set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
axis square

xlabel('I(v_1)')
ylabel('I(v_2)')
title('Relative Importance of 2dir')

ind_synergy=find(info_both>info1+info2);
ind_strong_synergy=find(info_both>(info1+info2)*1.5);
length(ind_strong_synergy);
ind_not_synergy=find(info_both<info1+info2);
Nsynergy=length(ind_synergy);
info_both(ind_strong_synergy)./(info1(ind_strong_synergy)+info2(ind_strong_synergy));


figure;
hold on;
xmax = max([max(info1+info2) max(info_both)]);
xmin = min([min(info1+info2) min(info_both)]);
plot(info1+info2, info_both,'ko', 'markersize', 3);
% plot(info1(ind_synergy)+info2(ind_synergy),info_both(ind_synergy),'ro', 'markersize', 3);
% hp = plot(info1(ind_not_synergy)+info2(ind_not_synergy),info_both(ind_not_synergy),'ko', 'markersize', 3);
fplot('x', [0.5*xmin 1.75*xmax], 'k-')
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);
set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
axis square;
xlabel('I(v_1) + I(v_2)')
ylabel('I(v_1,v_2)')


figure;
hold on;
xmax = max([max(info1) max(info_sta)]);
xin = min([min(info1) min(info_sta)]);
plot(info1, info_sta, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');

fplot('x', [0.5*xmin 1.75*xmax], 'k-')
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);

set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
axis square
xlabel('I(v_1)')
ylabel('I(sta)')




% %===========================================================
% % Plot the histograms using #units
% %===========================================================
% 
% % [info_sta info1 info2 (info1+info2) info_both]
% 
% info1_info12 = 100 * info1 ./ info_both;
% info1_info1_info2 = 100 * info1 ./ (info1 + info2);
% info12_info1_info2 = 100 * info_both ./ (info1 + info2);
% info2_info1 = 100 * info2 ./ info1;
% infosta_info1 = 100 * info_sta ./ info1;
% 
% infosta_info12 = 100 * info_sta ./ info_both;
% 
% 
% figure;
% % subplot(2,2,1);
% bins = linspace(0,100,11);
% count = histc(info1_info12, bins);
% hb = bar(bins, count, 'histc');
% set(hb, 'facecolor', [0.75 0.75 0.75])
% axis([-10 110 0 1.05*max(count)]);
% set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
% set(gca,'tickdir', 'out');
% xlabel('Percent(%)', 'fontsize', 10);
% ylabel('#units', 'fontsize', 10);
% title('100 * Info(1) / Info(1,2)', 'fontsize', 10);
% set(gca,'fontsize', 8);
% 
% 
% figure;
% % subplot(2,2,1);
% bins = linspace(0,100,11);
% count = histc(info1_info1_info2, bins);
% hb = bar(bins, count, 'histc');
% set(hb, 'facecolor', [0.75 0.75 0.75])
% axis([-10 110 0 1.05*max(count)]);
% set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
% set(gca,'tickdir', 'out');
% xlabel('Percent(%)', 'fontsize', 10);
% ylabel('#units', 'fontsize', 10);
% title('100 * Info(1) / ( Info(1) + Info(2) )', 'fontsize', 10);
% set(gca,'fontsize', 8);
% 
% 
% figure;
% % subplot(2,2,2);
% bins = linspace(0,300,13);
% count = histc(info12_info1_info2, bins);
% hb = bar(bins, count, 'histc');
% set(hb, 'facecolor', [0.75 0.75 0.75])
% axis([100*[-0.1 3.1] 0 1.05*max(count)]);
% set(gca,'xtick', [0:50:300], 'xticklabel', [0:50:300]);
% set(gca,'tickdir', 'out');
% xlabel('Percent(%)', 'fontsize', 10);
% ylabel('#units', 'fontsize', 10);
% title('100 * Info(1,2) / (Info(1)+Info(2))', 'fontsize', 10);
% set(gca,'fontsize', 8);
% 
% 
% figure;
% % subplot(2,2,2);
% bins = linspace(0,100,11);
% count = histc(info2_info1, bins);
% hb = bar(bins, count, 'histc');
% set(hb, 'facecolor', [0.75 0.75 0.75])
% axis([-10 110 0 1.05*max(count)]);
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
% set(gca,'tickdir', 'out');
% xlabel('Percent(%)', 'fontsize', 10);
% ylabel('#units', 'fontsize', 10);
% title('100 * Info(2) / Info(1)', 'fontsize', 10);
% set(gca,'fontsize', 8);
% 
% 
% figure;
% % subplot(2,2,2);
% bins = linspace(0,100,11);
% count = histc(infosta_info1, bins);
% hb = bar(bins, count, 'histc');
% set(hb, 'facecolor', [0.75 0.75 0.75])
% axis([-10 110 0 1.05*max(count)]);
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
% set(gca,'tickdir', 'out');
% xlabel('Percent(%)', 'fontsize', 10);
% ylabel('#units', 'fontsize', 10);
% title('100 * Info(STA) / Info(1)', 'fontsize', 10);
% set(gca,'fontsize', 8);
% 
% 
% close all;




%===========================================================
% Plot the histograms but using proportion instead of #units
%===========================================================


% [info_sta info1 info2 (info1+info2) info_both]

info1_info12 = 100 * info1 ./ info_both;
info2_info12 = 100 * info2 ./ info_both;
info1_info1_info2 = 100 * info1 ./ (info1 + info2);
info12_info1_info2 = 100 * info_both ./ (info1 + info2);
info2_info1 = 100 * info2 ./ info1;
infosta_info1 = 100 * info_sta ./ info1;

infosta_info12 = 100 * info_sta ./ info_both;


figure;
bins = linspace(0,100,11);
count = histc(info1_info12, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
set(gca,'tickdir', 'out');
xlabel('100 * I(v_1) / I(v_1,v_2)');
ylabel('Percent (%)');



mean(info1_info12)
std(info1_info12)

pause




figure;
bins = linspace(0,50,11);
count = histc(info2_info12, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-5 55 0 1.05*max(count)]);
set(gca,'xtick', 0:10:50, 'xticklabel', 0:10:50);
set(gca,'tickdir', 'out');
xlabel('100 * I(v_2) / I(v_1,v_2)');
ylabel('Percent (%)');


figure;
bins = linspace(0,100,11);
count = histc(info1_info1_info2, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
set(gca,'tickdir', 'out');
xlabel('100 * I(v_1) / ( I(v_1) + I(v_2) )');
ylabel('Percent (%)');


figure;
bins = linspace(0,300,13);
count = histc(info12_info1_info2, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([100*[-0.1 3.1] 0 1.05*max(count)]);
set(gca,'xtick', [0:50:300], 'xticklabel', [0:50:300]);
set(gca,'tickdir', 'out');
xlabel('100 * I(v_1,v_2) / (I(v_1)+I(v_2))');
ylabel('Percent (%)');


figure;
bins = linspace(0,100,11);
count = histc(info2_info1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'tickdir', 'out');
xlabel('100 * I(v_2) / I(v_1)');
ylabel('Percent (%)');




figure;
bins = linspace(0,100,11);
count = histc(infosta_info1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'tickdir', 'out');
xlabel('100 * I(sta) / I(v_1)');
ylabel('Percent (%)');



%===========================================================
% Plot the second filter index versus the synergy index
%===========================================================
second_filter_index = 100 * info2 ./ info_both;
synergy_index = 100 * info_both ./ (info1 + info2);

xmin = 0;
xmax = 50;
ymin = 70;
ymax = 210;

figure;
xfit = min(second_filter_index):max(second_filter_index);
yfit = -2.209 .* xfit + 176.1;
hold on;
plot(second_filter_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
plot(xfit, yfit, 'k-');
xlim([xmin xmax]);
ylim([ymin ymax]);
set(gca,'tickdir', 'out');
[r, p] = corrcoef(second_filter_index, synergy_index);
title(sprintf('r=%.3f, p=%.5f', r(1,2), p(1,2)))
xlabel('100 I(v_2) / I(v_1,v_2)')
ylabel('100 I(v_1,v_2) / I(v_1)+I(v_2)')




%===========================================================
% Plot the second filter index versus the synergy index
%===========================================================
second_filter_index = 100 * info2 ./ (info1 + info2);
synergy_index = 100 * info_both ./ (info1 + info2);

xmin = 0;
xmax = 50;
ymin = 70;
ymax = 210;

figure;
xfit = min(second_filter_index):max(second_filter_index);
yfit = -1.988 .* xfit + 181.9;
hold on;
plot(second_filter_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
plot(xfit, yfit, 'k-');
xlim([xmin xmax]);
ylim([ymin ymax]);
set(gca,'tickdir', 'out');
[r, p] = corrcoef(second_filter_index, synergy_index);
title(sprintf('r=%.3f, p=%.5f', r(1,2), p(1,2)))
xlabel('100 I(v_2) / I(v_1)+I(v_2)')
ylabel('100 I(v_1,v_2) / I(v_1)+I(v_2)')


%===================================================================
% Plot the second filter contribution index versus the synergy index
%===================================================================
second_filter_index = 100 * info2 ./ info1;
synergy_index = 100 * info_both ./ (info1 + info2);

xmin = 0;
xmax = 75;
ymin = 70;
ymax = 210;

figure;
xfit = min(second_filter_index):max(second_filter_index);
yfit = -1.184 .* xfit + 173.4;
hold on;
plot(second_filter_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
plot(xfit, yfit, 'k-');
xlim([xmin xmax]);
ylim([ymin ymax]);
set(gca,'tickdir', 'out');
[r, p] = corrcoef(second_filter_index, synergy_index);
title(sprintf('r=%.3f, p=%.5f', r(1,2), p(1,2)))
xlabel('100 I(v_2) / I(v_1)')
ylabel('100 I(v_1,v_2) / I(v_1)+I(v_2)')


% figure;
% % subplot(2,2,2)
% hold on;
% plot(linearity_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
% %plot([0 100], [100 100],'k-')
% xlim([xmin xmax]);
% ylim([ymin ymax]);
% set(gca,'yscale', 'log');
% set(gca,'tickdir', 'out');
% % axis square
% xlabel('Linearity Index')
% ylabel('Log Synergy Index')
% [r, p] = corrcoef(linearity_index, log(synergy_index));
% title(sprintf('r=%.3f, p=%.5f', r(1,2), p(1,2)))
%
% 
% figure;
% % subplot(2,2,3)
% hold on;
% plot(linearity_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
% xlim([xmin xmax]);
% ylim([ymin ymax]);
% set(gca,'xscale', 'log');
% set(gca,'tickdir', 'out');
% % axis square
% xlabel('Log Linearity Index')
% ylabel('Synergy Index')
% [r, p] = corrcoef(log(linearity_index), synergy_index);
% title(sprintf('r=%.3f, p=%.5f', r(1,2), p(1,2)))
% 
% 
% figure;
% % subplot(2,2,4)
% hold on;
% plot(linearity_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
% xlim([xmin xmax]);
% ylim([ymin ymax]);
% set(gca,'xscale', 'log', 'yscale', 'log');
% set(gca,'tickdir', 'out');
% % axis square;
% xlabel('Log Linearity Index');
% ylabel('Log Synergy Index');
% [r, p] = corrcoef(log(linearity_index), log(synergy_index));
% title(sprintf('r=%.3f, p=%.5f', r(1,2), p(1,2)))


% close('all');


%===========================================================
% Plot the linearity index versus the synergy index
%===========================================================
linearity_index = 100 * info1 ./ info_both;
synergy_index = 100 * info_both ./ (info1 + info2);

x = linearity_index;
y = synergy_index;

xmin = 35;
xmax = 105;
ymin = 70;
ymax = 210;

figure;
% xfit = 40:85;
% yfit = 7565 .* xfit .^ (-0.9528) + -20.92;
hold on;
plot(linearity_index, synergy_index, 'ko', 'markersize', 3); %, 'markerfacecolor', 'k');
plot(xfit, yfit, 'k-');
xlim([xmin xmax]);
ylim([ymin ymax]);
set(gca,'tickdir', 'out');
% axis square;
[r, p] = corrcoef(linearity_index, synergy_index);
title(sprintf('Corr between Synergy and Linear Processing, r=%.3f, p=%.5f', r(1,2), p(1,2)))
% xlabel('Linearity Index')
% ylabel('Synergy Index')
xlabel('100 I(v_1) / I(v_1,v_2)')
ylabel('100 I(v_1,v_2) / I(v_1)+I(v_2)')

% close('all');


xmin = 35;
xmax = 95;
ymin = 70;
ymax = 210;

[p,s] = polyfit(log10(x),y,1);
xfit = log10(min(x)):.01:log10(max(x));
yfit = polyval(p,xfit);
xfit = 10.^(xfit);

% xfit = 35:95;
% yfit = 7565 .* xfit .^ (-0.9528) + -20.92;

figure;
hold on;
plot(linearity_index, synergy_index, 'ko', 'markersize', 3);
plot(xfit, yfit, 'k-');
set(gca, 'xtick', 40:10:xmax, 'xticklabel', 40:10:xmax);
set(gca,'xscale', 'log');
xlim([xmin xmax]);
ylim([ymin ymax]);
set(gca,'tickdir', 'out');
% axis square;
[r, p] = corrcoef(log10(linearity_index), synergy_index);
title(sprintf('Corr between Synergy and Linear Processing, r=%.3f, p=%.5f', r(1,2), p(1,2)))
% xlabel('Linearity Index')
% ylabel('Synergy Index')
xlabel('100 I(v_1) / I(v_1,v_2)')
ylabel('100 I(v_1,v_2) / I(v_1)+I(v_2)')


return;






