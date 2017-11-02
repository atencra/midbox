function plot_mid_puretone_params(ptparams)
% plot_mid_puretone_params - plot cf, bw, latency from sta, mid filters
%
% plot_mid_puretone_params(ptparams)
% ---------------------------------------------------------
%
% ptparams : struct array holding parameters for  each filter. It is
% obtained from get_mid_puretone_params.m
%
% caa


sta_latency = [ptparams.sta_latency];
sta_cf = [ptparams.sta_cf];
sta_bw = [ptparams.sta_bw];
sta_q = [ptparams.sta_q];

mid1_latency = [ptparams.mid1_latency];
mid1_cf = [ptparams.mid1_cf];
mid1_bw = [ptparams.mid1_bw];
mid1_q = [ptparams.mid1_q];



close all;

markersize = 3;

% Compare latency

subplot(2,2,1);
minlat = 2.5;
maxlat = 10.5;
hold on;
plot(sta_latency, mid1_latency, 'ko', 'markerfacecolor', 'k', 'markersize', markersize);
plot([minlat maxlat],[minlat maxlat],'k-');
[r,pval] = corrcoef(sta_latency, mid1_latency);
r = r(1,2);
pval = pval(1,2);
text(1.1*minlat, 0.9*maxlat, sprintf('r = %.3f\np = %.3f', r, pval));
box off;
xlim([minlat maxlat]);
ylim([minlat maxlat]);
ticks = round(10*linspace(minlat, maxlat, 5)) ./ 10;
set(gca,'xtick', ticks, 'xticklabel', ticks);
set(gca,'ytick', ticks, 'yticklabel', ticks);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('STA Latency (ms)');
ylabel('MID1 Latency (ms)');

% Compare cf

subplot(2,2,2);
minmin = 0.5;
maxmax = 30;
hold on;
plot(sta_cf, mid1_cf, 'ko', 'markerfacecolor', 'k', 'markersize', markersize);
plot([minmin maxmax],[minmin maxmax],'k-');
[r,pval] = corrcoef(log10(sta_cf), log10(mid1_cf));
r = r(1,2);
pval = pval(1,2);
text(1.1*minmin, 0.6*maxmax, sprintf('r = %.3f\np = %.3f', r, pval));
box off;
xlim([0.5 30]);
ylim([0.5 30]);
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'xtick', [], 'ytick', []);
set(gca,'xminortick', 'off', 'yminortick', 'off');
ticks = round(10*logspace(log10(0.5), log10(30), 5)) ./ 10;
set(gca,'xtick', ticks, 'xticklabel', ticks);
set(gca,'ytick', ticks, 'yticklabel', ticks);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('STA CF (kHz)');
ylabel('MID1 CF (kHz)');

% Compare bandwidth

subplot(2,2,3);
hold on;
minmin = 0.3;
maxmax = 8;
plot(sta_bw, mid1_bw, 'ko', 'markerfacecolor', 'k', 'markersize', markersize);
plot([minmin maxmax],[minmin maxmax],'k-');
[r,pval] = corrcoef(log10(sta_bw), log10(mid1_bw));
r = r(1,2);
pval = pval(1,2);
text(1.1*minmin, 0.6*maxmax, sprintf('r = %.3f\np = %.3f', r, pval));
box off;
xlim([minmin maxmax]);
ylim([minmin maxmax]);
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'xtick', [], 'ytick', []);
set(gca,'xminortick', 'off', 'yminortick', 'off');
ticks = round(10*logspace(log10(minmin), log10(maxmax), 5)) ./ 10;
set(gca,'xtick', ticks, 'xticklabel', ticks);
set(gca,'ytick', ticks, 'yticklabel', ticks);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('STA BW (kHz)');
ylabel('MID1 BW (kHz)');


% Compare Q

subplot(2,2,4);
hold on;
minmin = 0.4;
maxmax = 7;
plot(sta_q, mid1_q, 'ko', 'markerfacecolor', 'k', 'markersize', markersize);
plot([minmin maxmax],[minmin maxmax],'k-');
[r,pval] = corrcoef(log10(sta_q), log10(mid1_q));
r = r(1,2);
pval = pval(1,2);
text(1.1*minmin, 0.6*maxmax, sprintf('r = %.3f\np = %.3f', r, pval));
box off;
xlim([minmin maxmax]);
ylim([minmin maxmax]);
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'xtick', [], 'ytick', []);
set(gca,'xminortick', 'off', 'yminortick', 'off');
ticks = round(10*logspace(log10(minmin), log10(maxmax), 5)) ./ 10;
set(gca,'xtick', ticks, 'xticklabel', ticks);
set(gca,'ytick', ticks, 'yticklabel', ticks);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('STA Q');
ylabel('MID1 Q');

print_mfilename(mfilename);
