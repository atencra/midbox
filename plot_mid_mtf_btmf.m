function plot_mid_mtf_btmf(mid, saveit)
%plot_mid_rtf_mtf - Plot modulation processing parameters for MIDs
%
% plot_mid_rtf_mtf(mid, saveit)
%
% mid : struct array holding MID data. It is saved in a file called:
%
%		mid_fio_info_database.mat                              
%
%   caa 10/28/08

if ( nargin == 1 )
   saveit = 0;
end

if ( nargin == 2 )
   if ( saveit ~= 1 )
      saveit = 0;
   end
end

ntot = 1;
fignum = 1;
figure;

max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];

numplots = ceil(length(mid)/4);

gaussian = fspecial('gaussian', [3 3], 1);

close all;

btmf1_pop = [];
btmf2_pop = [];

bxmf1_pop = [];
bxmf2_pop = [];

for i = 1:length(mid)

   % Get number of time and frequency bins used to make filters
   % Note: the time bins are 5 ms, and the freq bins are 1/6 octave
   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   % Get mean firing rate (probability) over the entire stimulus
   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   % Get STA
   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
   sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   fx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   pspx_sta = fx_sta;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxsp_sta = fx_sta ./ pspike .* px_sta;


   % Get MID1
   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   pspx_v1 = fx_v1;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxsp_v1 = fx_v1 ./ pspike .* px_v1;

   % Get rid of strange boundary effect in MID1
   v1(:,end) = 0;
   v1(:,1) = 0;
   v1(end,:) = 0;
   v1(end-1,:) = 0;


   % Get MID2
   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

   % Get rid of strange boundary effect in MID2
   v2(:,end) = 0;
   v2(:,1) = 0;
   v2(end,:) = 0;
   v2(end-1,:) = 0;


   % Get the 2D MID nonlinearity
   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   % Keep track of min/max values of 1D nonlinearities for axis bounds
   max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
   max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
   max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];

	% apply some minimal smoothing to the filters to make them look pretty
	v1 = imfilter(v1,gaussian);
	v2 = imfilter(v2,gaussian);

	% Process all the modulation functions/parameters
   [tmf, xmf, rtf1] = mid_rtf(v1);
   [tmf, xmf, rtf2] = mid_rtf(v2);

   [tmf1s, tmtf1, xmf1s, xmtf1, rtf1s1] = mid_mtf(tmf, xmf, rtf1);
   [tmf1s, tmtf2, xmf1s, xmtf2, rtf1s2] = mid_mtf(tmf, xmf, rtf2);

	[tbmf1, tbw6db1, tbw3db1, twc3db1] = mtf_bmf_bw_wc(tmf1s, tmtf1);
	[tbmf2, tbw6db2, tbw3db2, twc3db2] = mtf_bmf_bw_wc(tmf1s, tmtf2);

	[xbmf1, xbw6db1, xbw3db1, xwc3db1] = mtf_bmf_bw_wc(xmf1s, xmtf1);
	[xbmf2, xbw6db2, xbw3db2, xwc3db2] = mtf_bmf_bw_wc(xmf1s, xmtf2);

	btmf1_pop = [btmf1_pop tbmf1];
	btmf2_pop = [btmf2_pop tbmf2];

	bxmf1_pop = [bxmf1_pop xbmf1];
	bxmf2_pop = [bxmf2_pop xbmf2];


	% Get values for nice mtf tick marks
	xmf1s = round(10*xmf1s) / 10; % keep only 1 decimal place
	[temp, temp, indextmf] = intersect([0 10 20 30 40], tmf1s);
	[temp, temp, indexxmf] = intersect([0 1.0 2.0 2.9], xmf1s);

end % (for i)



%------------------------------------------------
% Best Temporal Modulation Frequency
%------------------------------------------------

% Make the scatter plot of btmf values
% ==================================================
[r,pval] = corrcoef(btmf1_pop, btmf2_pop);
r = r(1,2);
pval = pval(1,2);

[beta] = polyfit(btmf1_pop, btmf2_pop, 1);
xfit = 0:.01:25;
yfit = polyval(beta, xfit);


hscatter = subplot(2,2,4);
hold on;
hp = plot(btmf1_pop, btmf2_pop, 'ko');
xlim([0 25]);
ylim([0 25]);
plot(xfit, yfit, 'k-');
set(hp,'markeredgecolor', 'k');
set(hp,'markerfacecolor', 'k');
set(hp, 'markersize', 2);
set(gca, 'xtick', 0:5:25, 'xticklabel', 0:5:25);
set(gca, 'ytick', 0:5:25, 'yticklabel', 0:5:25);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlabel('MID 1 Best TMF (cyc/s)');
ylabel('MID 2 Best TMF (cyc/s)');
title(sprintf('r = %.2f, p = %.4f', r, pval));



% Make the marginal for the x-axis
% ==================================================
hx = subplot(2,2,2);
edges = linspace(0,25,16);
n1 = histc(btmf1_pop, edges);
hb = bar(edges, n1, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'xtick', 0:5:25, 'xticklabel', '');
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlim([0 25]);
ylabel('#Neurons');


% Make the marginal for the y-axis
% ==================================================
hy = subplot(2,2,3);
d = diff(edges);
d = d(1);
yctrs = d/2:d:max(edges);
[n2,cy] = hist(btmf2_pop,yctrs);
hb = barh(cy, -n2, 1); 

set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca,'ytick', 0:5:25, 'yticklabel', '');
xtick = get(gca,'xtick');
set(gca,'xticklabel', abs(xtick), 'xtick', xtick);
ylim([0 25]);
xlabel('#Neurons');


% Make scatter plot bigger, histograms smaller
% ==================================================
set(hscatter,'Position',[0.35 0.35 0.55 0.55],'tag','scatter');
set(hx,'Position',[.35 .075 .55 .15],'tag','xhist');
set(hy,'Position',[.075 .35 .15 .55],'tag','yhist');




%------------------------------------------------
% Best Spectral Modulation Frequency
%------------------------------------------------

% Make the scatter plot of bsmf values
% ==================================================
[r,pval] = corrcoef(bxmf1_pop, bxmf2_pop);
r = r(1,2);
pval = pval(1,2);

[beta] = polyfit(bxmf1_pop, bxmf2_pop, 1);
xfit = 0:.01:1.5;
yfit = polyval(beta, xfit);

xlim1 = 0;
xlim2 = 1.2;
dx = 0.4;

figure;

hscatter = subplot(2,2,4);
hold on;
hp = plot(bxmf1_pop, bxmf2_pop, 'ko');
xlim([xlim1 xlim2]);
ylim([xlim1 xlim2]);
plot(xfit, yfit, 'k-');
set(hp,'markeredgecolor', 'k');
set(hp,'markerfacecolor', 'k');
set(hp, 'markersize', 2);
set(gca, 'xtick', 0:dx:xlim2, 'xticklabel', 0:dx:xlim2);
set(gca, 'ytick', 0:dx:xlim2, 'yticklabel', 0:dx:xlim2);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlabel('MID 1 Best SMF (cyc/oct)');
ylabel('MID 2 Best SMF (cyc/oct)');
title(sprintf('r = %.2f, p = %.4f', r, pval));



% Make the marginal for the x-axis
% ==================================================
hx = subplot(2,2,2);
edges = linspace(0,xlim2,16);
n1 = histc(bxmf1_pop, edges);
hb = bar(edges, n1, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'xtick', 0:dx:xlim2, 'xticklabel', '');
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlim([xlim1 xlim2]);
ylabel('#Neurons');


% Make the marginal for the y-axis
% ==================================================
hy = subplot(2,2,3);
d = diff(edges);
d = d(1);
yctrs = d/2:d:max(edges);
[n2,cy] = hist(bxmf2_pop,yctrs);
hb = barh(cy, -n2, 1); 

set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca,'ytick', 0:dx:xlim2, 'yticklabel', '');
xtick = get(gca,'xtick');
set(gca,'xticklabel', abs(xtick), 'xtick', xtick);
ylim([xlim1 xlim2]);
xlabel('#Neurons');


% Make scatter plot bigger, histograms smaller
% ==================================================
set(hscatter,'Position',[0.35 0.35 0.55 0.55],'tag','scatter');
set(hx,'Position',[.35 .075 .55 .15],'tag','xhist');
set(hy,'Position',[.075 .35 .15 .55],'tag','yhist');



return;














