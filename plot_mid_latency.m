function plot_mid_latency(mid, saveit)
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

latency1 = [];
latency2 = [];

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



   % Plot the MIDs
   % ========================================================

	t = [0 25 50 75 100];
	indext = [1 6 11 16 20];

	if ( fbins == 30 )
		f = [1.25 2.5 5 10 20 38];
		indexf = [1 7 13 19 25 30];
	elseif ( fbins == 25 )
		f = [1.25 2.5 5 10 20];
		indexf = [1 7 13 19 25];
	end

   time = 5* (0:(size(v1,2)-1));

%    % Plot the first MID
%    %------------------------------
%    subplot(2,2,1);
%    minmin = min(min(v1));
%    maxmax = max(max(v1));
%    boundary = max([abs(minmin) abs(maxmax)]);
%    imagesc(fliplr(v1));
%    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    axis('xy');
% 	set(gca, 'xtick', indext, 'xticklabel', t);
% 	set(gca, 'ytick', indexf, 'yticklabel', f);
% 	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
% 	xlabel('Time (ms)');
% 	ylabel('Frequency (kHz)');
%    title('MID_1');
% 
%    % Plot the second MID
%    %------------------------------
%    subplot(2,2,2);
%    minmin = min(min(v2));
%    maxmax = max(max(v2));
%    boundary = max([abs(minmin) abs(maxmax)]);
%    imagesc(fliplr(v2));
%    axis('xy');
% 	set(gca, 'xtick', indext, 'xticklabel', t);
% 	set(gca, 'ytick', indexf, 'yticklabel', f);
% 	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
% 	xlabel('Time (ms)');
% 	title('MID_2');
% 
   trf1 = sum(fliplr(v1),1);
   trf2 = sum(fliplr(v2),1);
% 
%    subplot(2,2,3);
%    plot(time, trf1);
% 	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%    title('TRF1');
% 
% 
%    subplot(2,2,4);
%    plot(time, trf2);	
%    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%    title('TRF2');


   index1 = find( trf1 == max(trf1), 1 );
   delay1 = time( min( index1 ) );

   index2 = find( trf2 == max(trf2), 1 );
   delay2 = time( min( index2 ) );


   latency1 = [latency1 delay1];
   latency2 = [latency2 delay2];

% pause;
% 
% close all;

end % (for i)

% return;


index = find( latency1 ~= 0 & latency2 ~= 0 );

latency1 = latency1(index);
latency2 = latency2(index);

[h,p] = kstest2(latency1, latency2)

median(latency1)
median(latency2)
mean(latency1)
mean(latency2)

%------------------------------------------------
% Best Temporal Modulation Frequency
%------------------------------------------------

% Make the scatter plot of btmf values
% ==================================================
[r,pval] = corrcoef(latency1, latency2);
r = r(1,2);
pval = pval(1,2);

[beta] = polyfit(latency1, latency2, 1);
xfit = 0:.01:100;
yfit = polyval(beta, xfit);


hscatter = subplot(2,2,4);
hold on;
hp = plot(latency1, latency2, 'ko');
xlim([0 100]);
ylim([0 100]);
plot(xfit, yfit, 'k-');
set(hp,'markeredgecolor', 'k');
set(hp,'markerfacecolor', 'k');
set(hp, 'markersize', 2);
set(gca, 'xtick', 0:20:100, 'xticklabel', 0:20:100);
set(gca, 'ytick', 0:20:100, 'yticklabel', 0:20:100);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlabel('Latency MID 1 (ms)');
ylabel('Latency MID 2 (ms)');
title(sprintf('r = %.2f, p = %.4f', r, pval));



% Make the marginal for the x-axis
% ==================================================
hx = subplot(2,2,2);
edges = linspace(0,100,16);
n1 = histc(latency1, edges);
hb = bar(edges, n1, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'xtick', 0:20:100, 'xticklabel', '');
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlim([0 100]);
ylabel('#Neurons');


% Make the marginal for the y-axis
% ==================================================
hy = subplot(2,2,3);
d = diff(edges);
d = d(1);
yctrs = d/2:d:max(edges);
[n2,cy] = hist(latency2,yctrs);
hb = barh(cy, -n2, 1); 

set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca,'ytick', 0:20:100, 'yticklabel', '');
xtick = get(gca,'xtick');
set(gca,'xticklabel', abs(xtick), 'xtick', xtick);
ylim([0 100]);
xlabel('#Neurons');


% Make scatter plot bigger, histograms smaller
% ==================================================
set(hscatter,'Position',[0.35 0.35 0.55 0.55],'tag','scatter');
set(hx,'Position',[.35 .075 .55 .15],'tag','xhist');
set(hy,'Position',[.075 .35 .15 .55],'tag','yhist');


return;














