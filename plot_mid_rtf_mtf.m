function plot_mid_rtf_mtf(mid, saveit)
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

	v1 = imfilter(v1,gaussian);
	v2 = imfilter(v2,gaussian);

   [tmf, xmf, rtf1] = mid_rtf(v1);
   [tmf, xmf, rtf2] = mid_rtf(v2);

   [tmf1s, tmtf1, xmf1s, xmtf1, rtf1s1] = mid_mtf(tmf, xmf, rtf1);
   [tmf1s, tmtf2, xmf1s, xmtf2, rtf1s2] = mid_mtf(tmf, xmf, rtf2);

	[tbmf1, tbw6db1, tbw3db1, twc3db1] = mtf_bmf_bw_wc(tmf1s, tmtf1);
	[tbmf2, tbw6db2, tbw3db2, twc3db2] = mtf_bmf_bw_wc(tmf1s, tmtf2);



   figure;

	xmf1s = round(10*xmf1s) / 10; % keep only 1 decimal place

	[temp, temp, indextmf] = intersect([0 10 20 30 40], tmf1s);
	[temp, temp, indexxmf] = intersect([0 1.0 2.0 2.9], xmf1s);



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

   % Plot the first MID
   %------------------------------
   subplot(4,2,1);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(fliplr(v1));
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   axis('xy');
	set(gca, 'xtick', indext, 'xticklabel', t);
	set(gca, 'ytick', indexf, 'yticklabel', f);
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
	xlabel('Time (ms)');
	ylabel('Frequency (kHz)');
   title('MID_1');

   % Plot the second MID
   %------------------------------
   subplot(4,2,2);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(fliplr(v2));
   axis('xy');
	set(gca, 'xtick', indext, 'xticklabel', t);
	set(gca, 'ytick', indexf, 'yticklabel', f);
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
	xlabel('Time (ms)');
	title('MID_2');


   % Plot the RTFs
   % ========================================================

   subplot(4,2,3);
   imagesc(rtf1s1);
	set(gca, 'xtick', indextmf, 'xticklabel', tmf1s(indextmf));
	set(gca, 'ytick', indexxmf, 'yticklabel', round(xmf1s(indexxmf)));
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   axis xy;
	ylabel(sprintf('Spectral Modulation\nFreq (cyc/oct)'));


   subplot(4,2,4);
   imagesc(rtf1s2);
	set(gca, 'xtick', indextmf, 'xticklabel', tmf1s(indextmf));
	set(gca, 'ytick', indexxmf, 'yticklabel', round(xmf1s(indexxmf)));
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   axis xy;


   % Plot the tMTFs
   % ========================================================


   subplot(4,2,5);
   hold on;
   hp = plot(tmf1s, tmtf1, 'ko-');
	ylim([0 1]);
	set(hp,'markeredgecolor', 'k');
	set(hp,'markerfacecolor', 'k');
	set(hp, 'markersize', 2);
	set(gca, 'xtick', 0:10:40, 'xticklabel', 0:10:40);
	set(gca, 'ytick', 0:0.25:1, 'yticklabel', 0:0.25:1);
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   xlabel('Temporal Modulation Freq (cyc/s)');
	ylabel(sprintf('Normalized\nAmplitude'));


   subplot(4,2,6);
   hold on;
   hp = plot(tmf1s, tmtf2, 'ko-');
	ylim([0 1]);
	set(hp,'markeredgecolor', 'k');
	set(hp,'markerfacecolor', 'k');
	set(hp, 'markersize', 2);
	set(gca, 'xtick', 0:10:40, 'xticklabel', 0:10:40);
	set(gca, 'ytick', 0:0.25:1, 'yticklabel', 0:0.25:1);
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   xlabel('Temporal Modulation Freq (cyc/s)');


   % Plot the tMTFs
   % ========================================================


   subplot(4,2,7);
   hold on;
   hp = plot(xmf1s, xmtf1, 'ko-');
	ylim([0 1]);
   xlim([0 1.05*max(xmf1s)]);
	set(hp,'markeredgecolor', 'k');
	set(hp,'markerfacecolor', 'k');
	set(hp, 'markersize', 2);
	set(gca, 'xtick', 0:1:4, 'xticklabel', 0:1:4);
	set(gca, 'ytick', 0:0.25:1, 'yticklabel', 0:0.25:1);
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   xlabel('SM Frequency (cyc/s)');
	ylabel(sprintf('Normalized\nAmplitude'));


   subplot(4,2,8);
   hold on;
   hp = plot(xmf1s, xmtf2, 'ko-');
	ylim([0 1]);
   xlim([0 1.05*max(xmf1s)]);
	set(hp,'markeredgecolor', 'k');
	set(hp,'markerfacecolor', 'k');
	set(hp, 'markersize', 2);
	set(gca, 'xtick', 0:1:4, 'xticklabel', 0:1:4);
	set(gca, 'ytick', 0:0.25:1, 'yticklabel', 0:0.25:1);
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   xlabel('SM Frequency (cyc/s)');


   print_mfilename(mfilename);

   set(gcf,'position', [360   520   290   402]);

end % (for i)

return;















