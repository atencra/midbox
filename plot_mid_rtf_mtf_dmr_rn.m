function plot_mid_rtf_mtf_dmr_rn(dmrmid, rnmid)
%plot_mid_rtf_mtf_dmr_rn - Plot modulation parameters for DMR/RN MIDs
%
% plot_mid_rtf_mtf(mid, saveit)
%
% mid : struct array holding MID data. It is saved in a file called:
%
%		mid_fio_info_database.mat                              
%
%   caa 10/28/08

% if ( nargin == 1 )
%    saveit = 0;
% end
% 
% if ( nargin == 2 )
%    if ( saveit ~= 1 )
%       saveit = 0;
%    end
% end

ntot = 1;
fignum = 1;
figure;

numplots = ceil(length(dmrmid)/4);

close all;

rn_tbmf1 = zeros(1,length(dmrmid));
rn_tbmf2 = zeros(1,length(dmrmid));
rn_xbmf1 = zeros(1,length(dmrmid));
rn_xbmf2 = zeros(1,length(dmrmid));

dmr_tbmf1 = zeros(1,length(dmrmid));
dmr_tbmf2 = zeros(1,length(dmrmid));
dmr_xbmf1 = zeros(1,length(dmrmid));
dmr_xbmf2 = zeros(1,length(dmrmid));

for i = 1:length(dmrmid)


   exp = dmrmid(i).exp;
   site = dmrmid(i).site;

   dmrunit = dmrmid(i).unit;
   rnunit = rnmid(i).unit;

   [dmrdata] = get_mid_filters_rtf_mtf_for_plots( dmrmid(i) );

   t = dmrdata.t;
   indext = dmrdata.indext;
   f = dmrdata.f;
   indexf = dmrdata.indexf;
   v1 = dmrdata.v1;
   v2 = dmrdata.v2;
   indextmf = dmrdata.indextmf;
   indexxmf = dmrdata.indexxmf;
   tmf1s = dmrdata.tmf1s;
   xmf1s = dmrdata.xmf1s;
   rtf1s1 = dmrdata.rtf1s1;
   rtf1s2 = dmrdata.rtf1s2;
   tmtf1 = dmrdata.tmtf1;
   tmtf2 = dmrdata.tmtf2;
   xmtf1 = dmrdata.xmtf1;
   xmtf2 = dmrdata.xmtf2;
   tbmf1 = dmrdata.tbmf1;
   tbmf2 = dmrdata.tbmf2;
   xbmf1 = dmrdata.xbmf1;
   xbmf1 = dmrdata.xbmf2;

   dmr_tbmf1(i) = dmrdata.tbmf1;
   dmr_tbmf2(i) = dmrdata.tbmf2;
   dmr_xbmf1(i) = dmrdata.xbmf1;
   dmr_xbmf2(i) = dmrdata.xbmf2;


%    [t, indext, f, indexf, v1, v2, indextmf, indexxmf, ...
%    tmf1s, xmf1s, rtf1s1, rtf1s2, tmtf1, tmtf2, xmtf1, xmtf2 ] = ...
%       get_mid_filters_rtf_mtf_for_plots( dmrmid(i) );


   figure;

   % Plot the MIDs
   % ========================================================

% 	t = [0 25 50 75 100];
% 	indext = [1 6 11 16 20];
% 
% 	if ( fbins == 30 )
% 		f = [1.25 2.5 5 10 20 38];
% 		indexf = [1 7 13 19 25 30];
% 	elseif ( fbins == 25 )
% 		f = [1.25 2.5 5 10 20];
% 		indexf = [1 7 13 19 25];
% 	end

   % Plot the first MID
   %------------------------------
   subplot(4,4,1);
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
   title('DMR MID_1');

   % Plot the second MID
   %------------------------------
   subplot(4,4,2);
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
	title('DMR MID_2');


   % Plot the RTFs
   % ========================================================

   subplot(4,4,5);
   imagesc(rtf1s1);
	set(gca, 'xtick', indextmf, 'xticklabel', tmf1s(indextmf));
	set(gca, 'ytick', indexxmf, 'yticklabel', round(xmf1s(indexxmf)));
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   axis xy;
	ylabel(sprintf('Spectral Modulation\nFreq (cyc/oct)'));


   subplot(4,4,6);
   imagesc(rtf1s2);
	set(gca, 'xtick', indextmf, 'xticklabel', tmf1s(indextmf));
	set(gca, 'ytick', indexxmf, 'yticklabel', round(xmf1s(indexxmf)));
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   axis xy;


   % Plot the tMTFs
   % ========================================================


   subplot(4,4,9);
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


   subplot(4,4,10);
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


   subplot(4,4,13);
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


   subplot(4,4,14);
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


   [rndata] = get_mid_filters_rtf_mtf_for_plots( rnmid(i) );

   t = rndata.t;
   indext = rndata.indext;
   f = rndata.f;
   indexf = rndata.indexf;
   v1 = rndata.v1;
   v2 = rndata.v2;
   indextmf = rndata.indextmf;
   indexxmf = rndata.indexxmf;
   tmf1s = rndata.tmf1s;
   xmf1s = rndata.xmf1s;
   rtf1s1 = rndata.rtf1s1;
   rtf1s2 = rndata.rtf1s2;
   tmtf1 = rndata.tmtf1;
   tmtf2 = rndata.tmtf2;
   xmtf1 = rndata.xmtf1;
   xmtf2 = rndata.xmtf2;
   tbmf1 = rndata.tbmf1;
   tbmf2 = rndata.tbmf2;
   xbmf1 = rndata.xbmf1;
   xbmf1 = rndata.xbmf2;

   rn_tbmf1(i) = rndata.tbmf1;
   rn_tbmf2(i) = rndata.tbmf2;
   rn_xbmf1(i) = rndata.xbmf1;
   rn_xbmf2(i) = rndata.xbmf2;


   % Plot the first MID
   %------------------------------
   subplot(4,4,3);
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
% 	ylabel('Frequency (kHz)');
   title('RN MID_1');

   % Plot the second MID
   %------------------------------
   subplot(4,4,4);
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
	title('RN MID_2');


   % Plot the RTFs
   % ========================================================

   subplot(4,4,7);
   imagesc(rtf1s1);
	set(gca, 'xtick', indextmf, 'xticklabel', tmf1s(indextmf));
	set(gca, 'ytick', indexxmf, 'yticklabel', round(xmf1s(indexxmf)));
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   axis xy;
% 	ylabel(sprintf('Spectral Modulation\nFreq (cyc/oct)'));


   subplot(4,4,8);
   imagesc(rtf1s2);
	set(gca, 'xtick', indextmf, 'xticklabel', tmf1s(indextmf));
	set(gca, 'ytick', indexxmf, 'yticklabel', round(xmf1s(indexxmf)));
	set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
   axis xy;


   % Plot the tMTFs
   % ========================================================


   subplot(4,4,11);
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
% 	ylabel(sprintf('Normalized\nAmplitude'));


   subplot(4,4,12);
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


   subplot(4,4,15);
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
% 	ylabel(sprintf('Normalized\nAmplitude'));


   subplot(4,4,16);
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

   suptitle(sprintf('#%.0f: %s site %.0f\nDMR unit %.0f, RN unit %.0f', i, exp, site, dmrunit, rnunit));

   print_mfilename(mfilename);
   orient portrait;
   set(gcf,'position', [242   147   665   678]);

end % (for i)



figure;

for i = 1:length(dmr_tbmf1)
   
   % Plot best temporal/spectral modulation frequency for filter 1 and 2
   % color code dmr and rn stimuli
   subplot(2,1,1);
   hold on;
   x = [dmr_tbmf1(i) rn_tbmf1(i)];
   y = [dmr_xbmf1(i) rn_xbmf1(i)];
   plot(dmr_tbmf1(i), dmr_xbmf1(i), 'ko', 'markerfacecolor', 'k');
   plot(rn_tbmf1(i), rn_xbmf1(i), 'ro', 'markerfacecolor', 'r');
   plot(x,y,'k-');
   tickpref;
   xlim([0 40]);
   ylim([0 4]);
   xlabel('Temporal Modulation (cyc/s)');
   ylabel('Spectral Modulation (cyc/oct)');
   title('MID1');

   subplot(2,1,2);
   hold on;
   plot(dmr_tbmf2(i), dmr_xbmf2(i), 'ko', 'markerfacecolor', 'k');
   plot(rn_tbmf2(i), rn_xbmf2(i), 'ro', 'markerfacecolor', 'r');
   x = [dmr_tbmf2(i) rn_tbmf2(i)];
   y = [dmr_xbmf2(i) rn_xbmf2(i)];
   plot(x,y,'k-');
   tickpref;
   xlim([0 40]);
   ylim([0 4]);
   xlabel('Temporal Modulation (cyc/s)');
   ylabel('Spectral Modulation (cyc/oct)');
   title('MID2');

end % (for i)

legend('DMR', 'RN');








figure;

for i = 1:length(dmr_tbmf1)
   
   % Plot best temporal/spectral modulation frequency for filter 1 and 2
   % color code dmr and rn stimuli
   subplot(2,2,1);
   hold on;
   plot([1, 2],[dmr_tbmf1(i) rn_tbmf1(i)], 'k-');
   tickpref;
   xlim([0.5 2.5]);
   ylim([0 25]);
   set(gca,'xtick', [1 2], 'xticklabel', {'DMR', 'RN'});
   ylabel('Temporal Modulation (cyc/s)');
   title('MID1');

   subplot(2,2,2);
   hold on;
   plot([1, 2],[dmr_tbmf2(i) rn_tbmf2(i)], 'k-');
   tickpref;
   xlim([0.5 2.5]);
   ylim([0 25]);
   set(gca,'xtick', [1 2], 'xticklabel', {'DMR', 'RN'});
   ylabel('Temporal Modulation (cyc/s)');
   title('MID2');

   subplot(2,2,3);
   hold on;
   plot([1, 2],[dmr_xbmf1(i) rn_xbmf1(i)], 'k-');
   tickpref;
   xlim([0.5 2.5]);
   ylim([0 1.5]);
   set(gca,'xtick', [1 2], 'xticklabel', {'DMR', 'RN'});
   ylabel('Spectral Modulation (cyc/s)');
   title('MID1');

   subplot(2,2,4);
   hold on;
   plot([1, 2],[dmr_xbmf2(i) rn_xbmf2(i)], 'k-');
   tickpref;
   xlim([0.5 2.5]);
   ylim([0 1.5]);
   set(gca,'xtick', [1 2], 'xticklabel', {'DMR', 'RN'});
   ylabel('Spectral Modulation (cyc/s)');
   title('MID2');

end % (for i)



p = signrank(dmr_tbmf1, rn_tbmf1);
fprintf('bTMF: DMR vs RN MID1: p = %.4f\n', p);

p = signrank(dmr_tbmf2, rn_tbmf2);
fprintf('bTMF: DMR vs RN MID2: p = %.4f\n', p);

p = signrank(dmr_xbmf1, rn_xbmf1);
fprintf('bSMF: DMR vs RN MID1: p = %.4f\n', p);

p = signrank(dmr_xbmf2, rn_xbmf2);
fprintf('bSMF: DMR vs RN MID2: p = %.4f\n', p);


return;





function [data] = get_mid_filters_rtf_mtf_for_plots(mid)

% [t, indext, f, indexf, v1, v2, indextmf, indexxmf, ...
% tmf1s, xmf1s, rtf1s1, rtf1s2, tmtf1, tmtf2, xmtf1, xmtf2 ] = get_mid_filters_rtf_mtf_for_plots(mid)

max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];


% Get number of time and frequency bins used to make filters
% Note: the time bins are 5 ms, and the freq bins are 1/6 octave
tbins = mid.tbins;
fbins = mid.fbins;

t = [0 25 50 75 100];
indext = [1 6 11 16 20];

if ( fbins == 30 )
   f = [1.25 2.5 5 10 20 38];
   indexf = [1 7 13 19 25 30];
elseif ( fbins == 25 )
   f = [1.25 2.5 5 10 20];
   indexf = [1 7 13 19 25];
end


% Get mean firing rate (probability) over the entire stimulus
pspike = mid.rpx1pxpxt_sta.mean_firing;

% Get STA
sta = mid.rpsta.filter;
sta = reshape(sta, fbins, tbins);
%    sta(abs(sta)<0.75) = 0;
x_sta = mid.rpx1pxpxt_sta.x;
fx_sta = mid.rpx1pxpxt_sta.ior_mean;
pspx_sta = fx_sta;
px_sta = mid.rpx1pxpxt_sta.px_mean;
pxsp_sta = fx_sta ./ pspike .* px_sta;


% Get MID1
v1 = mid.rpdtest2_v1.filter;
v1 = reshape(v1, fbins, tbins);
%    v1(abs(v1)<0.75) = 0;
x_v1 = mid.rpdx1x2px_pxt_2.x1;
fx_v1 = mid.rpdx1x2px_pxt_2.ior1_mean;
pspx_v1 = fx_v1;
px_v1 = mid.rpdx1x2px_pxt_2.px1_mean;
pxsp_v1 = fx_v1 ./ pspike .* px_v1;

% % Get rid of strange boundary effect in MID1
% v1(:,end) = 0;
% v1(:,1) = 0;
% v1(end,:) = 0;
% v1(end-1,:) = 0;


% Get MID2
v2 = mid.rpdtest2_v2.filter;
v2 = reshape(v2, fbins, tbins);
%    v2(abs(v2)<0.75) = 0;
x_v2 = mid.rpdx1x2px_pxt_2.x2;
fx_v2 = mid.rpdx1x2px_pxt_2.ior2_mean;
pspx_v2 = fx_v2;
px_v2 = mid.rpdx1x2px_pxt_2.px2_mean;
pxsp_v2 = fx_v2 ./ pspike .* px_v2;

% % Get rid of strange boundary effect in MID2
% v2(:,end) = 0;
% v2(:,1) = 0;
% v2(end,:) = 0;
% v2(end-1,:) = 0;


% Get the 2D MID nonlinearity
x12 = mid.rpdx1x2px_pxt_2.x12;
y12 = mid.rpdx1x2px_pxt_2.y12;
fx_v1v2 = mid.rpdx1x2px_pxt_2.ior12_mean;

% Keep track of min/max values of 1D nonlinearities for axis bounds
max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];

v1 = smoothmat(v1);
v2 = smoothmat(v2);

[tmf, xmf, rtf1] = mid_rtf(v1);
[tmf, xmf, rtf2] = mid_rtf(v2);

[tmf1s, tmtf1, xmf1s, xmtf1, rtf1s1] = mid_mtf(tmf, xmf, rtf1);
[tmf1s, tmtf2, xmf1s, xmtf2, rtf1s2] = mid_mtf(tmf, xmf, rtf2);

[tbmf1, tbw6db1, tbw3db1, twc3db1] = mtf_bmf_bw_wc(tmf1s, tmtf1);
[tbmf2, tbw6db2, tbw3db2, twc3db2] = mtf_bmf_bw_wc(tmf1s, tmtf2);

[xbmf1, xbw6db1, xbw3db1, xwc3db1] = mtf_bmf_bw_wc(xmf1s, xmtf1);
[xbmf2, xbw6db2, xbw3db2, xwc3db2] = mtf_bmf_bw_wc(xmf1s, xmtf2);

xmf1s = round(10*xmf1s) / 10; % keep only 1 decimal place

[temp, temp, indextmf] = intersect([0 10 20 30 40], tmf1s);
[temp, temp, indexxmf] = intersect([0 1.0 2.0 2.9], xmf1s);


data.t = t;
data.indext = indext;
data.f = f;
data.indexf = indexf;
data.v1 = v1;
data.v2 = v2;
data.indextmf = indextmf;
data.indexxmf = indexxmf;
data.tmf1s = tmf1s;
data.xmf1s = xmf1s;
data.rtf1s1 = rtf1s1;
data.rtf1s2 = rtf1s2;
data.tmtf1 = tmtf1;
data.tmtf2 = tmtf2;
data.xmtf1 = xmtf1;
data.xmtf2 = xmtf2;
data.tbmf1 = tbmf1;
data.tbmf2 = tbmf2;
data.xbmf1 = xbmf1;
data.xbmf2 = xbmf2;

return;












