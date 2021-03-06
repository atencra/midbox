function plot_rat_sta_mid1_mid2_fx(mid)
% plot_sta_mid1_mid2_fx(mid) - 
%
% For illustrations in adobe the following is useful:
%
% exportfig(gcf,'sta_v1_v2_fx.eps', 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
%
% caa 6/01/06


ntot = 1;
% len = length(strf);


for i = 1:length(mid)

   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
   sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   pspkx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxspk_sta = pspkx_sta ./ pspike .* px_sta;

   pspkx_sta = pspkx_sta ./ 0.005;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   pspkx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxspk_v1 = pspkx_v1 ./ pspike .* px_v1;
% [x_v1(:) pspkx_v1(:)]

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   pspkx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxspk_v2 = pspkx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   pspkx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;


	gaussian = fspecial('gaussian', [3 3], 1);
	sta = imfilter(sta,gaussian);
%    v1(:,end) = 0.9 * (rand(length(v1(:,end)),1)-0.5);
	v1 = imfilter(v1,gaussian);
   v2 = imfilter(v2,gaussian);



   figure;

   % Plot STA and its nonlinearity
   %------------------------------
   subplot(1,4,1);

   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(sta);
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   xlabel('Time (ms)');
   ylabel('Frequency (kHz)');
%    ylabel(sprintf('%.0f : %.0f - %.0f', i, mid(i).location, mid(i).unit));
   title('STA', 'fontsize', 8);


   subplot(1,4,2);
   hold on;
   plot(x_sta, pspkx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7 -3.5 0 3.5 7], 'xticklabel', [-7 -3.5 0 3.5 7]);
%    set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   mx = str2num(sprintf('%.0f', 1.025*max(pspkx_sta)));
   mx2 = str2num(sprintf('%.0f', mx/2));
   ylim([0 mx]);
   set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca,'fontsize', 8);
   ylabel('Firing Rate (sp/s)');
   title('STA Nonlinearity', 'fontsize', 8);


   % Plot V1 and its nonlinearity
   %------------------------------
   subplot(1,4,3);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v1);
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   title('MID_1', 'fontsize', 8);


%    subplot(1,6,4);
%    hold on;
%    plot(x_v1, pspkx_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
%    plot([min(x_v1) max(x_v1)], [pspike pspike], 'k:');
%    xlim([-7.5 7.5]);
%    set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
%    mx = str2num(sprintf('%.3f', 1.025*max(pspkx_v1)));
%    mx2 = str2num(sprintf('%.3f', mx/2));
%    ylim([0 mx]);
%    set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);
% 
%    set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    set(gca,'fontsize', 8);
%    title('MID_1 P(spk|proj)', 'fontsize', 8);

   % Plot V2 and its nonlinearity
   %------------------------------
   subplot(1,4,4);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v2);
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   title('MID_2', 'fontsize', 8);


%    subplot(1,6,6);
%    hold on;
%    plot(x_v2, pspkx_v2, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
%    plot([min(x_v2) max(x_v2)], [pspike pspike], 'k:');
%    xlim([-6 6]);
%    set(gca,'xtick', [-5 -2.5 0 2.5 5], 'xticklabel', [-5 -2.5 0 2.5 5]);
%    mx = str2num(sprintf('%.3f', 1.025*max(pspkx_v2)));
%    mx2 = str2num(sprintf('%.3f', mx/2));
%    ylim([0 mx]);
%    set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);
%    set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    set(gca,'fontsize', 8);
%    title('MID_2 P(spk|proj)', 'fontsize', 8);


%    % Plot P(spk|x1,x2)
%    %------------------------------
%    subplot(1,7,7);
%    %imagesc(x12, y12, fx_v1v2);
%    %get(hi)
%    %set(hi, 'zscale', 'log')
%    imagesc(x12, y12, log10(pspkx_v1v2+0.001));
%    %cm = colormap(gray);
%    %colormap(gca, flipud(cm));
%    set(gca,'clim', [-3 0]);
%    set(gca,'tickdir', 'out');
%    xlabel('MID_1 Proj');
%    ylabel('MID_2 Proj');
%    axis('xy');
%    %axis('square');
%    hc = colorbar;
%    set(gca,'fontsize', 8);
%    set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
%    title('P(spk|x1,x2)', 'fontsize', 8);

      orient landscape;
      print_mfilename(mfilename);
      set(gcf,'position',[64 617 1087 148]);
end % (for i)












