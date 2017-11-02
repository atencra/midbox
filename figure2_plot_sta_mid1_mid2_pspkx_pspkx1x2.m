function figure2_plot_sta_mid1_mid2_pspkx_pspkx1x2(mid, saveit)
%figure2_plot_sta_v1_v2_pspx_pspx1x2 - Plot only the nonlinearities
%   for the sta, mid1, and mid2 estimates.
%
%   plot_sta_pspx_v1_pspx_v2_pspx_pspx1x2(mid, saveit)
%
% For illustrations in adobe the following is useful:
%
% exportfig(gcf,'sta_v1_v2_fx.eps', 'color', 'rgb', 'height', 8, 'width', 8, 'fontmode', 'fixed', 'fontsize', 6);
%
%   caa 6/28/02

if ( nargin == 1 )
   saveit = 0;
end

if ( nargin == 2 )
   if ( saveit ~= 1 )
      saveit = 0;
   end
end

if ( ~isfield(mid,'position') )
   error('You need to input the mid variable with the position field.');
end

set(0,'defaultAxesFontName', 'Palatino')
set(0,'defaultAxesFontSize', 10);

close all;

ntot = 1;
fignum = 1;
figure;

max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];

for i = 1:length(mid)

   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
   sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   fx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   pspx_sta = fx_sta;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxsp_sta = fx_sta ./ pspike .* px_sta;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   pspx_v1 = fx_v1;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxsp_v1 = fx_v1 ./ pspike .* px_v1;

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
   max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
   max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];

   % Plot the nonlinearity for the STA
   %----------------------------------
   subplot(length(mid), 4, (i-1)*4+1);
   hold on;
   plot(x_sta, fx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 1);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   if ( i~=length(mid) ), set(gca,'xticklabel', ''); end;
   maxfx = str2num(sprintf('%.3f',max(fx_sta)));
   ylim([0 max(fx_sta)]);
   set(gca,'ytick', [0 max(fx_sta)/2 max(fx_sta)], 'yticklabel', [0 maxfx/2 maxfx]);
   set(gca,'tickdir', 'out');
   hy = ylabel(sprintf('%.0fum', mid(i).position));
   set(hy, 'rotation', 0);
   if ( i == 1 )
      title('STA P(sp|x)');
   end;


   % Plot the nonlinearity for the MID1
   %-----------------------------------
   subplot(length(mid), 4, (i-1)*4+2);
   hold on;
   plot(x_v1, fx_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 1);
   plot([min(x_v1) max(x_v1)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   if ( i~=length(mid) ), set(gca,'xticklabel', ''); end;
   maxfx = str2num(sprintf('%.3f',max(fx_v1)));
   ylim([0 max(fx_v1)]);
   set(gca,'ytick', [0 max(fx_v1)/2 max(fx_v1)], 'yticklabel', [0 maxfx/2 maxfx]);
   set(gca,'tickdir', 'out');
   if ( i == 1 )
      title('V1 P(sp|x)');
   end;


   % Plot the nonlinearity for the MID2
   %-----------------------------------
   subplot(length(mid), 4, (i-1)*4+3);
   hold on;
   plot(x_v2, fx_v2, 'ko-', 'markerfacecolor', 'k', 'markersize', 1);
   plot([min(x_v2) max(x_v2)], [pspike pspike], 'k:');
   xlim([-6 6]);
   set(gca,'xtick', [-5 -2.5 0 2.5 5], 'xticklabel', [-5 -2.5 0 2.5]);
   if ( i~=length(mid) ), set(gca,'xticklabel', ''); end;
   maxfx = str2num(sprintf('%.3f',max(fx_v2)));
   ylim([0 max(fx_v2)]);
   set(gca,'ytick', [0 max(fx_v2)/2 max(fx_v2)], 'yticklabel', [0 maxfx/2 maxfx]);
   set(gca,'tickdir', 'out');
   if ( i == 1 )
      title('V2 P(sp|x)');
   end;


   % Plot the nonlinearity for MID1 and MID2
   %----------------------------------------
   subplot(length(mid), 4, (i-1)*4+4);
   imagesc(x12, y12, log10(fx_v1v2+0.001));
   set(gca,'clim', [-3 0]);
   set(gca,'tickdir', 'out');
   if ( i == length(mid) )
      xlabel('MID1 Proj');
      ylabel('MID2 Proj');
   end
   axis('xy');
   if ( i==length(mid) )
      hc = colorbar;
      set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
   end
   if ( i == 1 )
      title('V1/V2 F(x1,x2)');
   end

end % (for i)


options = struct('color', 'rgb', 'height', 11, 'width', 4, 'fontmode', 'fixed', 'fontsize', 4, 'bounds', 'tight');

if ( saveit )
   exportfig(gcf,sprintf('figure2_example_sta_pspx_v1_pspx_v2_pspx_pspx1x2.eps'), options); 
end

% max_min_sta
% max_min_v1
% max_min_v2






if ( nargin == 1 )
   saveit = 0;
end

if ( nargin == 2 )
   if ( saveit ~= 1 )
      saveit = 0;
   end
end

if ( ~isfield(mid,'position') )
   error('You need to input the mid variable with the position field.');
end

set(0,'defaultAxesFontName', 'Palatino')
set(0,'defaultAxesFontSize', 10);

close all;

ntot = 1;
fignum = 1;
figure;


max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];


for i = 1:length(mid)

   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
   sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   fx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   pspx_sta = fx_sta;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxsp_sta = fx_sta ./ pspike .* px_sta;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   pspx_v1 = fx_v1;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxsp_v1 = fx_v1 ./ pspike .* px_v1;

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
   max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
   max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];

   % Plot STA and its nonlinearity
   %------------------------------
   subplot(length(mid), 4, (i-1)*4+1);
   hold on;
   plot(x_sta, fx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 1);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   if ( i~=length(mid) ), set(gca,'xticklabel', ''); end;
   maxfx = str2num(sprintf('%.3f',max(fx_sta)));
   ylim([0 max(fx_sta)]);
   set(gca,'ytick', [0 max(fx_sta)/2 max(fx_sta)], 'yticklabel', [0 maxfx/2 maxfx]);
   set(gca,'tickdir', 'out');
   hy = ylabel(sprintf('%.0fum', mid(i).position));
   set(hy, 'rotation', 0);
   if ( i == 1 )
      title('STA P(sp|x)');
   end;


   % Plot V1 and its nonlinearity
   %------------------------------
   subplot(length(mid), 4, (i-1)*4+2);
   hold on;
   plot(x_v1, fx_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 1);
   plot([min(x_v1) max(x_v1)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   if ( i~=length(mid) ), set(gca,'xticklabel', ''); end;
   maxfx = str2num(sprintf('%.3f',max(fx_v1)));
   ylim([0 max(fx_v1)]);
   set(gca,'ytick', [0 max(fx_v1)/2 max(fx_v1)], 'yticklabel', [0 maxfx/2 maxfx]);
   set(gca,'tickdir', 'out');
   if ( i == 1 )
      title('V1 P(sp|x)');
   end;


   % Plot V2 and its nonlinearity
   %------------------------------
   subplot(length(mid), 4, (i-1)*4+3);
   hold on;
   plot(x_v2, fx_v2, 'ko-', 'markerfacecolor', 'k', 'markersize', 1);
   plot([min(x_v2) max(x_v2)], [pspike pspike], 'k:');
   xlim([-6 6]);
   set(gca,'xtick', [-5 -2.5 0 2.5 5], 'xticklabel', [-5 -2.5 0 2.5]);
   if ( i~=length(mid) ), set(gca,'xticklabel', ''); end;
   maxfx = str2num(sprintf('%.3f',max(fx_v2)));
   ylim([0 max(fx_v2)]);
   set(gca,'ytick', [0 max(fx_v2)/2 max(fx_v2)], 'yticklabel', [0 maxfx/2 maxfx]);
   set(gca,'tickdir', 'out');
   if ( i == 1 )
      title('V2 P(sp|x)');
   end;


   % Plot F(V1,V2)
   %------------------------------
   subplot(length(mid), 4, (i-1)*4+4);
   imagesc(x12, y12, log10(fx_v1v2+0.001));
   set(gca,'clim', [-3 0]);
   set(gca,'tickdir', 'out');
   if ( i == length(mid) )
      xlabel('MID1 Proj');
      ylabel('MID2 Proj');
   end
   axis('xy');
   if ( i==length(mid) )
      hc = colorbar;
      set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
   end
   if ( i == 1 )
      title('V1/V2 F(x1,x2)');
   end

end % (for i)


options = struct('color', 'rgb', 'height', 11, 'width', 4, 'fontmode', 'fixed', 'fontsize', 4, 'bounds', 'tight');

if ( saveit )
   exportfig(gcf,sprintf('figure2_example_sta_pspx_v1_pspx_v2_pspx_pspx1x2.eps'), options); 
end

% max_min_sta
% max_min_v1
% max_min_v2






