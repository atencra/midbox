function options = figure1_chapter3_plot_sta_pspkx_mid1_pspkx_mid2_pspkx_pspkx1x2(mid)
%figure1_chapter3_plot_sta_pspx_v1_pspx_v2_pspx_pspx1x2 - Plot the sta, its 
%  nonlinearity, the first mid, its nonlinearity, the second
%  mid, its nonlinearity, and then the 2D nonlinearity.
%
%  options = figure1_chapter3_plot_sta_pspx_v1_pspx_v2_pspx_pspx1x2(mid, saveit)
%
%  options is a struct you can use with the function exportfig to make nice
%  eps figures, which can later be used in Adobe Illustrator.
%
% options is defined as:
%
% options = struct('color', 'rgb', 'height', 11, 'width', 6, ...
%                  'fontmode', 'fixed', 'fontsize', 4, 'bounds', 'tight');
%
% And thus the command to make EPS figures is:
%
% exportfig(gcf,'<figure name>.eps', options);
%
%   caa 6/28/02

if ( nargin ~= 1 )
   error('You need 1 input argument.');
end

if ( ~isfield(mid,'position') )
   error('You need to input the mid variable with the position field.');
end

set(0,'defaultAxesFontName', 'Palatino')
set(0,'defaultAxesFontSize', 10);

close all;

figure;

tbins = mid(1).tbins;
fbins = mid(1).fbins;

ttick = [1 tbins/2 tbins];
tlabel = [100 50 0];

ftick = [1 12 24 36 48];
ftick = ftick(ftick<fbins);
f = [5000 10000 20000 40000];
flabel = f(1:length(ftick));

for i = 1:20 %length(mid)

   if ( i <= length(mid) & i<21 )

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


   % Plot the STA
   %------------------------------
   subplot(20, 7, (i-1)*7+1);
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(sta);
   axis('xy');
   set(gca,'ytick', ftick, 'yticklabel', flabel);
   if ( i == length(mid) )
      set(gca,'xtick', ttick, 'xticklabel', tlabel);
   else
      set(gca,'xtick', ttick, 'xticklabel', '');
   end
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   hy = ylabel(sprintf('%.0fum', mid(i).position));
   set(hy, 'rotation', 0);
   if ( i == 1 )
      title('STA');
   end;


   % Plot the STA's nonlinearity
   %------------------------------
   subplot(20, 7, (i-1)*7+2);
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
   if ( i == 1 )
      title('STA P(spike|x)');
   end;


   % Plot V1
   %------------------------------
   subplot(20, 7, (i-1)*7+3);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v1);
   axis('xy');
   set(gca,'ytick', ftick, 'yticklabel', '');
   if ( i == length(mid) )
      set(gca,'xtick', ttick, 'xticklabel', tlabel);
   else
      set(gca,'xtick', ttick, 'xticklabel', '');
   end
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   if ( i == 1 )
      title('V1');
   end;

   % Plot V1's nonlinearity
   %------------------------------
   subplot(20, 7, (i-1)*7+4);
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
      title('V1 P(spike|x)');
   end;


   % Plot V2
   %------------------------------
   subplot(20, 7, (i-1)*7+5);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v2);
   axis('xy');
   set(gca,'ytick', ftick, 'yticklabel', '');
   if ( i == length(mid) )
      set(gca,'xtick', ttick, 'xticklabel', tlabel);
   else
      set(gca,'xtick', ttick, 'xticklabel', '');
   end
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   if ( i == 1 )
      title('V2');
   end;

   % Plot V2's nonlinearity
   %------------------------------
   subplot(20, 7, (i-1)*7+6);
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
      title('V2 P(spike|x)');
   end;


   % Plot F(V1,V2)
   %------------------------------
   subplot(20, 7, (i-1)*7+7);
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
      title('P(spike|x1,x2)');
   end



   else

      subplot(20, 7, (i-1)*7+1);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

      subplot(20, 7, (i-1)*7+2);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

      subplot(20, 7, (i-1)*7+3);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

      subplot(20, 7, (i-1)*7+4);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

      subplot(20, 7, (i-1)*7+5);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

      subplot(20, 7, (i-1)*7+6);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

      subplot(20, 7, (i-1)*7+7);
      set(gca,'xtick', [], 'xticklabel', '');
      set(gca,'ytick', [], 'yticklabel', '');

   end

end % (for i)

options = struct('color', 'rgb', 'height', 11, 'width', 6, 'fontmode', 'fixed', 'fontsize', 4, 'bounds', 'tight');




















