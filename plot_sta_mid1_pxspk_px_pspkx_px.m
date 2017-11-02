function plot_sta_mid1_pxspk_px_pspkx_px(mid, saveit, celltype)
%plot_sta_mid1_pxspk_px_pspkx_px - 
%
%
% For illustrations in adobe the following is useful:
%
% exportfig(gcf,'sta_mid1_pxspk_px_pspkx_px.eps', 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
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


ntot = 1;
fignum = 1;
figure;

numplots = ceil(length(mid)/4);

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

% index = find( pxsp_v1>0 & px_v1>0 )
% sum( pxsp_v1(index) .* log2( pxsp_v1(index) ./ px_v1(index) ) )

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

% index = find( pxsp_v2>0 & px_v2>0 )
% sum( pxsp_v2(index) .* log2( pxsp_v2(index) ./ px_v2(index) ) )
% pause

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;


   % Plot STA
   %------------------------------
   subplot(5,7,(fignum-1)*7+1);
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(sta);
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( isfield(mid,'location') )
      ylabel(sprintf('%.0f : %.0f - %.0f', i, mid(i).location, mid(i).unit));
   else
      ylabel(sprintf('%s\n%.0f-%.0f', mid(i).exp, mid(i).site, mid(i).position));
   end
   if ( fignum == 1 )
      title('STA', 'fontsize', 8);
   end;


   % Plot MID1
   %------------------------------
   subplot(5,7,(fignum-1)*7+2);
   v1(:,end) = 0.9 * (rand(length(v1(:,end)),1)-0.5);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v1);
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1', 'fontsize', 8);
   end;



   % Plot P(x|spk) and P(x) for the STA
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+3);
   hold on;
   plot(x_sta, px_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot(x_sta, pxsp_sta, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
   if ( fignum == 1 )
      legend('P(x)', 'P(x|spk)', 0);
   end

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f', 1.025*max([max(px_sta) max(pxsp_sta)]) ));
   maxfx_div2 = str2num(sprintf('%.3f',maxfx/2));
   ylim([0 maxfx]);
   set(gca,'ytick', [0 maxfx/2 maxfx], 'yticklabel', [0 maxfx_div2 maxfx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('STA', 'fontsize', 8);
   end;


   % Plot P(x|sp) and P(x) for the STA
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+4);
   hold on;
   plot(x_v1, px_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot(x_v1, pxsp_v1, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f', 1.025*max([max(px_v1) max(pxsp_v1)]) ));
   maxfx_div2 = str2num(sprintf('%.3f',maxfx/2));
   ylim([0 maxfx]);
   set(gca,'ytick', [0 maxfx/2 maxfx], 'yticklabel', [0 maxfx_div2 maxfx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1', 'fontsize', 8);
   end;


   % Plot P(spike|x) for the STA and the MID1 on the same axes
   %----------------------------------------------------------
   subplot(5,7,(fignum-1)*7+5);
   hold on;
   plot(x_sta, fx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot(x_v1, fx_v1, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
   plot([min([x_sta(:); x_v1(:)]) max([x_sta(:); x_v1(:)])], [pspike pspike], 'k:');

   if ( fignum == 1 )
      legend('STA', 'MID_1', 0);
   end

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f', 1.025*max([max(fx_sta) max(fx_v1)]) ));
   maxfx_div2 = str2num(sprintf('%.3f',maxfx/2));
   ylim([0 maxfx]);
   set(gca,'ytick', [0 maxfx/2 maxfx], 'yticklabel', [0 maxfx_div2 maxfx]);

   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('P(spk|x)', 'fontsize', 8);
   end;


   % Plot P(spike|x), P(x) for the STA
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+6);
   hold on;
   plot(x_sta, px_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot(x_sta, pspx_sta, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
   if ( fignum == 1 )
      legend('P(x)', 'P(spk|x)', 0);
   end

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f', 1.025*max([max(pspx_sta) max(px_sta)]) ));
   maxfx_div2 = str2num(sprintf('%.3f',maxfx/2));
   ylim([0 maxfx]);
   set(gca,'ytick', [0 maxfx/2 maxfx], 'yticklabel', [0 maxfx_div2 maxfx]);

   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('STA', 'fontsize', 8);
   end;


   % Plot P(x|spike), P(x) for the MID1
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+7);
   hold on;
   plot(x_v1, px_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot(x_v1, pspx_v1, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f', 1.025*max([max(pspx_v1) max(px_v1)]) ));
   maxfx_div2 = str2num(sprintf('%.3f',maxfx/2));
   ylim([0 maxfx]);
   set(gca,'ytick', [0 maxfx/2 maxfx], 'yticklabel', [0 maxfx_div2 maxfx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1', 'fontsize', 8);
   end;

   if ( mod(fignum,5) )
      fignum = fignum + 1;
      if ( saveit )
         if ( i == length(mid) )
            orient landscape;
            print_mfilename(mfilename);
            if ( nargin==3 )
               filename = sprintf('%s_sta_mid1_pxspk_px_pspkx_px_%.0f.eps', celltype, i);
               exportfig(gcf,filename, 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
               pause(2);
            else
               exportfig(gcf,sprintf('sta_mid1_pxspk_px_pspkx_px_%.0f.eps',i), 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
            end
            pause(2);

         end
      end
   else
      orient landscape;
      print_mfilename(mfilename);
      set(gcf,'position',[50 100 700 800]);
      if ( saveit )
         if ( nargin==3 )
            filename = sprintf('%s_sta_mid1_pxspk_px_pspkx_px_%.0f.eps', celltype, i);
            exportfig(gcf,filename, 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
            pause(1);
            close;
         else
            exportfig(gcf,sprintf('sta_mid1_pxspk_px_pspkx_px_%.0f.eps',i), 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
            pause(1);
            close;
         end
      else
         pause;
      end
      fignum = 1;
      if ( i ~= numplots ), figure; end;
   end

end % (for i)







