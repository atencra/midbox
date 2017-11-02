function plot_sta_pspkx_mid1_pspkx_px(mid, saveit, celltype)
%plot_sta_v1_fx_px - Plot the STA, the STA nonlinearity, the MID1, the
% MID1 nonlinearity, both nonlinearities on the same axis, then
% then the nonlinearity and the prior distribution for the STA and the
% MID1.
%
% plot_sta_pspkx_v1_pspkx_px(mid, saveit, celltype)
% -------------------------------------------------
% Each figure has 5 rows, one row per neuron. The figure has 7 columns,
% corresponding to the plots mentioned above.
%
% celltype is an optional string argument used to give appropriate titles
% to the output figures, and may be either 'fsu' or 'rsu'
%
% For illustrations in adobe the following is useful:
%
% exportfig(gcf,'sta_pspkx_v1_pspkx_px.eps', 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
%
% caa 6/28/02


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
   pspkx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxspk_sta = pspkx_sta ./ pspike .* px_sta;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   pspkx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxspk_v1 = pspkx_v1 ./ pspike .* px_v1;

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

   % Plot the STA nonlinearity
   %------------------------------
   subplot(5,7,(fignum-1)*7+2);
   hold on;
   plot(x_sta, pspkx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   mx = str2num(sprintf('%.3f', 1.025*max(pspkx_sta)));
   mx2 = str2num(sprintf('%.3f', mx/2));
   ylim([0 mx]);
   set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('STA P(spk|x)', 'fontsize', 8);
   end;


   % Plot MID1
   %------------------------------
   subplot(5,7,(fignum-1)*7+3);
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

   % Plot the MID1 nonlinearity
   %------------------------------
   subplot(5,7,(fignum-1)*7+4);
   hold on;
   plot(x_v1, pspkx_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot([min(x_v1) max(x_v1)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   mx = str2num(sprintf('%.3f', 1.025*max(pspkx_v1)));
   mx2 = str2num(sprintf('%.3f', mx/2));
   ylim([0 mx]);
   set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1 P(spk|x)', 'fontsize', 8);
   end;

   % Plot P(spike|x) for the STA and MID1
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+5);
   hold on;
   plot(x_sta, pspkx_sta, 'k-', 'linewidth', 2); %'markerfacecolor', 'k', 'markersize', 2);
   plot(x_v1, pspkx_v1, 'r-', 'linewidth', 2); %'markerfacecolor', 'r', 'markersize', 2);

   if ( fignum == 1 )
      legend('STA', 'MID_1', 0);
   end

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   mx = str2num(sprintf('%.3f', 1.025*max([pspkx_sta(:); pspkx_v1(:)]) ));
   mx2 = str2num(sprintf('%.3f',mx/2));
   ylim([0 mx]);
   set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);

   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('P(spk|x)', 'fontsize', 8);
   end


   % Plot P(spike|x), P(x) for the STA
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+6);
   hold on;
   plot(x_sta, px_sta, 'k-', 'linewidth', 2); %'markerfacecolor', 'r', 'markersize', 2);
   plot(x_sta, pspkx_sta, 'r-', 'linewidth', 2); %'markerfacecolor', 'k', 'markersize', 2);
%    plot(x_sta, pspx_sta, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
%    plot(x_sta, px_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);

   if ( fignum == 1 )
      legend('P(x)', 'P(spk|x)', 0);
   end

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   mx = str2num(sprintf('%.3f', 1.025*max([pspkx_sta(:); px_sta(:)]) ));
   mx2 = str2num(sprintf('%.3f',mx/2));
   ylim([0 mx]);
   set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);

   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('STA', 'fontsize', 8);
   end;


   % Plot P(spike|x), P(x) for the V1
   %----------------------------------------------
   subplot(5,7,(fignum-1)*7+7);
   hold on;
   plot(x_v1, px_v1, 'k-', 'linewidth', 2); %'markerfacecolor', 'r', 'markersize', 2);
   plot(x_v1, pspkx_v1, 'r-', 'linewidth', 2); %'markerfacecolor', 'k', 'markersize', 2);
%    plot(x_v1, pspx_v1, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
%    plot(x_v1, px_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);

   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   mx = str2num(sprintf('%.3f', 1.025*max([pspkx_v1(:); px_v1(:)]) ));
   mx2 = str2num(sprintf('%.3f',mx/2));
   ylim([0 mx]);
   set(gca,'ytick', [0 mx2 mx], 'yticklabel', [0 mx2 mx]);
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
               filename = sprintf('%s_sta_pspkx_v1_pspkx_px_%.0f.eps', celltype, i);
               exportfig(gcf,filename, 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
               pause(2);
            else
               exportfig(gcf,sprintf('sta_pspkx_v1_pspkx_px_%.0f.eps',i), 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
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
            filename = sprintf('%s_sta_pspkx_v1_pspkx_px_%.0f.eps', celltype, i);
            exportfig(gcf,filename, 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
            pause(1);
            close;
         else
            exportfig(gcf,sprintf('sta_pspkx_v1_pspkx_px_%.0f.eps',i), 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
            pause(1);
            close;
         end
      else
         pause
      end
      fignum = 1;
      if ( i ~= numplots ), figure; end;
   end

end % (for i)









