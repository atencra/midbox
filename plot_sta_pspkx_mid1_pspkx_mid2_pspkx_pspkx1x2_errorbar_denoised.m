function plot_sta_pspkx_mid1_pspkx_mid2_pspkx_pspkx1x2_errorbar_denoised(mid, saveit)
%plot_sta_pspx_v1_pspx_v2_pspx_pspx1x2 - Plots tuning curve data
%   obtained using the Michigan 16 channel silicon probe.
%
%
% For illustrations in adobe the following is useful:
%
% exportfig(gcf,'sta_v1_v2_fx.eps', 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
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

max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];

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
   fx_se_sta = mid(i).rpx1pxpxt_sta.ior_std ./ 4;
   pspx_sta = fx_sta;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxsp_sta = fx_sta ./ pspike .* px_sta;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   fx_se_v1 = mid(i).rpdx1x2px_pxt_2.ior1_std ./ 4;
   pspx_v1 = fx_v1;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxsp_v1 = fx_v1 ./ pspike .* px_v1;

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   fx_se_v2 = mid(i).rpdx1x2px_pxt_2.ior2_std ./ 4;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
   max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
   max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];

   % Plot the STA
   %------------------------------
   subplot(5,7,(fignum-1)*7+1);
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(sta);
%    imagesc(sta(12:24,:));
   axis('xy');
   %set(gca,'ydir', 'normal');
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   hc = colorbar;
   ylm = get(hc, 'ylim');
   set(hc, 'ytick', [ylm(1) 0 ylm(2)]);
   set(hc, 'yticklabel', [ylm(1) 0 ylm(2)]);
   set(hc,'tickdir', 'out');
   set(hc,'fontsize', 8);
   set(gca,'fontsize', 8);
   ylabel(sprintf('%.0f : %.0f - %.0f', i, mid(i).location, mid(i).unit));
   if ( fignum == 1 )
      title('STA', 'fontsize', 8);
   end;

   % Plot the STA nonlinearity
   %------------------------------
   subplot(5,7,(fignum-1)*7+2);
   hold on;

   he = errorbar(x_sta, fx_sta, fx_se_sta, 'ko-');
   set(he(2), 'color', [0 0 0], 'linewidth', 0.5, 'markersize', 2, 'markerfacecolor', 'k');
   set(he(1), 'color', [0 0 0], 'linewidth', 0.5);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f', max(1.025 * [fx_sta + fx_se_sta])));
   maxfx_div2 = str2num(sprintf('%.3f', max(1.025 * [fx_sta + fx_se_sta]) / 2));
   ylim([0 max(1.05*[fx_sta + fx_se_sta])]);
   set(gca,'ytick', [0 max(fx_sta)/2 max(fx_sta)], 'yticklabel', [0 maxfx_div2 maxfx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('STA P(spk|proj)', 'fontsize', 8);
   end;


%    plot(x_sta, fx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
%    plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
%    xlim([-7.5 7.5]);
%    set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
%    maxfx = str2num(sprintf('%.3f',max(fx_sta)));
%    maxfx_div2 = str2num(sprintf('%.3f',max(fx_sta)/2));
%    ylim([0 max(fx_sta)]);
%    set(gca,'ytick', [0 max(fx_sta)/2 max(fx_sta)], 'yticklabel', [0 maxfx_div2 maxfx]);
%    set(gca,'tickdir', 'out');
%    set(gca,'fontsize', 8);
%    if ( fignum == 1 )
%       title('STA P(spk|proj)', 'fontsize', 8);
%    end




   % Plot the MID1
   %------------------------------
   subplot(5,7,(fignum-1)*7+3);
%    v1(:,end) = 0.9 * (rand(length(v1(:,end)),1)-0.5);
   v1(:,end) = 0; %0.5 * (rand(length(v2(:,end)),1)-0.5);
   v1(:,1) = 0; %0.5 * (rand(length(v2(:,end)),1)-0.5);
   v1(end,:) = 0; %0.5 * (rand(1,length(v2(1,:)))-0.5);
   v1(end-1,:) = 0; %0.5 * (rand(1,length(v2(2,:)))-0.5);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v1);
%    imagesc(v1(12:24,:));
   axis('xy');
   %set(gca,'ydir', 'normal');
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   hc = colorbar;
   ylm = get(hc, 'ylim');
   set(hc, 'ytick', [ylm(1) 0 ylm(2)]);
   set(hc, 'yticklabel', [ylm(1) 0 ylm(2)]);
   set(hc,'tickdir', 'out');
   set(hc,'fontsize', 8);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1', 'fontsize', 8);
   end;

   % Plot the MID1 nonlinearity
   %------------------------------
   subplot(5,7,(fignum-1)*7+4);
   hold on;
   fx_se_v1(isinf(fx_se_v1)) = max(fx_se_v1(~isinf(fx_se_v1)));
   he = errorbar(x_v1, fx_v1, fx_se_v1, 'ko-');
   set(he(2), 'color', [0 0 0], 'linewidth', 0.5, 'markersize', 2, 'markerfacecolor', 'k');
   set(he(1), 'color', [0 0 0], 'linewidth', 0.5);
   plot([min(x_v1) max(x_v1)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   no_inf_fx = fx_v1(~isinf(fx_v1));
   no_inf_fx_se = fx_se_v1(~isinf(fx_se_v1));
   maxfx = max([ max(no_inf_fx) max(no_inf_fx_se) ]);
   maxfx = str2num(sprintf('%.3f', maxfx ) );
   maxfx_div2 = str2num(sprintf('%.3f', maxfx / 2));
   ylim([0 max(1.05*[fx_v1 + fx_se_v1])]);
   set(gca,'ytick', [0 max(fx_v1)/2 max(fx_v1)], 'yticklabel', [0 maxfx_div2 maxfx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);

%    plot(x_v1, fx_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
%    plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
%    xlim([-7.5 7.5]);
%    set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
%    maxfx = str2num(sprintf('%.3f',max(fx_v1)));
%    maxfx_div2 = str2num(sprintf('%.3f',max(fx_v1)/2));
%    ylim([0 max(fx_v1)]);
%    set(gca,'ytick', [0 max(fx_v1)/2 max(fx_v1)], 'yticklabel', [0 maxfx_div2 maxfx]);
%    set(gca,'tickdir', 'out');
%    set(gca,'fontsize', 8);
%    if ( fignum == 1 )
%       title('MID_1 P(spk|proj)', 'fontsize', 8);
%    end;

   if ( fignum == 1 )
      title('MID_1 P(spk|proj)', 'fontsize', 8);
   end;





   % Plot the MID2
   %------------------------------
   subplot(5,7,(fignum-1)*7+5);
   v2(:,end) = 0; %0.5 * (rand(length(v2(:,end)),1)-0.5);
   v2(:,1) = 0; %0.5 * (rand(length(v2(:,end)),1)-0.5);
   v2(end,:) = 0; %0.5 * (rand(1,length(v2(1,:)))-0.5);
   v2(end-1,:) = 0; %0.5 * (rand(1,length(v2(2,:)))-0.5);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v2);
%    imagesc(v2(12:24,:));
   axis('xy');
   %set(gca,'ydir', 'normal');
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   hc = colorbar;
   ylm = get(hc, 'ylim');
   set(hc, 'ytick', [ylm(1) 0 ylm(2)]);
   set(hc, 'yticklabel', [ylm(1) 0 ylm(2)]);
   set(hc,'tickdir', 'out');
   set(hc,'fontsize', 8);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_2', 'fontsize', 8);
   end;

   % Plot the MID2 nonlinearity
   %------------------------------
   subplot(5,7,(fignum-1)*7+6);
   hold on;
   fx_se_v2(isinf(fx_se_v2)) = max(fx_se_v2(~isinf(fx_se_v2)));
   he = errorbar(x_v2, fx_v2, fx_se_v2, 'ko-');
   set(he(2), 'color', [0 0 0], 'linewidth', 0.5, 'markersize', 2, 'markerfacecolor', 'k');
   set(he(1), 'color', [0 0 0], 'linewidth', 0.5);
   plot([min(x_v2) max(x_v2)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   no_inf_fx = fx_v2(~isinf(fx_v2));
   no_inf_fx_se = fx_se_v2(~isinf(fx_se_v2));
   maxfx = max([ max(no_inf_fx) max(no_inf_fx_se) ]);
   maxfx = str2num(sprintf('%.3f', maxfx ) );
   maxfx_div2 = str2num(sprintf('%.3f', maxfx / 2));
   ylim([0 max(1.05*[fx_v2 + fx_se_v2])]);
   set(gca,'ytick', [0 max(fx_v2)/2 max(fx_v2)], 'yticklabel', [0 maxfx_div2 maxfx]);
   set(gca,'tickdir', 'out');
   set(gca,'fontsize', 8);


%    plot(x_v2, fx_v2, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
%    plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
%    xlim([-6 6]);
%    set(gca,'xtick', [-5 -2.5 0 2.5 5], 'xticklabel', [-5 -2.5 0 2.5]);
%    maxfx = str2num(sprintf('%.3f',max(fx_v2)));
%    maxfx_div2 = str2num(sprintf('%.3f',max(fx_v2)/2));
%    ylim([0 max(fx_v2)]);
%    set(gca,'ytick', [0 max(fx_v2)/2 max(fx_v2)], 'yticklabel', [0 maxfx_div2 maxfx]);
%    set(gca,'tickdir', 'out');
%    set(gca,'fontsize', 8);

   if ( fignum == 1 )
      title('MID_2 P(spk|proj)', 'fontsize', 8);
   end;


   % Plot the MID1 and MID2 nonlinearity
   %------------------------------------
   subplot(5,7,(fignum-1)*7+7);
   %imagesc(x12, y12, fx_v1v2);
   %get(hi)
   %set(hi, 'zscale', 'log')
   dx = diff(x12);
   dx = dx(1);
   dy = diff(y12);
   dy = dy(1);
   x12 = [-dx+x12(1) x12(1:7) 0 x12(8:end)];
   y12 = [-dy+y12(1) y12(1:7) 0 y12(8:end)];
   nonlinearity = 0.001 + zeros( size(fx_v1v2,1)+1, size(fx_v1v2,1)+1 );
   nonlinearity(2:end, 2:end) = nonlinearity(2:end, 2:end) + fx_v1v2; 
   
   h = imagesc(x12, y12, log10(nonlinearity));

%   h = imagesc(x12, y12, log10(fx_v1v2+0.001));


% x12
% y12
% size(fx_v1v2)
% pause


%    get(h,'xdata')
%    pause
   %cm = colormap(gray);
   %colormap(gca, flipud(cm));
   set(gca,'clim', [-3 0]);
   set(gca,'tickdir', 'out');
   if ( fignum == 5 )
      xlabel('MID1 Proj','fontsize',8);
      ylabel('MID2 Proj', 'fontsize',8);
   end
   axis('xy');
   %axis('square');
   hc = colorbar;
   set(hc,'fontsize', 8);
   set(gca,'fontsize', 8);
   set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
   set(hc,'tickdir', 'out');
   if ( fignum == 1 )
      title('P(spk|x1,x2)', 'fontsize', 8);
   end


   if ( mod(fignum,5) )
      fignum = fignum + 1;
      if ( saveit )
         if ( i == length(mid) )
            exportfig(gcf,sprintf('sta_pspx_v1_pspx_v2_pspx_pspx1x2_%.0f.eps',i), 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
            pause(2);
         end
      end
   else
      orient landscape;
      print_mfilename(mfilename);
      %set(gcf,'position',[50 100 700 800]);
      if ( saveit )
         exportfig(gcf,sprintf('sta_pspx_v1_pspx_v2_pspx_pspx1x2_%.0f.eps',i), 'color', 'rgb', 'height', 6, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
         pause(1);
         close;
      end

      fignum = 1;
      if ( i ~= numplots ), eval(['f' num2str(i) '=figure;']); end;
   end

   pause(0.5)

end % (for i)






