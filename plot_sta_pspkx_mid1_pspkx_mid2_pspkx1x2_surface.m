function plot_sta_pspx_v1_pspx_v2_pspx1x2_surface(mid, saveit)
%plot_sta_pspx_v1_pspx_v2_pspx_pspx1x2 - Plots tuning curve data
%   obtained using the Michigan 16 channel silicon probe.
%
%   plot_sta_pspx_v1_pspx_v2_pspx_pspx1x2(mid, saveit)
%
%   strf is a struct array holding the receptive field data
%
%   trigger is a vector of trigger times associated with the
%     ripple stimlui
%
%   flim is an optional 1x2 vector, holding the frequency range 
%      over which the strfs are plotted. No input gives plots
%      over the full frequency range.
%
%   tlim is the same as flim, except it over the time axis.
%
%   printit == 1 if you want to print the figures. This 
%   argument is optional. The default is no printing.
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




   % Plot V1
   %------------------------------
   subplot(5,7,(fignum-1)*7+1);
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
      title('V1', 'fontsize', 10);
   end;


   % Plot V2
   %------------------------------
   subplot(5,7,(fignum-1)*7+2);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(v2);
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   set(gca, 'tickdir', 'out');
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('V2', 'fontsize', 10);
   end;


   % Plot F(V1,V2)
   %------------------------------
   subplot(5,7,(fignum-1)*7+3);
   nbins = size(fx_v1v2,1);
   xedges = linspace(min(x12), max(x12), nbins);
   yedges = linspace(min(y12), max(y12), nbins);
   [xmat, ymat] = meshgrid(xedges, yedges);
   surface(xmat, ymat, fx_v1v2);
   axis([min(x12) max(x12) min(y12) max(y12) 0 max(max(fx_v1v2))]);
   %colormap(gray);
   set(gca,'tickdir', 'out');

   % Plot F(V1,V2)
   %------------------------------
   subplot(5,7,(fignum-1)*7+4);
   nbins = size(fx_v1v2,1);
   xedges = linspace(min(x12), max(x12), nbins);
   yedges = linspace(min(y12), max(y12), nbins);
   [xmat, ymat] = meshgrid(xedges, yedges);
   surface(xmat, ymat, fx_v1v2);
   set(gca,'zscale', 'log');
   axis([min(x12) max(x12) min(y12) max(y12) 0 max(max(fx_v1v2))]);
   %colormap(gray);
   set(gca,'tickdir', 'out');



   % Plot F(V1,V2)
   %------------------------------
   subplot(5,7,(fignum-1)*7+5);
   nbins = size(fx_v1v2,1);
   xedges = linspace(min(x12), max(x12), nbins);
   yedges = linspace(min(y12), max(y12), nbins);
   [xmat, ymat] = meshgrid(xedges, yedges);
   mesh(xmat, ymat, fx_v1v2);
   axis([min(x12) max(x12) min(y12) max(y12) 0 max(max(fx_v1v2))]);
   set(gca,'tickdir', 'out');
   if ( fignum == 5 )
      xlabel('MID1 Proj');
      ylabel('MID2 Proj');
   end


   % Plot F(V1,V2)
   %------------------------------
   subplot(5,7,(fignum-1)*7+6);
   nbins = size(fx_v1v2,1);
   xedges = linspace(min(x12), max(x12), nbins);
   yedges = linspace(min(y12), max(y12), nbins);
   [xmat, ymat] = meshgrid(xedges, yedges);
   mesh(xmat, ymat, fx_v1v2);
   axis([min(x12) max(x12) min(y12) max(y12) 0 max(max(fx_v1v2))]);
   set(gca,'zscale', 'log');   
   set(gca,'tickdir', 'out');


   % Plot F(V1,V2)
   %------------------------------
   subplot(5,7,(fignum-1)*7+7);
   %imagesc(x12, y12, fx_v1v2);
   %get(hi)
   %set(hi, 'zscale', 'log')
   imagesc(x12, y12, log10(fx_v1v2+0.001));
   %cm = colormap(gray);
   %colormap(gca, jet);
   set(gca,'clim', [-3 0]);
   set(gca,'tickdir', 'out');
   if ( fignum == 5 )
      xlabel('MID1 Proj');
      ylabel('MID2 Proj');
   end
   axis('xy');
   %axis('square');
   hc = colorbar;
   set(gca,'fontsize', 8);
   set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
   if ( fignum == 1 )
      title('V1/V2 F(x1,x2)', 'fontsize', 10);
   end


cm = colormap(gray);
colormap([0 0 0]);

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
      set(gcf,'position',[50 100 700 800]);
      if ( saveit )
         exportfig(gcf,sprintf('sta_pspx_v1_pspx_v2_pspx_pspx1x2_%.0f.eps',i), 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
         pause(2);
         close;
      end
pause
      fignum = 1;
      if ( i ~= numplots ), figure; end;
   end

end % (for i)






