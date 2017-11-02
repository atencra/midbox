function plot_fio_sta_mid1_mid2(fio, sr)
%plot_fio - plots nonlinearity and probability distributions for fio data.
%
% plot_fio(fio)
% -------------------------------------------------------------------
% Nonlinearities are plotted as p(spk|x) vs. projection (sd).
% 
% fio : struct array holding nonlinearity data. Each element of fio
% holds the data for 1 neuron. Usually saved in a file such as:
%
%    load 2003-11-24-site7-2380um-40db-dmr1-fs18115-mid-fio-projinfo.mat
%
% sr : spike train sampling rate for the analysis. If the spikes were
% binned with 1 ms resolution, then sr = 1000. If the resolution was
% 5 ms, then sr = 200. When a value for sr is supplied, a second figure
% of nonlinearities is plotted, this time with the ordinate in units
% of spikes / second.
%
% caa 2/1/10


close all;

ntot = 1;
len = length(fio);
fignum = 1;
figure;

for i = 1:length(fio)

   location = fio(i).location;
   unit = fio(i).unit;

   xbins = fio(i).x0bins;
   pspk = fio(i).pspk;
   pspkx0 = fio(i).pspkx0;
   pspkx1 = fio(i).pspkx1;
   pspkx2 = fio(i).pspkx2;


   pspkx0_mtx = [];
   pspkx1_mtx = [];
   pspkx2_mtx = [];
   pspk_mtx = [];
   for j = 1:length(pspkx0)
      pspkx0_mtx = [pspkx0_mtx pspkx0{j}];
      pspkx1_mtx = [pspkx1_mtx pspkx1{j}];
      pspkx2_mtx = [pspkx2_mtx pspkx2{j}];
      pspk_mtx = [pspk_mtx pspk{j}];
   end

   pspk_mn = nanmean(pspk_mtx,2);
   pspkx0_mn = nanmean(pspkx0_mtx,2);
   pspkx0_std = nanstd(pspkx0_mtx,0,2);
   pspkx1_mn = nanmean(pspkx1_mtx,2);
   pspkx1_std = nanstd(pspkx1_mtx,0,2);
   pspkx2_mn = nanmean(pspkx2_mtx,2);
   pspkx2_std = nanstd(pspkx2_mtx,0,2);

   index = 1:length(pspkx0_mn);
   xbins = xbins(index);
   pspkx0_mn = pspkx0_mn(index);
   pspkx0_std = pspkx0_std(index);
   pspkx1_mn = pspkx1_mn(index);
   pspkx1_std = pspkx1_std(index);
   pspkx2_mn = pspkx2_mn(index);
   pspkx2_std = pspkx2_std(index);



   max0 = sr * max( [ pspkx0_mn + pspkx0_std] );
   max1 = sr * max( [ pspkx1_mn + pspkx1_std] );
   max2 = sr * max( [ pspkx2_mn + pspkx2_std] );
   maxmax = max([max0 max1]);

   ytick = [0 round(maxmax/2) 2*round(maxmax/2)];
   ytick2 = [0 round(max2/2) 2*round(max2/2)];


   subplot(5,3,fignum);
   hold on;
   plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
   hp = ebar(xbins, sr * pspkx0_mn, sr * pspkx0_std);
   set(hp, 'markersize', 2);
   set(hp,'markerfacecolor', 'k', 'markeredgecolor', 'k');
   xlim([-8 8]);
   ylim([0 1.05 * maxmax]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
   set(gca,'ytick', ytick', 'yticklabel', ytick);

   ylabel(sprintf('%.0f - %.0f', location, unit));

   if ( fignum == 1 )
      title('STA Nonlinearity');
   end
   if ( fignum ~= 13 )
      set(gca, 'xticklabel', []);
   end
%    if ( fignum == 9 )
%       xlabel('Projection (SD)');
%    end


   subplot(5,3,fignum+1);
   hold on;
   plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
   hp = ebar(xbins, sr * pspkx1_mn, sr * pspkx1_std);
   set(hp, 'markersize', 2);
   set(hp,'markerfacecolor', 'k', 'markeredgecolor', 'k');
   xlim([-8 8]);
   ylim([0 1.05 * maxmax]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
   set(gca,'ytick', ytick', 'yticklabel', ytick);

   if ( fignum == 1 )
      title('MID1 Nonlinearity');
   end
   if ( fignum ~= 13 )
      set(gca, 'xticklabel', []);
   end
   if ( fignum == 13 )
      xlabel('Projection (SD)');
   end


   subplot(5,3,fignum+2);
   hold on;
   plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
   hp = ebar(xbins, sr * pspkx2_mn, sr * pspkx2_std);
   set(hp, 'markersize', 2);
   set(hp,'markerfacecolor', 'k', 'markeredgecolor', 'k');
   xlim([-8 8]);
   ylim([0 1.05 * maxmax]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
   set(gca,'ytick', ytick', 'yticklabel', ytick);

   if ( fignum == 1 )
      title('MID2 Nonlinearity');
   end
   if ( fignum ~= 13 )
      set(gca, 'xticklabel', []);
   end
%    if ( fignum == 13 )
%       xlabel('Projection (SD)');
%    end





   if ( mod(fignum,13) )
      fignum = fignum + 3;
      figuretitle = 0;
   else
      orient tall;
      suptitle(sprintf('ICC Data'));
      print_mfilename(mfilename);
      % subplotspace('vertical', 20);
      set(gcf,'position', [725   319   484   599]);
      figuretitle = 1;

      if ( i ~= len )
         fignum = 1;
         figure;
      end
   end

   ntot = ntot + 1;

end


if ( figuretitle == 0 )
   orient tall;
   suptitle(sprintf('ICC Data'));
   print_mfilename(mfilename);
   set(gcf,'position', [725   319   484   599]);
end



return;































