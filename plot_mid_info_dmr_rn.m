function plot_mid_info_dmr_rn(dmrMid, rnMid)
% plot_mid_info_dmr_rn DMR/RN STA, MID1, and MID2 filters
%
%    plot_mid_info_dmr_rn(dmrMid, rnMid) displays the STA,
%    MID1, and MID2 for dynamic moving ripple and ripple noise
%    stimuli. Each figure holds the results for one neuron.
%    Each column in each figure holds the data from DMR or RN.
% 
%    dmrMid is a struct array holding MID data for dmr stimulation.
%    rnMid is a struct array that holds data for ripple noise.
% 
%    caa 8/1/2012
% 
%    plot_filtstr_dmr_rn(dmrMid, rnMid)

mapscheme = 'rdbu';

if ( nargin ~= 2 )
   error('You need 2 input args.');
end

if ( length(dmrMid) ~= length(rnMid) )
   error('Input args must have the same length.');
end


% time = 1000 .* dmrFiltStr(1).time; % change to ms
% freq = dmrFiltStr(1).freq / 1000;
% 
% ytick = [1 floor(length(freq)/2) length(freq)];
% yticklabel = round(freq(ytick)*10)/10;
% 
% xtick = [1 floor(length(time)/2)+1 length(time)];
% xticklabel = round( time(xtick) );

close all;

for i = 1:length(dmrMid)

   exp = dmrMid(i).exp;
   site = dmrMid(i).site;

   dmrunit = dmrMid(i).unit;
   rnunit = rnMid(i).unit;


   % Dynamic Moving Ripple data:
   % ====================================

   % Filters
   dmrSTA = dmrMid(i).rpsta.filter;
   dmrSTA = fliplr(dmrSTA);
   dmrSTA = smoothmat(dmrSTA);

   dmrV1 = dmrMid(i).rpdtest2_v1.filter;
   dmrV1 = fliplr(dmrV1);
   dmrV1 = smoothmat(dmrV1);

   dmrV2 = dmrMid(i).rpdtest2_v2.filter;
   dmrV2 = fliplr(dmrV2);
   dmrV2 = smoothmat(dmrV2);

   % Nonlinearities
   xDmrSta = dmrMid(i).rpx1pxpxt_sta.x;
   fDmrSta = dmrMid(i).rpx1pxpxt_sta.ior_mean;

   xDmrV1 = dmrMid(i).rpdx1x2px_pxt_2.x1;
   fDmrV1 = dmrMid(i).rpdx1x2px_pxt_2.ior1_mean;

   xDmrV2 = dmrMid(i).rpdx1x2px_pxt_2.x2;
   fDmrV2 = dmrMid(i).rpdx1x2px_pxt_2.ior2_mean;

   xDmrSta = xDmrSta(2:end-1);
   fDmrSta = fDmrSta(2:end-1) ./ 0.005;

   xDmrV1 = xDmrV1(2:end-1);
   fDmrV1 = fDmrV1(2:end-1) ./ 0.005;

   xDmrV2 = xDmrV2(3:end-2);
   fDmrV2 = fDmrV2(3:end-2) ./ 0.005;



   % Ripple Noise data:
   % ====================================

   % Filters
   rnSTA = rnMid(i).rpsta.filter;
   rnSTA = fliplr(rnSTA);
   rnV1 = rnMid(i).rpdtest2_v1.filter;
   rnV1 = fliplr(rnV1);
   rnV2 = rnMid(i).rpdtest2_v2.filter;
   rnV2 = fliplr(rnV2);

   % Nonlinearities
   xRnSta = rnMid(i).rpx1pxpxt_sta.x;
   fRnSta = rnMid(i).rpx1pxpxt_sta.ior_mean;

   xRnV1 = rnMid(i).rpdx1x2px_pxt_2.x1;
   fRnV1 = rnMid(i).rpdx1x2px_pxt_2.ior1_mean;

   xRnV2 = rnMid(i).rpdx1x2px_pxt_2.x2;
   fRnV2 = rnMid(i).rpdx1x2px_pxt_2.ior2_mean;

   xRnSta = xRnSta(2:end-1);
   fRnSta = fRnSta(2:end-1) ./ 0.005;

   xRnV1 = xRnV1(2:end-1);
   fRnV1 = fRnV1(2:end-1) ./ 0.005;

   xRnV2 = xRnV2(3:end-3);
   fRnV2 = fRnV2(3:end-3) ./ 0.005;

   

   figure;

   subplot(3,3,1);
   imagesc( dmrSTA );
%    imagesc(time, freq, dmrSTA );
   axis xy;
   minmin = min(min(dmrSTA));
   maxmax = max(max(dmrSTA));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
%    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
%    tickpref;
   set(gca,'xtick', [], 'xticklabel', []);
   set(gca,'ytick', [], 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   h = title(sprintf('DMR\n%s site %.0f unit %.0f', exp, site, dmrunit));
   set(h, 'fontweight', 'bold');
   h = ylabel('STA');
   pos = get(h, 'position');
   pos(1) = -2.5;
   set(h, 'rotation', 0, 'fontweight', 'bold', 'position', pos);


   subplot(3,3,4);
   imagesc( dmrV1 );
   axis xy;
   minmin = min(min(dmrV1));
   maxmax = max(max(dmrV1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
%    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
%    tickpref;
   set(gca,'xtick', [], 'xticklabel', []);
   set(gca,'ytick', [], 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   colormap(cmap);
   h = ylabel('MID1');
   pos = get(h, 'position');
   pos(1) = -2.5;
   set(h, 'rotation', 0, 'fontweight', 'bold', 'position', pos);


   subplot(3,3,7);
   imagesc( dmrV2 );
   axis xy;
   minmin = min(min(dmrV2));
   maxmax = max(max(dmrV2));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
%    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
%    tickpref;
   set(gca,'xtick', [], 'xticklabel', []);
   set(gca,'ytick', [], 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   colormap(cmap);
   h = ylabel('MID2');
   pos = get(h, 'position');
   pos(1) = -2.5;
   set(h, 'rotation', 0, 'fontweight', 'bold', 'position', pos);
%    xlabel('Time (ms)');



   subplot(3,3,2);
   imagesc( rnSTA );
%    imagesc(time, freq, rnSTA );
   axis xy;
   minmin = min(min(rnSTA));
   maxmax = max(max(rnSTA));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
%    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
%    tickpref;
   set(gca,'xtick', [], 'xticklabel', []);
   set(gca,'ytick', [], 'yticklabel', []);
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   h = title(sprintf('RN\n%s site %.0f unit %.0f', exp, site, rnunit));
   set(h, 'fontweight', 'bold');


   subplot(3,3,5);
   imagesc( rnV1 );
   axis xy;
   minmin = min(min(rnV1));
   maxmax = max(max(rnV1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
%    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
%    tickpref;
   set(gca,'xtick', [], 'xticklabel', []);
   set(gca,'ytick', [], 'yticklabel', []);
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   colormap(cmap);

   subplot(3,3,8);
   imagesc( rnV2 );
   axis xy;
   minmin = min(min(rnV2));
   maxmax = max(max(rnV2));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
%    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
%    tickpref;
   set(gca,'xtick', [], 'xticklabel', []);
   set(gca,'ytick', [], 'yticklabel', []);
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   colormap(cmap);
%    xlabel('Time (ms)');



   % Nonlinearities
   subplot(3,3,3);
   hold on;
   plot(xDmrSta, fDmrSta, 'ko-', 'markerfacecolor', 0.6*ones(1,3), ...
      'markersize', 3);
   plot(xRnSta, fRnSta, 'ks--', 'markerfacecolor', 0*ones(1,3), ...
      'markersize', 3);
   tickpref;
   legend('DMR', 'RN', 'location', 'northwest');
   title('Nonlinearities', 'fontweight', 'bold');


   subplot(3,3,6);
   hold on;
   plot(xDmrV1, fDmrV1, 'ko-', 'markerfacecolor', 0.6*ones(1,3), ...
      'markersize', 3);
   plot(xRnV1, fRnV1, 'ks--', 'markerfacecolor', 0*ones(1,3), ...
      'markersize', 3);
   tickpref;
   ylabel('Firing Rate (sp/s)');

   subplot(3,3,9);
   hold on;
   plot(xDmrV2, fDmrV2, 'ko-', 'markerfacecolor', 0.6*ones(1,3), ...
      'markersize', 3);
   plot(xRnV2, fRnV2, 'ks--', 'markerfacecolor', 0*ones(1,3), ...
      'markersize', 3);
   tickpref;
   xlabel('Projection (SD)');


   orient portrait;
   print_mfilename(mfilename);
   set(gcf,'position',[171   396   665   514]);

end

cascade;


return;



% 3 mid2 structure, symmetric nonlinearity; structure a little speckled
% 4 similar mid features, just smaller bandwidth for rn
% 5 huge mid2 changes, small mid1 changes
% 9 mid2 fio modulation for dmr, but not for rn
% 11 same as 9






















