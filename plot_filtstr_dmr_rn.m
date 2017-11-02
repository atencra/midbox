function plot_filtstr_dmr_rn(dmrFiltStr, rnFiltStr)
% plot_filtstr_dmr_rn DMR/RN STA, MID1, and MID2 filters
%
%    plot_filtstr_dmr_rn(dmrFiltStr, rnFiltStr) displays the STA,
%    MID1, and MID2 for dynamic moving ripple and ripple noise
%    stimuli. Each figure holds the results for one neuron.
%    Each column in each figure holds the data from DMR or RN.
% 
%    dmrFiltStr is a struct array holding MID data for dmr stimulation.
%    rnFiltStr is a struct array that holds data for ripple noise.
% 
%    caa 8/1/2012
% 
%    plot_filtstr_dmr_rn(dmrFiltStr, rnFiltStr)

mapscheme = 'rdbu';

if ( nargin ~= 2 )
   error('You need 2 input args.');
end

if ( length(dmrFiltStr) ~= length(rnFiltStr) )
   error('Input args must have the same length.');
end


time = 1000 .* dmrFiltStr(1).time; % change to ms
freq = dmrFiltStr(1).freq / 1000;

ytick = [1 floor(length(freq)/2) length(freq)];
yticklabel = round(freq(ytick)*10)/10;

xtick = [1 floor(length(time)/2)+1 length(time)];
xticklabel = round( time(xtick) );

close all;

for i = 1:length(dmrFiltStr)

   exp = dmrFiltStr(i).exp;
   site = dmrFiltStr(i).site;

   dmrunit = dmrFiltStr(i).unit;
   rnunit = rnFiltStr(i).unit;

   dmrSTA = dmrFiltStr(i).v_sta;
   dmrSTA = fliplr(dmrSTA);
   dmrV1 = dmrFiltStr(i).v1;
   dmrV1 = fliplr(dmrV1);
   dmrV2 = dmrFiltStr(i).v2;
   dmrV2 = fliplr(dmrV2);

   rnSTA = rnFiltStr(i).v_sta;
   rnSTA = fliplr(rnSTA);
   rnV1 = rnFiltStr(i).v1;
   rnV1 = fliplr(rnV1);
   rnV2 = rnFiltStr(i).v2;
   rnV2 = fliplr(rnV2);

   figure;

   subplot(3,2,1);
   imagesc( dmrSTA );
%    imagesc(time, freq, dmrSTA );
   axis xy;
   minmin = min(min(dmrSTA));
   maxmax = max(max(dmrSTA));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   tickpref;
   cmap = cschemes(mapscheme,15);
   h = title(sprintf('DMR\n%s site %.0f unit %.0f', exp, site, dmrunit));
   set(h, 'fontweight', 'bold');
   h = ylabel('STA');
   set(h, 'rotation', 0, 'fontweight', 'bold');


   subplot(3,2,3);
   imagesc( dmrV1 );
   axis xy;
   minmin = min(min(dmrV1));
   maxmax = max(max(dmrV1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   tickpref;
   cmap = cschemes(mapscheme,15);
   colormap(cmap);
   h = ylabel('MID1');
   set(h, 'rotation', 0, 'fontweight', 'bold');


   subplot(3,2,5);
   imagesc( dmrV2 );
   axis xy;
   minmin = min(min(dmrV2));
   maxmax = max(max(dmrV2));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   tickpref;
   cmap = cschemes(mapscheme,15);
   colormap(cmap);
   h = ylabel('MID2');
   set(h, 'rotation', 0, 'fontweight', 'bold');
   xlabel('Time (ms)');



   subplot(3,2,2);
   imagesc( rnSTA );
%    imagesc(time, freq, rnSTA );
   axis xy;
   minmin = min(min(rnSTA));
   maxmax = max(max(rnSTA));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   tickpref;
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   h = title(sprintf('RN\n%s site %.0f unit %.0f', exp, site, rnunit));
   set(h, 'fontweight', 'bold');


   subplot(3,2,4);
   imagesc( rnV1 );
   axis xy;
   minmin = min(min(rnV1));
   maxmax = max(max(rnV1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   tickpref;
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   colormap(cmap);

   subplot(3,2,6);
   imagesc( rnV2 );
   axis xy;
   minmin = min(min(rnV2));
   maxmax = max(max(rnV2));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   tickpref;
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
   colormap(cmap);
   xlabel('Time (ms)');

   orient tall;
   print_mfilename(mfilename);
   set(gcf,'position',[50 100 450 750]);

end



return;


























