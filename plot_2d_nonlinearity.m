function plot_2d_nonlinearity(xbins, pspk, px1x2, px1x2spk, pspkx1x2)
%plot_2d_nonlinearity - Display combined MID1 and MID2 nonlinearity
%
% plot_2d_nonlinearity(xbins, pspk, px1x2, px1x2spk, pspkx1x2)
% ------------------------------------------------------------
%
% xbins : bin centers at which the probability values were calculated
%
% pspk : 1x4 cell array. Probability of a spike for each training set.
%
% px1x2 : 1x4 cell array. Prior projection probability distribution.
%
% px1x2spk : 1x4 cell array. Probability of a projection for each 
% training set. 
%
% pspkx1x2 : 1x4 cell array. Probability of a spike given a projection 
% for each training set. 
%
% caa 3/17/09


% Plot the MID1 and MID2 nonlinearity
%------------------------------------
figure;

subplot(2,2,1);
% h = imagesc(xbins, xbins, pspkx1x2{1});
h = imagesc(xbins, xbins, log10(pspkx1x2{1}+0.0001));
set(gca,'clim', [-3 0]);
set(gca,'tickdir', 'out');
set(gca,'xtick', xbins);
set(gca,'ytick', xbins);
xlabel('X1 (SD)','fontsize',8);
ylabel('X2 (SD)', 'fontsize',8);
axis('xy');
%axis('square');
% hc = colorbar('eastoutside');
% set(gca,'fontsize', 8);
% set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
% set(hc,'tickdir', 'out'); 
title('Training Set 1', 'fontsize', 8);


subplot(2,2,2);
% h = imagesc(xbins, xbins, pspkx1x2{2});
h = imagesc(xbins, xbins, log10(pspkx1x2{2}+0.0001));
set(gca,'clim', [-3 0]);
set(gca,'tickdir', 'out');
set(gca,'xtick', xbins);
set(gca,'ytick', xbins);
xlabel('X1 (SD)','fontsize',8);
ylabel('X2 (SD)', 'fontsize',8);
axis('xy');
%axis('square');
title('Training Set 2', 'fontsize', 8);


subplot(2,2,3);
% h = imagesc(xbins, xbins, pspkx1x2{3});
h = imagesc(xbins, xbins, log10(pspkx1x2{3}+0.0001));
set(gca,'clim', [-3 0]);
set(gca,'tickdir', 'out');
set(gca,'xtick', xbins);
set(gca,'ytick', xbins);
xlabel('X1 (SD)','fontsize',8);
ylabel('X2 (SD)', 'fontsize',8);
axis('xy');
%axis('square');
title('Training Set 3', 'fontsize', 8);


subplot(2,2,4);
% h = imagesc(xbins, xbins, pspkx1x2{4});
h = imagesc(xbins, xbins, log10(pspkx1x2{4}+0.0001));
set(gca,'clim', [-3 0]);
set(gca,'tickdir', 'out');
set(gca,'xtick', xbins);
set(gca,'ytick', xbins);
xlabel('X1 (SD)','fontsize',8);
ylabel('X2 (SD)', 'fontsize',8);
axis('xy');
%axis('square');
title('Training Set 4', 'fontsize', 8);





