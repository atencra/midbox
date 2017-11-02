function [filtstr] = mid_plot_sta_dat_files(exp, site, stim) 
% get_filters_from_dat_files - all filters from mid analysis in one data
%
% [filtstr] = get_filters_from_dat_files(files, prefix, paramfile)
% -----------------------------------------------------------------------------
%
% files : struct array holding files from the mid analysis. Has the 
%              following fields:
%
%     location
%     unit
%     x0
%     nh
%     nv
%     nlags
%     tbins
%     fbins
%     rpsta
%     rpx1pxpxt_sta
%     rpdbest2_v1
%     rpdbest2_v2
%     rpdtest2_v1
%     rpdtest2_v2
%     rpdx1x2px_pxt_2
%
%     location, unit, x0, nh, nb, nlags, tbins, and fbins are scalars.
%     The other fields are 1x4 cell arrays, with each element a string that
%     specifies the *.dat files for the corresponding data type.
%
%     In some cases *sta files do not exist. In this case the current function
%     calculates the STAs.
%
%     files is obtained by:
%
%     files = get_rippledir_detailed_mid_files(location, x0, nh, nlags)
%
% prefix : path to the folder containing the *.dat files specified in
% files. If not included, then the default is 'dat_files\'
%
% paramfile : dmr parameter file. Needed because it holds the frequency and
% time vectors for the filters. If not included then the default is
%
% D:\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat
%
% caa 1/22/10


library('midbox');


% Get the STA files
dfile = dir(['rpsta_' num2str(exp) '_' num2str(site) '_' num2str(stim) '_*.dat']);
if ( ~isempty(dfile) & length(dfile) == 4 )
   files.rpsta{1} = dfile(1).name;
   files.rpsta{2} = dfile(2).name;
   files.rpsta{3} = dfile(3).name;
   files.rpsta{4} = dfile(4).name;
end


% Get the STA input/output function files
dfile = dir(['rpx1pxpxt_sta_' num2str(exp) '_' ...
    num2str(site) '_' ...
    num2str(stim) '_*.dat']);

if ( ~isempty(dfile) & length(dfile) == 4 )
   files.rpx1pxpxt_sta{1} = dfile(1).name;
   files.rpx1pxpxt_sta{2} = dfile(2).name;
   files.rpx1pxpxt_sta{3} = dfile(3).name;
   files.rpx1pxpxt_sta{4} = dfile(4).name;
end


exfile = files.rpsta{1};

index = findstr(exfile, 'x');

Nh = str2double(exfile(index(1)+1:index(2)-1));
nlags = str2double(exfile(index(2)+1:index(2)+2));
Nv = 1;





if ( length(files.rpsta)==4 )
   [v_sta_mean, coeff_sta, projection_sta, mtx_sta] = mid_auditory_filter(files.rpsta, Nv, Nh, nlags);
else
    error('Need 4 filter files.');
end

if ( length(files.rpx1pxpxt_sta) == 4 )
    fio = mid_sta_fio_from_dat_files(files.rpx1pxpxt_sta, coeff_sta);
else
    error('Need 4 nonlinearity files.');
end


vsta{1} = reshape(mtx_sta(:,1),Nh,nlags);
vsta{2} = reshape(mtx_sta(:,2),Nh,nlags);
vsta{3} = reshape(mtx_sta(:,3),Nh,nlags);
vsta{4} = reshape(mtx_sta(:,4),Nh,nlags);

pspkx{1} = fio.pspkx_mtx(:,1);
pspkx{2} = fio.pspkx_mtx(:,2);
pspkx{3} = fio.pspkx_mtx(:,3);
pspkx{4} = fio.pspkx_mtx(:,4);

pspk{1} = fio.pspk_mtx(:,1);
pspk{2} = fio.pspk_mtx(:,2);
pspk{3} = fio.pspk_mtx(:,3);
pspk{4} = fio.pspk_mtx(:,4);

xbins{1} = fio.x_mtx(:,1);
xbins{2} = fio.x_mtx(:,2);
xbins{3} = fio.x_mtx(:,3);
xbins{4} = fio.x_mtx(:,4);



figure;

subplot(2,4,1);
imagesc(vsta{1});
minmin = min(min(vsta{1}));
maxmax = max(max(vsta{1}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
ylabel('Filter');
title('Train Set 1');

subplot(2,4,2);
imagesc(vsta{2});
minmin = min(min(vsta{2}));
maxmax = max(max(vsta{2}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
title('Train Set 2');


subplot(2,4,3);
imagesc(vsta{3});
minmin = min(min(vsta{3}));
maxmax = max(max(vsta{3}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
title('Train Set 3');


subplot(2,4,4);
imagesc(vsta{4});
minmin = min(min(vsta{4}));
maxmax = max(max(vsta{4}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
title('Train Set 4');


% Plot the nonlinearities
% ----------------------------------------------------------
maxmax = max([max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

subplot(2,4,5);
hold on;
plot(xbins{1}, pspkx{1}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{1} pspk{1}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
tickpref;
xlabel('Projection (SD)');
ylabel('P(spk|x)');

subplot(2,4,6);
hold on;
plot(xbins{2}, pspkx{2}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{2} pspk{2}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
tickpref;
xlabel('Projection (SD)');

subplot(2,4,7);
hold on;
plot(xbins{3}, pspkx{3}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{3} pspk{3}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
tickpref;
xlabel('Projection (SD)');

subplot(2,4,8);
hold on;
plot(xbins{4}, pspkx{4}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{4} pspk{4}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
tickpref;
xlabel('Projection (SD)');

set(gcf,'position', [370 556 1056 363]);



return;









