function [filtstr] = mid_plot_mid1_dat_files(exp, site, stim, optrun) 
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

narginchk(3,4);

if nargin == 3
    optrun = 2;
end

library('midbox');
library('mid_dat_box');

% Get MID1 filter files
fileform = sprintf('rpdtest%.0f_v1_%.0f_%.0f_%.0f_*.dat', optrun, exp, site, stim);

dfile = dir(fileform);
if ( length(dfile) == 4 )
   files.mid1{1} = dfile(1).name;
   files.mid1{2} = dfile(2).name;
   files.mid1{3} = dfile(3).name;
   files.mid1{4} = dfile(4).name;
else
   error('Need 4 filter files.');
end



% Get MID1 nonlinearity files
fileform = sprintf('rpdx1x2px_pxt_1_%.0f_%.0f_%.0f_1x25x20_1_*.dat', exp, site, stim);

% Get the STA input/output function files
dfile = dir([fileform]);
if ( length(dfile) == 4 )
   files.fio{1} = dfile(1).name;
   files.fio{2} = dfile(2).name;
   files.fio{3} = dfile(3).name;
   files.fio{4} = dfile(4).name;
else
    error('Need 4 nonlinearity files.');
end


exfile = files.mid1{1};
index = findstr(exfile, 'x');
Nh = str2double(exfile(index(end-1)+1:index(end)-1));
nlags = str2double(exfile(index(end)+1:index(end)+2));
Nv = 1;


[mid1_mean, coeff_mid1, projection_mid1, mtx_mid1] = ...
    mid_auditory_filter(files.mid1, Nv, Nh, nlags);

fio = mid_mid1_fio_from_dat_files(files.fio, coeff_mid1);


vmid1{1} = reshape(mtx_mid1(:,1),Nh,nlags);
vmid1{2} = reshape(mtx_mid1(:,2),Nh,nlags);
vmid1{3} = reshape(mtx_mid1(:,3),Nh,nlags);
vmid1{4} = reshape(mtx_mid1(:,4),Nh,nlags);

pspkx{1} = fio.pspkx1_mtx(:,1);
pspkx{2} = fio.pspkx1_mtx(:,2);
pspkx{3} = fio.pspkx1_mtx(:,3);
pspkx{4} = fio.pspkx1_mtx(:,4);

pspk{1} = fio.pspk_mtx(:,1);
pspk{2} = fio.pspk_mtx(:,2);
pspk{3} = fio.pspk_mtx(:,3);
pspk{4} = fio.pspk_mtx(:,4);

xbins{1} = fio.x1_mtx(:,1);
xbins{2} = fio.x1_mtx(:,2);
xbins{3} = fio.x1_mtx(:,3);
xbins{4} = fio.x1_mtx(:,4);



figure;

subplot(2,4,1);
imagesc(vmid1{1});
minmin = min(min(vmid1{1}));
maxmax = max(max(vmid1{1}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
ylabel('Filter');
title('Train Set 1');


subplot(2,4,2);
imagesc(vmid1{2});
minmin = min(min(vmid1{2}));
maxmax = max(max(vmid1{2}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
title('Train Set 2');


subplot(2,4,3);
imagesc(vmid1{3});
minmin = min(min(vmid1{3}));
maxmax = max(max(vmid1{3}));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
cmap = brewmaps('rdbu',21);
colormap(cmap);
title('Train Set 3');


subplot(2,4,4);
imagesc(vmid1{4});
minmin = min(min(vmid1{4}));
maxmax = max(max(vmid1{4}));
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










