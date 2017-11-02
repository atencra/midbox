function plot_fio_single(fio, sr)
%plot_fio - plots nonlinearity and probability distributions for fio data.
%
% plot_fio_single(fio, sr)
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

if ( nargin == 1 )
   sr = 1000;
end


xbins = fio(1).x0bins;
pspk = fio(1).pspk;
pspkx = fio(1).pspkx0;

pspkx_mtx = [];
pspk_mtx = [];
for i = 1:length(pspkx)
   pspkx_mtx = [pspkx_mtx pspkx{i}];
   pspk_mtx = [pspk_mtx pspk{i}];
end

pspkx_mn = nanmean(pspkx_mtx,2);
pspkx_std = nanstd(pspkx_mtx,0,2);
pspk_mn = nanmean(pspk_mtx,2);





location = fio(1).location;
unit = fio(1).unit;

xbins = fio(1).x0bins;
pspk = fio(1).pspk;
pspkx0 = fio(1).pspkx0;
pspkx1 = fio(1).pspkx1;

pspkx0_mtx = [];
pspkx1_mtx = [];
pspk_mtx = [];
for j = 1:length(pspkx0)
   pspkx0_mtx = [pspkx0_mtx pspkx0{j}];
   pspkx1_mtx = [pspkx1_mtx pspkx1{j}];
   pspk_mtx = [pspk_mtx pspk{j}];
end

pspk_mn = nanmean(pspk_mtx,2);
pspkx0_mn = nanmean(pspkx0_mtx,2);
pspkx0_std = nanstd(pspkx0_mtx,0,2);
pspkx1_mn = nanmean(pspkx1_mtx,2);
pspkx1_std = nanstd(pspkx1_mtx,0,2);

max0 = sr * max( [ pspkx0_mn + pspkx0_std] );
max1 = sr * max( [ pspkx1_mn + pspkx1_std] );
maxmax = max([max0 max1]);

ytick = [0 round(maxmax/2) 2*round(maxmax/2)];



close all;

figure;

maxmax = sr * nanmax([ max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

subplot(5,2,1);
hold on;
plot(xbins, sr * pspkx{1}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], sr * [pspk{1} pspk{1}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% ylabel('FR (sp/s)');
ylabel('Train Set 1');
set(gca,'xticklabel', []);


subplot(5,2,2);
hold on;
plot(xbins, sr * pspkx{2}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], sr * [pspk{2} pspk{2}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
ylabel('Train Set 2');
set(gca,'xticklabel', []);


subplot(5,2,3);
hold on;
plot(xbins, sr * pspkx{3}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], sr * [pspk{3} pspk{3}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% ylabel('FR (sp/s)');
ylabel('Train Set 3');
set(gca,'xticklabel', []);


subplot(5,2,4);
hold on;
plot(xbins, sr * pspkx{4}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], sr * [pspk{4} pspk{4}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
ylabel('Train Set 4');
set(gca,'xticklabel', []);


subplot(5,2,5);
hold on;
plot(xbins, sr * pspkx_mn, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
% xlabel('Projection (SD)');
% ylabel('FR (sp/s)');
ylabel('Mean');
set(gca,'xticklabel', []);



subplot(5,2,6);
hold on;
plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
hp = ebar(xbins, sr * pspkx_mn, sr * pspkx_std);
xlim([-8 8]);
maxmax = sr * max( [ pspkx_mn + pspkx_std ] );
ylim([0 1.05 * maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
% xlabel('Projection (SD)');
% title('Nonlinearity');
set(gca,'xticklabel', []);



tempx = xbins(1:end-1);

temp0_mn = sr * pspkx0_mn(1:end-1);
temp0_std = sr * pspkx0_std(1:end-1);

temp1_mn = sr * pspkx1_mn(1:end-1);
temp1_std = sr * pspkx1_std(1:end-1);

maxmax = nanmax([ temp0_mn+temp0_std;  temp1_mn+temp1_std]);

ytick = [0 round(maxmax/2) 2*round(maxmax/2)];



subplot(5,2,7);
hold on;
plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
hp = ebar(tempx, temp0_mn, temp0_std);
% hp = ebar(xbins, sr * pspkx0_mn, sr * pspkx0_std);
set(hp, 'markersize', 2);
xlim([-8 8]);
ylim([0 1.05 * maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
set(gca,'ytick', ytick', 'yticklabel', ytick);
% xlabel('Projection (SD)');
ylabel(sprintf('%.0f - %.0f', location, unit));
% title('STA Nonlinearity');
set(gca,'xticklabel', []);


subplot(5,2,8);
hold on;
plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
hp = ebar(tempx, temp1_mn, temp1_std);
% hp = ebar(xbins, sr * pspkx1_mn, sr * pspkx1_std);
set(hp, 'markersize', 2);
xlim([-8 8]);
ylim([0 1.05 * maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
set(gca,'ytick', ytick', 'yticklabel', ytick);
% xlabel('Projection (SD)');
% title('MID1 Nonlinearity');
set(gca,'xticklabel', []);



subplot(5,2,9);
hold on;
plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
hp = ebar(tempx, temp0_mn, temp0_std);
% hp = ebar(xbins, sr * pspkx0_mn, sr * pspkx0_std);
set(hp, 'markersize', 2);
xlim([-8 8]);
ylim([0 1.05 * maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
set(gca,'ytick', ytick', 'yticklabel', ytick);
xlabel('Projection (SD)');
ylabel(sprintf('%.0f - %.0f', location, unit));
% title('STA Nonlinearity');


subplot(5,2,10);
hold on;
plot([-7 7], sr * [pspk_mn pspk_mn], 'k--');
hp = ebar(tempx, temp1_mn, temp1_std);
% hp = ebar(xbins, sr * pspkx1_mn, sr * pspkx1_std);
set(hp, 'markersize', 2);
xlim([-8 8]);
ylim([0 1.05 * maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.0255]);
set(gca,'ytick', ytick', 'yticklabel', ytick);
xlabel('Projection (SD)');
% title('MID1 Nonlinearity');



print_mfilename(mfilename);
set(gcf, 'position', [908   403   361   536]);


return;












