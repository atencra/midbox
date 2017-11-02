function plot_mid_icc_fio_similarity(fio)
%plot_fio - plots nonlinearity and probability distributions.
%
% plot_mid_icc_fio_similarity(fio)
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
% opts = struct('color', 'gray', 'FontMode','fixed','FontSize',8,'width',4.5);
%
% caa 2/1/10



close all;


si = zeros(1,length(fio));

for i = 1:length(fio)

   pspk = fio(i).pspk;
   pspkx0 = fio(i).pspkx0;
   pspkx1 = fio(i).pspkx1;

   pspkx0_mtx = [];
   pspkx1_mtx = [];
   for j = 1:length(pspkx0)
      pspkx0_mtx = [pspkx0_mtx pspkx0{j}];
      pspkx1_mtx = [pspkx1_mtx pspkx1{j}];
   end


   pspkx0_mn = nanmean(pspkx0_mtx,2);
   pspkx1_mn = nanmean(pspkx1_mtx,2);

   index = ~isnan(pspkx0_mn) & ~isnan(pspkx1_mn);

   r = corrcoef(pspkx0_mn(index), pspkx1_mn(index));

   si(i) = r(1,2);

end % (for i)


figure;


median(si)
mad(si,1)

bins = 0:0.05:1;

xtick = 0:0.2:1;

bins = linspace(0,1.00,11);
count = histc(si, bins);
total_count = sum(count);
count = 100 * count ./ total_count;

hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
set(gca,'fontsize', 8);
xlabel('STA/MID1 Similarity', 'fontsize', 10);
ylabel('Percent of Neurons');
suptitle('Correlation between STA and MID1 Nonlinearities');

print_mfilename(mfilename);


return;





