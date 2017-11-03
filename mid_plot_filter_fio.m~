function mid_plot_filter_fio(data)
%mid_plot_filter_fio Disply MID Filters/Nonlinearities
%    
%    mid_plot_filter_fio(data) plots the filters and nonlinearities stored
%    in the struct 'data'. data is stored in files having the form *-filter-fio.mat
%    and is obtained with mid_dir_file_struct_to_filters.m
%




iskfile = data.iskfile;
nspk = sum(data.locator);
nf_filter = data.nf_filter;
nt_filter = data.nt_filter;


sta_filt = data.filter_mean_sta;
sta_fio_x = data.fio_sta.x_mean;
sta_fio_ior = data.fio_sta.ior_mean;
sta_fio_ior_std = data.fio_sta.ior_std;

pspk = data.fio_sta.pspk_mean;

mid1_filt = data.filter_mean_test2_v1;
mid1_fio_x = data.fio_mid1.x1_mean;
mid1_fio_ior = data.fio_mid1.pspkx1_mean;
mid1_fio_ior_std = data.fio_mid1.pspkx1_std;


mid2_filt = data.filter_mean_test2_v2;
mid2_fio_x = data.fio_mid2.x2_mean;
mid2_fio_ior = data.fio_mid2.pspkx2_mean;
mid2_fio_ior_std = data.fio_mid2.pspkx2_std;



subplot(3,2,1);
imagesc(sta_filt);
tickpref;
title(sprintf('STA; nspk=%.0f, ndim=%.0f, nspk/ndim=%.2f', ...
    nspk, nf_filter*nt_filter, nspk / (nf_filter*nt_filter)));

subplot(3,2,2);
hold on;
plot(sta_fio_x, sta_fio_ior, 'ko-', 'markerfacecolor', 'k');
xlim([1.1*min(sta_fio_x) 1.1*max(sta_fio_x)]);
ylim([0 1.1*max(sta_fio_ior)]);
plot(xlim, [pspk pspk], 'k--');
tickpref;
title('STA Nonlinearity');
    

subplot(3,2,3);
imagesc(mid1_filt);
tickpref;
title('MID1');

subplot(3,2,4);
hold on;
plot(mid1_fio_x, mid1_fio_ior, 'ko-', 'markerfacecolor', 'k');
xlim([1.1*min(mid1_fio_x) 1.1*max(mid1_fio_x)]);
ylim([0 1.1*max(mid1_fio_ior)]);
plot(xlim, [pspk pspk], 'k--');
tickpref;
title('MID1 Nonlinearity');


subplot(3,2,5);
imagesc(mid2_filt);
tickpref;
title('MID2');

subplot(3,2,6);
hold on;
plot(mid2_fio_x, mid2_fio_ior, 'ko-', 'markerfacecolor', 'k');
xlim([1.1*min(mid2_fio_x) 1.1*max(mid2_fio_x)]);
ylim([0 1.1*max(mid2_fio_ior)]);
plot(xlim, [pspk pspk], 'k--');
tickpref;
title('MID2 Nonlinearity');

suptitle(iskfile);

ss = get(0,'screensize');
set(gcf,'position', [ss(3)-1.05*560 ss(4)-1.2*720 560 720]);

print_mfilename(mfilename);

return;





