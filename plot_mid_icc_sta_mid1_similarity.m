function plot_mid_icc_sta_mid1_similarity(filtstr)
%plot_mid_icc_sta_mid1_similarity - Plots tuning curve data
%   obtained using the Michigan 16 channel silicon probe.
%
% plot_mid_icc_sta_mid1_similarity(filtstr)
% ---------------------------------------------------------------------
%
% filtstr is saved in files ending with *-mid-filters.mat, for example:
%
%     2003-11-24-site7-2380um-40db-dmr1-fs18115-mid-filters.mat       
%
%   caa 1/28/10


si_sta_mid1 = []; %zeros(1,length(filtstr));

for i = 1:length(filtstr)

   sta = filtstr(i).v_sta;
   sta = sta(:);

   mid1 = filtstr(i).v1;
   mid1 = mid1(:);

   if ( ~isempty(sta) && ~isempty(mid1) )

      r = corrcoef(sta, mid1);

      si_sta_mid1 = [si_sta_mid1 r(1,2)];

   end

end % (for i)



close all;

figure;

bins = 0:0.05:1;
xtick = 0:0.2:1;

% subplot(4,2,2);
count = histc(si_sta_mid1, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
set(gca,'fontsize', 8);
xlabel('STA/MID1 Similarity', 'fontsize', 10);
ylabel('#Units', 'fontsize', 10);
suptitle('Correlation between STA and MID1 Filters');

print_mfilename(mfilename);

return;


