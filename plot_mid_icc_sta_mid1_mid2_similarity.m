function [si_sta_mid1, si_mid1_mid2] = plot_mid_icc_sta_mid1_mid2_similarity(filtstr)
%plot_mid_icc_sta_mid1_mid2_similarity - STA, MID1,2 filter correlations
%
% plot_mid_icc_sta_mid1_similarity(filtstr)
% ---------------------------------------------------------------------
%
% filtstr : struct array holding the filters. 
%
% Two histograms are created. One for the correlation between the
% STA and MID1, and another for the correlation between MID1 and MID2.
%
%   caa 1/28/10


si_sta_mid1 = [];
si_mid1_mid2 = [];

for i = 1:length(filtstr)

   sta = filtstr(i).v_sta;
   sta = sta(:);

   mid1 = filtstr(i).v1;
   mid1 = mid1(:);

   mid2 = filtstr(i).v2;
   mid2 = mid2(:);

   if ( ~isempty(sta) && ~isempty(mid1) && ~isempty(mid2) )

      r = corrcoef(sta, mid1);
      si_sta_mid1 = [si_sta_mid1 r(1,2)];

      r = corrcoef(mid1, mid2);
      si_mid1_mid2 = [si_mid1_mid2 r(1,2)];

   end

end % (for i)



close all;

figure;


subplot(2,1,1);
bins = 0:0.05:1;
xtick = 0:0.2:1;
ytick = 0:0.1:0.4;
count = histc(si_sta_mid1, bins);
count = count ./ sum(count);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
set(gca,'fontsize', 8);
xlabel('STA/MID1 Similarity', 'fontsize', 10);
ylabel('#Units', 'fontsize', 10);


subplot(2,1,2);
bins = -0.5:0.05:0.5;
xtick = -0.5:0.25:0.5;
ytick = 0:0.1:0.4;
count = histc(si_mid1_mid2, bins);
count = count ./ sum(count);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.6 0.6 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
set(gca,'fontsize', 8);
xlabel('MID1/MID2 Similarity', 'fontsize', 10);
ylabel('#Units', 'fontsize', 10);
suptitle(sprintf('Correlation between STA/MID1\n and MID1/MID2 Filters'));

set(gcf,'position', [753   226   411   668]);
print_mfilename(mfilename);


return;








