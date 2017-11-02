function get_sta_mid1_mid2_prob_data(mid)
% get_sta_mid1_mid2_prob_data(mid) - 
%
% For illustrations in adobe the following is useful:
%
% exportfig(gcf,'sta_v1_v2_fx.eps', 'color', 'rgb', 'height', 8, 'width', 11, 'fontmode', 'fixed', 'fontsize', 6);
%
% caa 6/01/06


sta_peak_prob = [];
mid1_peak_prob = [];

for i = 1:length(mid)

   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
   sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   pspkx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxspk_sta = pspkx_sta ./ pspike .* px_sta;

   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
   v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   pspkx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxspk_v1 = pspkx_v1 ./ pspike .* px_v1;

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
   v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   pspkx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxspk_v2 = pspkx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   pspkx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   sta_peak_prob = [sta_peak_prob max(pspkx_sta)];
   mid1_peak_prob = [mid1_peak_prob max(pspkx_v1)];

end % (for i)

% mean(sta_peak_prob)
% std(sta_peak_prob)
% length(sta_peak_prob)

fprintf('STA:  mean peak prob=%.3f, std = %.3f, N=%.0f\n',  mean(sta_peak_prob), ...
std(sta_peak_prob), length(sta_peak_prob) );

% mean(mid1_peak_prob)
% std(mid1_peak_prob)
% length(mid1_peak_prob)

fprintf('MID1: mean peak prob=%.3f, std = %.3f, N=%.0f\n',  mean(mid1_peak_prob), ...
std(mid1_peak_prob), length(mid1_peak_prob) );

[ha, pval, ci, stats] = ttest2(sta_peak_prob, mid1_peak_prob, 0.05, 0);
ha
pval

peak_prob = [sta_peak_prob mid1_peak_prob];
fprintf('\nCombined: mean peak prob=%.3f, std = %.3f, N=%.0f\n',  mean(peak_prob), ...
std(peak_prob), length(peak_prob) );






