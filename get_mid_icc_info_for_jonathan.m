function [ista, iv1, iv2, iv12] = get_mid_icc_info_for_jonathan(projinfo)
%figure7_mid_icc_projinfo_compare_with_a1  Compare ICC/AI sta and mid filter information
%
% plot_mid_icc_projinfo_compare_with_a1(projinfo, infodata)
% -------------------------------------------------------------------
% Mutual information for STA and MID1 and MID2 filters.
%
% The function plots I(STA) vs I(MID1), I(MID1) vs I(MID2), and I(MID1) vs I(MID12).
% 
% projinfo : struct array holding the information data for the ICC. 
% Each element of projinfo holds the data for 1 neuron. 
% Usually saved in a file such as:
%
%    load 2003-11-24-site7-2380um-40db-dmr1-fs18115-mid-fio-projinfo.mat
%
% infodata : struct array holding extrapolated information for the cortex.
%
% In place of infodata you may also use extrap_info, which is located
% in the file: extrapolation_information_database.mat           
%
% caa 1/21/11

close all;

% First get the ICC information
% ---------------------------------------------------------
ista = []; %zeros(1,length(projinfo));
iv1 = []; %zeros(1,length(projinfo));
iv2 = []; %zeros(1,length(projinfo));
iv12 = []; %zeros(1,length(projinfo));

for i = 1:length(projinfo)

   info0 = projinfo(i).info0_extrap_test;
   info1 = projinfo(i).info1_extrap_test;
   info2 = projinfo(i).info2_extrap_test;
   info12 = projinfo(i).info12_extrap_test;

   if ( ~isempty(info0) && ~isempty(info1) && ~isempty(info2) && ~isempty(info12) )

      i0 = info0(end);
      i1 = info1(end);
      i2 = info2(end);
      i12 = info12(end);

      ista = [ista i0];
      iv1 = [iv1 i1];
      iv2 = [iv2 i2];
      iv12 = [iv12 i12];

   end

end % (for i)


% Can't have negative information values, so get rid of these points
index = (ista>0) & (iv1>0) & (iv2>0) & (iv12>0);
ista = ista(index);
iv1 = iv1(index);
iv2 = iv2(index);
iv12 = iv12(index);

% Can't have sta > iv1, or iv1 > iv12, so get rid of these points
index = (ista<=iv1) & (iv1<=iv12) & (iv2<=iv1);
ista_icc = ista(index);
iv1_icc = iv1(index);
iv2_icc = iv2(index);
iv12_icc = iv12(index);


ista_iv1_icc = 100 * ista_icc ./ iv1_icc;
iv1_iv12_icc = 100 * iv1_icc ./ iv12_icc;
iv2_iv1_icc = 100 * iv2_icc ./ iv1_icc;




% =======================================================================
% ==================== Function Declarations ============================
% =======================================================================




function bar_cdf_plots_for_icc_a1_info_compare(ista_iv1_icc, ...
iv1_iv12_icc, ista_iv1_a1, iv1_iv12_a1, position)


% Get indices to the granular and nongranular layer data
index = 1:length(position); % all indicies
indexGran = find(position >= 600 & position <= 1100);
indexNonGran = setdiff(index,indexGran);


% Compare STA and MID1
% =======================================================================

% Get the granular and nongranular data
granData = ista_iv1_a1(indexGran); % cdf for granular layer AI data
nonGranData = ista_iv1_a1(indexNonGran); % cdf for non-granular AI data


% Get cumulative distribution functions
[xIcc, yIcc] = cumdist(ista_iv1_icc); % cdf for icc data
[xGran, yGran] = cumdist(granData);
[xNonGran, yNonGran] = cumdist(nonGranData);


% Compute KS-tests for icc vs. all AI, granular, nongranular data
[h,pval] = kstest2(ista_iv1_icc, ista_iv1_a1); % total icc/a1 comparison
[h,pGran] = kstest2(ista_iv1_icc, granData); % icc/granular layer
[h,pNonGran] = kstest2(ista_iv1_icc, nonGranData); % icc/nongranular layer


% Plot the STA and MID1 comparisons
subplot(2,2,1);
hold on;
plot(xNonGran, yNonGran, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3);
plot(xGran, yGran, '-', 'color', [0.4 0.4 0.4], 'linewidth', 3);
plot(xIcc, yIcc, 'k-', 'linewidth', 3);
xlim([0 100]);
ylim([0 1.05]);
set(gca,'xtick', 0:25:100, 'xticklabel', 0:25:100);
set(gca,'ytick', 0:0.2:1.05, 'yticklabel', 0:0.2:1.05);
set(gca,'tickdir','out', 'ticklength', [0.025 0.025]);
xlabel('STA Sufficiency');
ylabel('Cumulative Proportion');
title(sprintf('p = %.4f, pGran = %.4f', pval, pGran));
legend('NonGran', 'Gran', 'ICC', 0);



% Compare MID1 and MID12 contributions
% =======================================================================

% Get the granular and nongranular data
granData = iv1_iv12_a1(indexGran); % cdf for granular layer AI data
nonGranData = iv1_iv12_a1(indexNonGran); % cdf for non-granular AI data


% Get cumulative distribution functions
[xIcc, yIcc] = cumdist(iv1_iv12_icc); % cdf for icc data
[xGran, yGran] = cumdist(granData);
[xNonGran, yNonGran] = cumdist(nonGranData);


% Compute KS-tests for icc vs. all AI, granular, nongranular data
[h,pval] = kstest2(iv1_iv12_icc, iv1_iv12_a1); % total icc/a1 comparison
[h,pGran] = kstest2(iv1_iv12_icc, granData); % icc/granular layer
[h,pNonGran] = kstest2(iv1_iv12_icc, nonGranData); % icc/nongranular layer



% Plot the MID1 and MID12 comparisons
subplot(2,2,2);
hold on;
plot(xNonGran, yNonGran, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3);
plot(xGran, yGran, '-', 'color', [0.4 0.4 0.4], 'linewidth', 3);
plot(xIcc, yIcc, 'k-', 'linewidth', 3);
xlim([0 100]);
ylim([0 1.05]);
set(gca,'xtick', 0:25:100, 'xticklabel', 0:25:100);
set(gca,'ytick', 0:0.2:1.05, 'yticklabel', 0:0.2:1.05);
set(gca,'tickdir','out', 'ticklength', [0.025 0.025]);
xlabel('MID1 Contribution');
ylabel('Cumulative Proportion');
title(sprintf('p = %.4f, pGran = %.4f', pval, pGran));
legend('NonGran', 'Gran', 'ICC', 0);




% ======================================================================
% =========================== Bar Plots ================================
% ======================================================================

% calculate median ratios for icc
md_ista_iv1_icc = median(ista_iv1_icc);
md_iv1_iv12_icc = median(iv1_iv12_icc);

% calculate Median Absolute Deviation for icc
mad_ista_iv1_icc = mad(ista_iv1_icc(:), 1);
mad_iv1_iv12_icc = mad(iv1_iv12_icc(:), 1);

% calculate mean ratios for icc
mn_ista_iv1_icc = mean(ista_iv1_icc);
mn_iv1_iv12_icc = mean(iv1_iv12_icc);


% calculate median ratios for a1
md_ista_iv1_a1 = median(ista_iv1_a1);
md_iv1_iv12_a1 = median(iv1_iv12_a1);

% calculate Median Absolute Deviation for a1
mad_ista_iv1_a1 = mad(ista_iv1_a1(:),1);
mad_iv1_iv12_a1 = mad(iv1_iv12_a1(:),1);

% calculate mean ratios for a1
mn_ista_iv1_a1 = mean(ista_iv1_a1);
mn_iv1_iv12_a1 = mean(iv1_iv12_a1);


subplot(2,2,3);
ytick = linspace(0, 100, 5);
hb = bar([md_ista_iv1_icc 0 md_ista_iv1_a1]);
set(hb, 'facecolor', [0.6 0.6 0.6])
set(hb, 'edgecolor', [0.6 0.6 0.6])
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'xticklabel', {'ICC', '', 'AI', ''});
xlim([0 4]);
ylim([0 100]);
ylabel('Median STA sufficiency');

subplot(2,2,4);
hb = bar([md_iv1_iv12_icc 0 md_iv1_iv12_a1]);
set(hb, 'edgecolor', [0.6 0.6 0.6])
set(hb, 'facecolor', [0.6 0.6 0.6])
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'xticklabel', {'ICC', '', 'AI', ''});
xlim([0 4]);
ylim([0 100]);
ylabel('Median MID1 contribution');

set(gcf,'position', [1226         485         560         420]);
print_mfilename(mfilename);

return;



function cdf_plots_for_icc_a1_info_compare(ista_iv1_icc, iv1_iv12_icc, ...
ista_iv1_a1, iv1_iv12_a1)

figure;

% Compare STA and MID1
subplot(2,1,1);
[yy, xx, nn] = cdfcalc(ista_iv1_icc);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
icxcdf = [-Inf; xx(nn); Inf];
icycdf = [0; 0; yy(1+nn)];

[yy, xx, nn] = cdfcalc(ista_iv1_a1);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
a1xcdf = [-Inf; xx(nn); Inf];
a1ycdf = [0; 0; yy(1+nn)];

[h,pval] = kstest2(ista_iv1_icc, ista_iv1_a1);

hold on;
plot(icxcdf, icycdf, 'k-', 'linewidth', 3);
plot(a1xcdf, a1ycdf, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3);
xlim([ 0 100 ] );
ylim([0 1.05]);
set(gca,'ytick', 0:0.2:1.05, 'yticklabel', 0:0.2:1.05);
set(gca,'tickdir','out', 'ticklength', [0.025 0.025]);
ylabel('Cumulative Proportion');
title(sprintf('STA and MID1 Comparison; p = %.4f',pval));
legend('ICC', 'AI', 0);



% Compare MID1 contribution
subplot(2,1,2);
[yy, xx, nn] = cdfcalc(iv1_iv12_icc);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
icxcdf = [-Inf; xx(nn); Inf];
icycdf = [0; 0; yy(1+nn)];

[yy, xx, nn] = cdfcalc(iv1_iv12_a1);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
a1xcdf = [-Inf; xx(nn); Inf];
a1ycdf = [0; 0; yy(1+nn)];

[h,pval] = kstest2(iv1_iv12_icc, iv1_iv12_a1);

hold on;
plot(icxcdf, icycdf, 'k-', 'linewidth', 3);
plot(a1xcdf, a1ycdf, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3);
xlim([ 0 100 ] );
ylim([0 1.05]);
set(gca,'ytick', 0:0.2:1.05, 'yticklabel', 0:0.2:1.05);
set(gca,'tickdir','out', 'ticklength', [0.025 0.025]);
ylabel('Cumulative Proportion');
title(sprintf('MID1 Contribution; p = %.4f',pval));

set(gcf,'position', [144   418   341   453]);

print_mfilename(mfilename);


return;



function [infosta_test, infomid1_test, infomid2_test, infomid12_test, ...
position] = get_extrapolated_info_for_icc_a1_compare(infodata)
% plot_infodata - plot extrapolated mutual information for
%    the sta and the two maximally informative dimensions
%
% plot_infodata(infodata)
% ======================================================
%
% infodata - struct array holding the information data. It hold
%    the mutual information values for the training data and the test
%    data for the sta and the two mid's. The data are organized into
%    matrices where the rows correspond to a specific bin size
%    and the column specifies the data fraction size.
%
% The variable infodata is obtained by running the following
% command:
%
%     infodata = plot_proj_info(info);
%
% The input variable infodata is stored in the following file:
%
%     load extrapolation_information_database.mat
%
% As an example, infodata has the following fields:
% 
%             exp: '2003-3-5'
%            site: 32
%            chan: 8
%           model: 32
%           depth: 2357
%        position: 207
%            stim: 'dmr2'
%           atten: 30
%             spl: 75
%              sm: 4
%              tm: 40
%             mdb: 40
%      datalength: [50 60 70 80 90 92.5000 95 97.5000 100]
%     binboundary: 6
%         numbins: [11 13 15 17 19 21 23 25]
%         binsize: [1.2000 1 0.8571 0.7500 0.6667 0.6000 0.5455 0.5000]
%         numreps: 1
%             sta: [1x1 struct]
%            mid1: [1x1 struct]
%            mid2: [1x1 struct]
%           mid12: [1x1 struct]
%
% The sta and mid fields are structs with the following fields:
%
%                  part1_train: [8x9 double]
%                  part2_train: [8x9 double]
%                  part3_train: [8x9 double]
%                  part4_train: [8x9 double]
%           part1_train_extrap: [8x1 double]
%           part2_train_extrap: [8x1 double]
%           part3_train_extrap: [8x1 double]
%           part4_train_extrap: [8x1 double]
%     part1_normr_train_extrap: [8x1 double]
%     part2_normr_train_extrap: [8x1 double]
%     part3_normr_train_extrap: [8x1 double]
%     part4_normr_train_extrap: [8x1 double]
%                   part1_test: [8x9 double]
%                   part2_test: [8x9 double]
%                   part3_test: [8x9 double]
%                   part4_test: [8x9 double]
%            part1_test_extrap: [8x1 double]
%            part2_test_extrap: [8x1 double]
%            part3_test_extrap: [8x1 double]
%            part4_test_extrap: [8x1 double]
%      part1_normr_test_extrap: [8x1 double]
%      part2_normr_test_extrap: [8x1 double]
%      part3_normr_test_extrap: [8x1 double]
%      part4_normr_test_extrap: [8x1 double]
%
%
% infodata is stored in .mat files that are named according to
% the convention
% 
%        site*_extrapolation_information_data.mat
%
%
% where the wildcard is something like 515, 532, 603, etc.
%
%
% The figure(s) can be exported as '.eps' figures using the following
% command:
%
%
%    options = struct('color', 'bw', 'width', 11/2.54, 'height', 9/2.54, 'fontmode', 'fixed', 'fontsize', 8);
%    exportfig(gcf, 'mid_info_compare_scatter.eps', options);
%
%    options = struct('color', 'gray', 'width', 11/2.54, 'height', 9/2.54, 'fontmode', 'fixed', 'fontsize', 8);
%    exportfig(gcf, 'mid_info_compare_hist.eps', options);
%
%
% caa 2/4/07



% The following neurons need to be removed - they are too noisy
% or there are not enough spikes
%
% 515 - none
% 516 - none
% 517 - 12, 13, 16
% 519 - 5, 9, 10, 12
% 532 - 13, 14, 18, 19
% 534 - none
% 537 - 1, 2
% 602 - 1, 2
% 603 - 11
% 604 - 7
% 608 - 1, 12
% 609 - 17
% 611 - none
% 614 - 5
% 616 - none

numbins = infodata(1).numbins;
index1d = find(numbins == 15);
index2d = find(numbins == 15);
% index2d = find(numbins >= 15 & numbins <= 19);

infosta = zeros(1,length(infodata));
infosta_test = zeros(1,length(infodata));

infomid1 = zeros(1,length(infodata));
infomid1_test = zeros(1,length(infodata));

infomid2 = zeros(1,length(infodata));
infomid2_test = zeros(1,length(infodata));

infomid12 = zeros(1,length(infodata));
infomid12_test = zeros(1,length(infodata));

info_mean_pop = [];
info_mean_pop_test = [];
position = [];

for i = 1:length(infodata)

   pos = infodata(i).position;

   % Training data set:

   % Part 1 data:

   sta_info_part1 = infodata(i).sta.part1_train_extrap(index1d);
   mid1_info_part1 = infodata(i).mid1.part1_train_extrap(index1d);
   mid2_info_part1 = infodata(i).mid2.part1_train_extrap(index1d);
   mid12_info_part1 = mean( infodata(i).mid12.part1_train_extrap(index2d) );

   part1_data = [sta_info_part1 mid1_info_part1 mid2_info_part1 mid12_info_part1];
%    [part1_data, err1] = error_check_info_values(sta_info_part1, mid1_info_part1, mid2_info_part1, mid12_info_part1);


   % Part 2 data:

   sta_info_part2 = infodata(i).sta.part2_train_extrap(index1d);
   mid1_info_part2 = infodata(i).mid1.part2_train_extrap(index1d);
   mid2_info_part2 = infodata(i).mid2.part2_train_extrap(index1d);
   mid12_info_part2 = mean( infodata(i).mid12.part2_train_extrap(index2d) );


   part2_data = [sta_info_part2 mid1_info_part2 mid2_info_part2 mid12_info_part2];
%    [part2_data, err2] = error_check_info_values(sta_info_part2, mid1_info_part2, mid2_info_part2, mid12_info_part2);


   % Part 3 data:

   sta_info_part3 = infodata(i).sta.part3_train_extrap(index1d);
   mid1_info_part3 = infodata(i).mid1.part3_train_extrap(index1d);
   mid2_info_part3 = infodata(i).mid2.part3_train_extrap(index1d);
   mid12_info_part3 = mean( infodata(i).mid12.part3_train_extrap(index2d) );


   part3_data = [sta_info_part3 mid1_info_part3 mid2_info_part3 mid12_info_part3];
%    [part3_data, err3] = error_check_info_values(sta_info_part3, mid1_info_part3, mid2_info_part3, mid12_info_part3);


   % Part 4 data:

   sta_info_part4 = infodata(i).sta.part4_train_extrap(index1d);
   mid1_info_part4 = infodata(i).mid1.part4_train_extrap(index1d);
   mid2_info_part4 = infodata(i).mid2.part4_train_extrap(index1d);
   mid12_info_part4 = mean( infodata(i).mid12.part4_train_extrap(index2d) );


   part4_data = [sta_info_part4 mid1_info_part4 mid2_info_part4 mid12_info_part4];
%    [part4_data, err4] = error_check_info_values(sta_info_part4, mid1_info_part4, mid2_info_part4, mid12_info_part4);



   info_values = [part1_data; part2_data; part3_data; part4_data];

   if ( ~isempty(info_values) )
      info_mean = mean(info_values,1);
%       [info_mean, err] = error_check_info_values(info_mean(1), info_mean(2), info_mean(3), info_mean(4));
      info_mean_pop = [info_mean_pop; info_mean];
   end

   clear('part1_data', 'part2_data', 'part3_data', 'part4_data', 'info_values', 'info_mean');

   clear('sta_info_part1', 'mid1_info_part1', 'mid2_info_part1', 'mid12_info_part1');
   clear('sta_info_part2', 'mid1_info_part2', 'mid2_info_part2', 'mid12_info_part2');
   clear('sta_info_part3', 'mid1_info_part3', 'mid2_info_part3', 'mid12_info_part3');
   clear('sta_info_part4', 'mid1_info_part4', 'mid2_info_part4', 'mid12_info_part4');



   % Test data set:

   % Part 1 data:

   sta_info_part1 = infodata(i).sta.part1_test_extrap(index1d);
   mid1_info_part1 = infodata(i).mid1.part1_test_extrap(index1d);
   mid2_info_part1 = infodata(i).mid2.part1_test_extrap(index1d);
   mid12_info_part1 = mean( infodata(i).mid12.part1_test_extrap(index2d) );


%    part1_data = [sta_info_part1 mid1_info_part1 mid2_info_part1 mid12_info_part1];

   [part1_data, err1] = error_check_info_values(sta_info_part1, ...
      mid1_info_part1, mid2_info_part1, mid12_info_part1);


   % Part 2 data:

   sta_info_part2 = infodata(i).sta.part2_test_extrap(index1d);
   mid1_info_part2 = infodata(i).mid1.part2_test_extrap(index1d);
   mid2_info_part2 = infodata(i).mid2.part2_test_extrap(index1d);
   mid12_info_part2 = mean( infodata(i).mid12.part2_test_extrap(index2d) );


%    part2_data = [sta_info_part2 mid1_info_part2 mid2_info_part2 mid12_info_part2];

   [part2_data, err2] = error_check_info_values(sta_info_part2, ...
      mid1_info_part2, mid2_info_part2, mid12_info_part2);


   % Part 3 data:

   sta_info_part3 = infodata(i).sta.part3_test_extrap(index1d);
   mid1_info_part3 = infodata(i).mid1.part3_test_extrap(index1d);
   mid2_info_part3 = infodata(i).mid2.part3_test_extrap(index1d);
   mid12_info_part3 = mean( infodata(i).mid12.part3_test_extrap(index2d) );


%    part3_data = [sta_info_part3 mid1_info_part3 mid2_info_part3 mid12_info_part3];

   [part3_data, err3] = error_check_info_values(sta_info_part3, ...
      mid1_info_part3, mid2_info_part3, mid12_info_part3);


   % Part 4 data:

   sta_info_part4 = infodata(i).sta.part4_test_extrap(index1d);
   mid1_info_part4 = infodata(i).mid1.part4_test_extrap(index1d);
   mid2_info_part4 = infodata(i).mid2.part4_test_extrap(index1d);
   mid12_info_part4 = mean( infodata(i).mid12.part4_test_extrap(index2d) );


%    part4_data = [sta_info_part4 mid1_info_part4 mid2_info_part4 mid12_info_part4];

   [part4_data, err4] = error_check_info_values(sta_info_part4, ...
      mid1_info_part4, mid2_info_part4, mid12_info_part4);


   info_values = [part1_data; part2_data; part3_data; part4_data];


   if ( ~isempty(info_values) )
      info_mean = mean(info_values,1);
      [info_mean, err] = error_check_info_values(info_mean(1), info_mean(2), info_mean(3), info_mean(4));
      info_mean_pop_test = [info_mean_pop_test; info_mean];
      if ( ~err )
         position = [position; pos];
      end
   end

   clear('part1_data', 'part2_data', 'part3_data', 'part4_data', 'info_values', 'info_mean');
   clear('sta_info_part1', 'mid1_info_part1', 'mid2_info_part1', 'mid12_info_part1');
   clear('sta_info_part2', 'mid1_info_part2', 'mid2_info_part2', 'mid12_info_part2');
   clear('sta_info_part3', 'mid1_info_part3', 'mid2_info_part3', 'mid12_info_part3');
   clear('sta_info_part4', 'mid1_info_part4', 'mid2_info_part4', 'mid12_info_part4');

end % (for i)

close('all');



% length(infosta_test)

infosta = info_mean_pop(:,1);
infomid1 = info_mean_pop(:,2);
infomid2 = info_mean_pop(:,3);
infomid12 = info_mean_pop(:,4);

infosta_test = info_mean_pop_test(:,1);
infomid1_test = info_mean_pop_test(:,2);
infomid2_test = info_mean_pop_test(:,3);
infomid12_test = info_mean_pop_test(:,4);

% length(position)
% length(infosta_test)

% [length(infosta) length(infosta_test)]
% pause

return; % end of get_extrapolated_info_for_icc_a1_compare(infodata)




function [data, err] = error_check_info_values(sta_info, mid1_info, mid2_info, mid12_info)

data = [sta_info mid1_info mid2_info mid12_info];

% normr = [sta_info_normr mid1_info_normr mid2_info_normr mid12_info_normr];

if ( sta_info>mid12_info | sta_info>mid1_info | mid1_info>mid12_info | ...
   mid2_info>mid1_info | min(data)<0 | mid12_info<0.01 )
%      mid12_info<0.01 | sum( normr>0.05 ) )
   data = [];
   err = 1;
else
   err = 0;
end


if ( min(data) <= 0 ) % can't have negative or zero extrapolated information values
   data = [];
   err = 1;
else
   err = 0;
end



return;







function calc_stats_for_icc_a1_info_compare(ista_iv1_icc, iv1_iv12_icc, ...
ista_iv1_a1, iv1_iv12_a1)
% ==================================================================================
% ========================= Statistics for ICC and A1 ==============================
% ==================================================================================

% calculate median ratios for icc
md_ista_iv1_icc = median(ista_iv1_icc);
md_iv1_iv12_icc = median(iv1_iv12_icc);

% calculate Median Absolute Deviation for icc
mad_ista_iv1_icc = mad(ista_iv1_icc(:), 1);
mad_iv1_iv12_icc = mad(iv1_iv12_icc(:), 1);

% calculate mean ratios for icc
mn_ista_iv1_icc = mean(ista_iv1_icc);
mn_iv1_iv12_icc = mean(iv1_iv12_icc);


% calculate median ratios for a1
md_ista_iv1_a1 = median(ista_iv1_a1);
md_iv1_iv12_a1 = median(iv1_iv12_a1);

% calculate Median Absolute Deviation for a1
mad_ista_iv1_a1 = mad(ista_iv1_a1(:),1);
mad_iv1_iv12_a1 = mad(iv1_iv12_a1(:),1);

% calculate mean ratios for a1
mn_ista_iv1_a1 = mean(ista_iv1_a1);
mn_iv1_iv12_a1 = mean(iv1_iv12_a1);



fprintf('ICC Medians\n');
fprintf('-----------------------------\n');
fprintf('100 Info(sta) / Info(mid1) = %.3f\n', md_ista_iv1_icc);
fprintf('100 Info(mid1) / Info(mid12) = %.3f\n', md_iv1_iv12_icc);

fprintf('\n');

fprintf('ICC Median absolute deviation\n');
fprintf('-----------------------------\n');
fprintf('100 Info(sta) / Info(mid1) = %.3f\n', mad_ista_iv1_icc);
fprintf('100 Info(mid1) / Info(mid12) = %.3f\n', mad_iv1_iv12_icc);


fprintf('\n');

fprintf('ICC Means\n');
fprintf('-----------------------------\n');
fprintf('100 Info(sta) / Info(mid1) = %.3f\n', mn_ista_iv1_icc);
fprintf('100 Info(mid1) / Info(mid12) = %.3f\n', mn_iv1_iv12_icc);

fprintf('\n\n');

fprintf('A1 Medians\n');
fprintf('-----------------------------\n');
fprintf('100 Info(sta) / Info(mid1) = %.3f\n', md_ista_iv1_a1);
fprintf('100 Info(mid1) / Info(mid12) = %.3f\n', md_iv1_iv12_a1);

fprintf('\n');

fprintf('A1 Median absolute deviation\n');
fprintf('-----------------------------\n');
fprintf('100 Info(sta) / Info(mid1) = %.3f\n', mad_ista_iv1_a1);
fprintf('100 Info(mid1) / Info(mid12) = %.3f\n', mad_iv1_iv12_a1);

fprintf('\n');

fprintf('A1 Means\n');
fprintf('-----------------------------\n');
fprintf('100 Info(sta) / Info(mid1) = %.3f\n', mn_ista_iv1_a1);
fprintf('100 Info(mid1) / Info(mid12) = %.3f\n', mn_iv1_iv12_a1);


%p = ranksum(iv1_a1, ista_a1);
%p = ranksum(iv1_icc, ista_icc);

%p = ranksum(iv1_a1, iv12_a1);
%p = ranksum(iv1_icc, iv12_icc);


p = ranksum(ista_iv1_icc, ista_iv1_a1);
p = ranksum(iv1_iv12_icc, iv1_iv12_a1);

return;




function bar_plots_for_icc_a1_info_compare(ista_iv1_icc, iv1_iv12_icc, ...
ista_iv1_a1, iv1_iv12_a1)


% calculate median ratios for icc
md_ista_iv1_icc = median(ista_iv1_icc);
md_iv1_iv12_icc = median(iv1_iv12_icc);

% calculate Median Absolute Deviation for icc
mad_ista_iv1_icc = mad(ista_iv1_icc(:), 1);
mad_iv1_iv12_icc = mad(iv1_iv12_icc(:), 1);

% calculate mean ratios for icc
mn_ista_iv1_icc = mean(ista_iv1_icc);
mn_iv1_iv12_icc = mean(iv1_iv12_icc);


% calculate median ratios for a1
md_ista_iv1_a1 = median(ista_iv1_a1);
md_iv1_iv12_a1 = median(iv1_iv12_a1);

% calculate Median Absolute Deviation for a1
mad_ista_iv1_a1 = mad(ista_iv1_a1(:),1);
mad_iv1_iv12_a1 = mad(iv1_iv12_a1(:),1);

% calculate mean ratios for a1
mn_ista_iv1_a1 = mean(ista_iv1_a1);
mn_iv1_iv12_a1 = mean(iv1_iv12_a1);

% ==================================================================================
% =========================== Bar Plots ================================
% ==================================================================================

figure;
ytick = linspace(0, 100, 5);
subplot(1,2,1);
hb = bar([md_ista_iv1_icc mad_ista_iv1_icc md_ista_iv1_a1 mad_ista_iv1_a1]);
set(hb, 'facecolor', 1*[0 0 0])
set(hb, 'edgecolor', [0.6 0.6 0.6])
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'xticklabel', {'Midbrain', 'MAD', 'Cortex', 'MAD'});
xlim([0 5]);
ylim([0 100]);
ylabel('STA sufficiency');

subplot(1,2,2);
hb = bar([md_iv1_iv12_icc mad_iv1_iv12_icc md_iv1_iv12_a1 mad_iv1_iv12_a1]);
set(hb, 'facecolor', 1*[0 0 0])
set(hb, 'edgecolor', [0.6 0.6 0.6])
box off;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'xticklabel', {'Midbrain', 'MAD', 'Cortex', 'MAD'});
xlim([0 5]);
ylim([0 100]);
ylabel('MID1 contribution');

return;


function scatter_plots_for_icc_a1_info_compare(ista_icc, iv1_icc, ...
iv2_icc, iv12_icc, ista_iv1_icc, iv2_iv1_icc, iv1_iv12_icc, ista_a1, ...
iv1_a1, iv2_a1, iv12_a1, ista_iv1_a1, iv2_iv1_a1, iv1_iv12_a1)

% ==================================================================================
% =========================== Scatter Plots ================================
% ==================================================================================

markersize = 3;

figure;

subplot(3,2,1);
hold on;
xmax = max([max(iv1_icc) max(ista_icc)]);
xmin = min([min(iv1_icc) min(ista_icc)]);

hp = plot(iv1_a1, ista_a1, 'o');
set(hp, 'markersize', markersize); %, 'markerfacecolor', 'k');
set(hp, 'markerfacecolor', 0.6*[1 1 1]);
set(hp, 'markeredgecolor', 0.6*[1 1 1]);

hp = plot(iv1_icc, ista_icc, 'ko', 'markersize', markersize); %, 'markerfacecolor', 'k');
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
% plot([xmin xmax],[xmin xmax],'k-');


fplot('x', [0.5*xmin 1.75*xmax], 'k-')

xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);
logticks = logspace(log10(0.51*xmin), log10(1.74*xmax), 4);
set(gca,'xscale', 'log');
set(gca, 'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', logticks, 'xticklabel', sprintf('%.2g|', logticks)  );
set(gca,'ytick', logticks, 'yticklabel', sprintf('%.2g|', logticks)  );
set(gca,'xminortick', 'off', 'yminortick', 'off');
% set(gca,'xtick', logticks, 'xticklabel', logticks);
% set(gca,'ytick', logticks, 'yticklabel', logticks);
% set(gca,'minorticks', 'off');
%axis square
box off;
legend('Cortex', 'Midbrain');
ylabel('I(STA) (bits/sp)')
title('STA vs. MID1');



subplot(3,2,2);
hold on;

bins = linspace(0,100,11);
count = histc(ista_iv1_a1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', 1*[1 1 1])
set(hb, 'edgecolor', [0.6 0.6 0.6])

bins = linspace(0,100,11);
count = histc(ista_iv1_icc, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', 0*[1 1 1])
set(hb, 'edgecolor', 0*[1 1 1])

axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:20:100], 'yticklabel', [0:20:100]);
ticks = round( 10 * linspace(0,max(count),5) ) / 10;
ticks = round( 10 * linspace(0,40,5) ) / 10;
set(gca,'ytick', ticks, 'yticklabel', ticks);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
xlabel('100 * I(STA) / I(MID1)');
ylabel('Percent of Neurons');


subplot(3,2,3);
hold on
xmax = max([max(iv1_icc) max(iv2_icc)]);
xmin = min([min(iv1_icc) min(iv2_icc)]);

hp = plot(iv1_a1, iv2_a1, 'o');
set(hp, 'markersize', markersize); %, 'markerfacecolor', 'k');
set(hp, 'markerfacecolor', 0.6*[1 1 1]);
set(hp, 'markeredgecolor', 0.6*[1 1 1]);



hp = plot(iv1_icc, iv2_icc, 'ko', 'markersize', markersize); %, 'markerfacecolor', 'k');
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
fplot('x', [0.5*xmin 1.75*xmax], 'k-')
% axis square
% xlim([0.01 1.75*xmax]);
% ylim([0.01 1.75*xmax]);
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);
logticks = logspace(log10(0.51*xmin), log10(1.74*xmax), 4);
set(gca,'xscale', 'log');
set(gca, 'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', logticks, 'xticklabel', sprintf('%.2g|', logticks)  );
set(gca,'ytick', logticks, 'yticklabel', sprintf('%.2g|', logticks)  );
set(gca,'xminortick', 'off', 'yminortick', 'off');
% set(gca,'xtick', logticks, 'xticklabel', logticks);
% set(gca,'ytick', logticks, 'yticklabel', logticks);
% axis square
box off;
ylabel('I(MID2) (bits/sp)')
title('MID2 vs. MID1');


subplot(3,2,4);
hold on;

bins = linspace(0,100,11);
count = histc(iv2_iv1_a1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', 1*[1 1 1])
set(hb, 'edgecolor', [0.6 0.6 0.6])

bins = linspace(0,100,11);
count = histc(iv2_iv1_icc, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', 0*[1 1 1])
set(hb, 'edgecolor', 0*[1 1 1])

axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:30:100], 'yticklabel', [0:30:100]);

ticks = round( 10 * linspace(0,56,5) ) / 10;
set(gca,'ytick', ticks, 'yticklabel', ticks);

set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;

legend('Cortex', 'Midbrain');

xlabel('100 * I(MID2) / I(MID1)');
ylabel('Percent of Neurons');



subplot(3,2,5);
hold on;
xmax = max([max(iv1_icc) max(iv12_icc)]);
xmin = min([min(iv1_icc) min(iv12_icc)]);

hp = plot(iv1_a1, iv12_a1, 'o');
set(hp, 'markersize', markersize); %, 'markerfacecolor', 'k');
set(hp, 'markerfacecolor', 0.6*[1 1 1]);
set(hp, 'markeredgecolor', 0.6*[1 1 1]);

hp = plot(iv1_icc, iv12_icc, 'ko', 'markersize', markersize); %, 'markerfacecolor', 'k');
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
% plot([xmin xmax],[xmin xmax],'k-');
fplot('x', [0.5*xmin 1.75*xmax], 'k-')
xlim([0.5*xmin 1.75*xmax]);
ylim([0.5*xmin 1.75*xmax]);
logticks = logspace(log10(xmin), log10(1.74*xmax), 4);
set(gca,'xscale', 'log');
set(gca, 'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', logticks, 'xticklabel', sprintf('%.2g|', logticks)  );
set(gca,'ytick', logticks, 'yticklabel', sprintf('%.2g|', logticks)  );
set(gca,'xminortick', 'off', 'yminortick', 'off');
% set(gca,'xtick', logticks, 'xticklabel', logticks);
% set(gca,'ytick', logticks, 'yticklabel', logticks);
%axis square
box off;
xlabel('I(MID1) (bits/sp)')
ylabel('I(MID12) (bits/sp)')
title('MID12 vs. MID1');


subplot(3,2,6);
hold on;

bins = linspace(0,100,11);
count = histc(iv1_iv12_a1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', 1*[1 1 1])
set(hb, 'edgecolor', [0.6 0.6 0.6])


bins = linspace(0,100,11);
count = histc(iv1_iv12_icc, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', 0*[1 1 1])
set(hb, 'edgecolor', 0*[1 1 1])

axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:30:100], 'yticklabel', [0:30:100]);

ticks = round( 10 * linspace(0,60,5) ) / 10;
set(gca,'ytick', ticks, 'yticklabel', ticks);

set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
xlabel('100 * I(MID1) / I(MID12)');
ylabel('Percent of Neurons');

set(gcf,'position', [360   -25   710   723]);

print_mfilename(mfilename);

return;








