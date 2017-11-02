function [infosta_test, infomid1_test, infomid2_test, infomid12_test] = plot_infodata(infodata)
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

for i = 1:length(infodata)


   % Training data set:

   % Part 1 data:

   sta_info_part1 = infodata(i).sta.part1_train_extrap(index1d);
   mid1_info_part1 = infodata(i).mid1.part1_train_extrap(index1d);
   mid2_info_part1 = infodata(i).mid2.part1_train_extrap(index1d);
   mid12_info_part1 = mean( infodata(i).mid12.part1_train_extrap(index2d) );

   sta_info_part1_normr = infodata(i).sta.part1_normr_train_extrap(index1d);
   mid1_info_part1_normr = infodata(i).mid1.part1_normr_train_extrap(index1d);
   mid2_info_part1_normr = infodata(i).mid2.part1_normr_train_extrap(index1d);
   mid12_info_part1_normr = mean( infodata(i).mid12.part1_normr_train_extrap(index2d) );

   part1_data = [sta_info_part1 mid1_info_part1 mid2_info_part1 mid12_info_part1];
%    [part1_data, err1] = error_check_info_values(sta_info_part1, mid1_info_part1, mid2_info_part1, mid12_info_part1);


   % Part 2 data:

   sta_info_part2 = infodata(i).sta.part2_train_extrap(index1d);
   mid1_info_part2 = infodata(i).mid1.part2_train_extrap(index1d);
   mid2_info_part2 = infodata(i).mid2.part2_train_extrap(index1d);
   mid12_info_part2 = mean( infodata(i).mid12.part2_train_extrap(index2d) );

   sta_info_part2_normr = infodata(i).sta.part2_normr_train_extrap(index1d);
   mid1_info_part2_normr = infodata(i).mid1.part2_normr_train_extrap(index1d);
   mid2_info_part2_normr = infodata(i).mid2.part2_normr_train_extrap(index1d);
   mid12_info_part2_normr = mean( infodata(i).mid12.part2_normr_train_extrap(index2d) );


   part2_data = [sta_info_part2 mid1_info_part2 mid2_info_part2 mid12_info_part2];
%    [part2_data, err2] = error_check_info_values(sta_info_part2, mid1_info_part2, mid2_info_part2, mid12_info_part2);


   % Part 3 data:

   sta_info_part3 = infodata(i).sta.part3_train_extrap(index1d);
   mid1_info_part3 = infodata(i).mid1.part3_train_extrap(index1d);
   mid2_info_part3 = infodata(i).mid2.part3_train_extrap(index1d);
   mid12_info_part3 = mean( infodata(i).mid12.part3_train_extrap(index2d) );

   sta_info_part3_normr = infodata(i).sta.part3_normr_train_extrap(index1d);
   mid1_info_part3_normr = infodata(i).mid1.part3_normr_train_extrap(index1d);
   mid2_info_part3_normr = infodata(i).mid2.part3_normr_train_extrap(index1d);
   mid12_info_part3_normr = mean( infodata(i).mid12.part3_normr_train_extrap(index2d) );

   part3_data = [sta_info_part3 mid1_info_part3 mid2_info_part3 mid12_info_part3];
%    [part3_data, err3] = error_check_info_values(sta_info_part3, mid1_info_part3, mid2_info_part3, mid12_info_part3);


   % Part 4 data:

   sta_info_part4 = infodata(i).sta.part4_train_extrap(index1d);
   mid1_info_part4 = infodata(i).mid1.part4_train_extrap(index1d);
   mid2_info_part4 = infodata(i).mid2.part4_train_extrap(index1d);
   mid12_info_part4 = mean( infodata(i).mid12.part4_train_extrap(index2d) );

   sta_info_part4_normr = infodata(i).sta.part4_normr_train_extrap(index1d);
   mid1_info_part4_normr = infodata(i).mid1.part4_normr_train_extrap(index1d);
   mid2_info_part4_normr = infodata(i).mid2.part4_normr_train_extrap(index1d);
   mid12_info_part4_normr = mean( infodata(i).mid12.part4_normr_train_extrap(index2d) );

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

   clear('sta_info_part1_normr', 'mid1_info_part1_normr', 'mid2_info_part1_normr', 'mid12_info_part1_normr');
   clear('sta_info_part2_normr', 'mid1_info_part2_normr', 'mid2_info_part2_normr', 'mid12_info_part2_normr');
   clear('sta_info_part3_normr', 'mid1_info_part3_normr', 'mid2_info_part3_normr', 'mid12_info_part3_normr');
   clear('sta_info_part4_normr', 'mid1_info_part4_normr', 'mid2_info_part4_normr', 'mid12_info_part4_normr');



   % Test data set:

   % Part 1 data:

   sta_info_part1 = infodata(i).sta.part1_test_extrap(index1d);
   mid1_info_part1 = infodata(i).mid1.part1_test_extrap(index1d);
   mid2_info_part1 = infodata(i).mid2.part1_test_extrap(index1d);
   mid12_info_part1 = mean( infodata(i).mid12.part1_test_extrap(index2d) );

   sta_info_part1_normr = infodata(i).sta.part1_normr_test_extrap(index1d);
   mid1_info_part1_normr = infodata(i).mid1.part1_normr_test_extrap(index1d);
   mid2_info_part1_normr = infodata(i).mid2.part1_normr_test_extrap(index1d);
   mid12_info_part1_normr = mean( infodata(i).mid12.part1_normr_test_extrap(index2d) );

%    part1_data = [sta_info_part1 mid1_info_part1 mid2_info_part1 mid12_info_part1];

   [part1_data, err1] = error_check_info_values(sta_info_part1, mid1_info_part1, mid2_info_part1, mid12_info_part1, ...
                           sta_info_part1_normr, mid1_info_part1_normr, mid2_info_part1_normr, mid12_info_part1_normr);


   % Part 2 data:

   sta_info_part2 = infodata(i).sta.part2_test_extrap(index1d);
   mid1_info_part2 = infodata(i).mid1.part2_test_extrap(index1d);
   mid2_info_part2 = infodata(i).mid2.part2_test_extrap(index1d);
   mid12_info_part2 = mean( infodata(i).mid12.part2_test_extrap(index2d) );

   sta_info_part2_normr = infodata(i).sta.part2_normr_test_extrap(index1d);
   mid1_info_part2_normr = infodata(i).mid1.part2_normr_test_extrap(index1d);
   mid2_info_part2_normr = infodata(i).mid2.part2_normr_test_extrap(index1d);
   mid12_info_part2_normr = mean( infodata(i).mid12.part2_normr_test_extrap(index2d) );


%    part2_data = [sta_info_part2 mid1_info_part2 mid2_info_part2 mid12_info_part2];

   [part2_data, err2] = error_check_info_values(sta_info_part2, mid1_info_part2, mid2_info_part2, mid12_info_part2, ...
                              sta_info_part2_normr, mid1_info_part2_normr, mid2_info_part2_normr, mid12_info_part2_normr);


   % Part 3 data:

   sta_info_part3 = infodata(i).sta.part3_test_extrap(index1d);
   mid1_info_part3 = infodata(i).mid1.part3_test_extrap(index1d);
   mid2_info_part3 = infodata(i).mid2.part3_test_extrap(index1d);
   mid12_info_part3 = mean( infodata(i).mid12.part3_test_extrap(index2d) );

   sta_info_part3_normr = infodata(i).sta.part3_normr_test_extrap(index1d);
   mid1_info_part3_normr = infodata(i).mid1.part3_normr_test_extrap(index1d);
   mid2_info_part3_normr = infodata(i).mid2.part3_normr_test_extrap(index1d);
   mid12_info_part3_normr = mean( infodata(i).mid12.part3_normr_test_extrap(index2d) );

%    part3_data = [sta_info_part3 mid1_info_part3 mid2_info_part3 mid12_info_part3];

   [part3_data, err3] = error_check_info_values(sta_info_part3, mid1_info_part3, mid2_info_part3, mid12_info_part3, ...
                              sta_info_part3_normr, mid1_info_part3_normr, mid2_info_part3_normr, mid12_info_part3_normr);


   % Part 4 data:

   sta_info_part4 = infodata(i).sta.part4_test_extrap(index1d);
   mid1_info_part4 = infodata(i).mid1.part4_test_extrap(index1d);
   mid2_info_part4 = infodata(i).mid2.part4_test_extrap(index1d);
   mid12_info_part4 = mean( infodata(i).mid12.part4_test_extrap(index2d) );

   sta_info_part4_normr = infodata(i).sta.part4_normr_test_extrap(index1d);
   mid1_info_part4_normr = infodata(i).mid1.part4_normr_test_extrap(index1d);
   mid2_info_part4_normr = infodata(i).mid2.part4_normr_test_extrap(index1d);
   mid12_info_part4_normr = mean( infodata(i).mid12.part4_normr_test_extrap(index2d) );

%    part4_data = [sta_info_part4 mid1_info_part4 mid2_info_part4 mid12_info_part4];

   [part4_data, err4] = error_check_info_values(sta_info_part4, mid1_info_part4, mid2_info_part4, mid12_info_part4, ...
                              sta_info_part4_normr, mid1_info_part4_normr, mid2_info_part4_normr, mid12_info_part4_normr);


   info_values = [part1_data; part2_data; part3_data; part4_data];


   if ( ~isempty(info_values) )
      info_mean = mean(info_values,1);
      [info_mean, err] = error_check_info_values(info_mean(1), info_mean(2), info_mean(3), info_mean(4));
      info_mean_pop_test = [info_mean_pop_test; info_mean];
   end

   clear('part1_data', 'part2_data', 'part3_data', 'part4_data', 'info_values', 'info_mean');
   clear('sta_info_part1', 'mid1_info_part1', 'mid2_info_part1', 'mid12_info_part1');
   clear('sta_info_part2', 'mid1_info_part2', 'mid2_info_part2', 'mid12_info_part2');
   clear('sta_info_part3', 'mid1_info_part3', 'mid2_info_part3', 'mid12_info_part3');
   clear('sta_info_part4', 'mid1_info_part4', 'mid2_info_part4', 'mid12_info_part4');

   clear('sta_info_part1_normr', 'mid1_info_part1_normr', 'mid2_info_part1_normr', 'mid12_info_part1_normr');
   clear('sta_info_part2_normr', 'mid1_info_part2_normr', 'mid2_info_part2_normr', 'mid12_info_part2_normr');
   clear('sta_info_part3_normr', 'mid1_info_part3_normr', 'mid2_info_part3_normr', 'mid12_info_part3_normr');
   clear('sta_info_part4_normr', 'mid1_info_part4_normr', 'mid2_info_part4_normr', 'mid12_info_part4_normr');

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

% [length(infosta) length(infosta_test)]
% pause

% Plot the data

plot_extrap_info_scatter_test_data(infosta_test, infomid1_test, infomid2_test, infomid12_test);

% plot_extrapolated_information(infosta, infomid1, infomid2, infomid12);


% plot_observed_information(infosta_real, infomid1_real, infomid2_real, infomid12_real)
% 
% 
% plot_information_comparison(infosta_real, infomid1_real, infomid2_real, infomid12_real, ...
%                             infosta,      infomid1,      infomid2,      infomid12);

% plot_extrap_info_scatter(infosta, infomid1, infomid2, infomid12)

% plot_extrap_info_hist(infosta, infomid1, infomid2, infomid12)

%plot_extrap_info_synergy_scatter(infosta, infomid1, infomid2, infomid12)

% plot_synergy_versus_mids(infosta_test, infomid1_test, infomid2_test, infomid12_test);


% plot_extrap_info_hist_test_data(infosta_test, infomid1_test, infomid2_test, infomid12_test);
 



return; % end of Main




% ==================================================================================
% =========================== Function Declarations ================================
% ==================================================================================

function plot_extrapolated_information(infosta, infomid1, infomid2, infomid12)

% Plot the Extrapolated Information Values
% =======================================================
% =======================================================

figure;

% Plot the first mid against the sta
% -----------------------------------------
% [infomid1' infomid1_real' max([infomid1(:)'; infomid1_real(:)'])']
subplot(2,2,1);
plot(infosta, infomid1, 'ko');
mxmx = 1.05 * max([ max(infosta) max(infomid1) ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{STA}');
ylabel('I_{MID1}');
title('MID_{1} versus STA');



% Plot the first mid against the second mid
% -----------------------------------------
subplot(2,2,2);
plot(infomid1, infomid2, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
title('MID_{2} versus MID_{1}');



% Plot the first mid against the first and second mid
% ---------------------------------------------------
subplot(2,2,3);
plot(infomid1, infomid12, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{MID1}');
ylabel('I_{MID12}');
title('MID_{12} versus MID_{1}');




% Plot the synergy
% -----------------------------------------
subplot(2,2,4);
plot(infomid1+infomid2, infomid12, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{MID1} + I_{MID2}');
ylabel('I_{MID12}');
title('MID_{1} + MID_{2} versus MID_{12}');

suptitle('Extrapolation to Infinite Data Set Size');
print_mfilename(mfilename);
orient('landscape');

return; % end of function plot_extrapolated_information



function plot_observed_information(infosta_real, infomid1_real, infomid2_real, infomid12_real)

% Plot the Actual Information Values
% =======================================================
% =======================================================

figure;

% Plot the first mid against the sta
% -----------------------------------------
subplot(2,2,1);
plot(infosta_real, infomid1_real, 'ko');
mxmx = 1.05 * max([ max(infosta_real) max(infomid1_real) ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{STA}');
ylabel('I_{MID1}');
title('MID_{1} versus STA');



% Plot the first mid against the second mid
% -----------------------------------------
subplot(2,2,2);
plot(infomid1_real, infomid2_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
title('MID_{2} versus MID_{1}');



% Plot the first mid against the first and second mid
% ---------------------------------------------------
subplot(2,2,3);
plot(infomid1_real, infomid12_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{MID1}');
ylabel('I_{MID12}');
title('MID_{12} versus MID_{1}');


% Plot the synergy
% -----------------------------------------
subplot(2,2,4);
plot(infomid1_real+infomid2_real, infomid12_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('I_{MID1} + I_{MID2}');
ylabel('I_{MID12}');
title('MID_{1} + MID_{2} versus MID_{12}');

suptitle('Actual Data Set Size');
print_mfilename(mfilename);
orient('landscape');

return; % end of function plot_observed_information




function plot_information_comparison(infosta_real, infomid1_real, infomid2_real, infomid12_real, ...
                                     infosta,      infomid1,      infomid2,      infomid12)

% Compare the Information Values
% =======================================================
% =======================================================

figure;

% Plot the sta data
% -----------------------------------------
subplot(2,2,1);
plot(infosta, infosta_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('STA extrap');
ylabel('STA real');
[r,p]=corrcoef(infosta,infosta_real);
title(sprintf('STA: r=%.3f, p=%.6f', r(1,2), p(1,2)));


% Plot the mid1 data
% -----------------------------------------
subplot(2,2,2);
plot(infomid1, infomid1_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('MID1 extrap');
ylabel('MID1 real');
[r,p]=corrcoef(infomid1,infomid1_real);
title(sprintf('MID1: r=%.3f, p=%.6f', r(1,2), p(1,2)));


% Plot the mid2 data
% -----------------------------------------
subplot(2,2,3);
plot(infomid2, infomid2_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('MID2 extrap');
ylabel('MID2 real');
[r,p]=corrcoef(infomid2,infomid2_real);
title(sprintf('MID2: r=%.3f, p=%.6f', r(1,2), p(1,2)));


% Plot the mid12 data
% -----------------------------------------
subplot(2,2,4);
plot(infomid12, infomid12_real, 'ko');
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out');
axis('square');
xlabel('MID12 extrap');
ylabel('MID12 real');
[r,p]=corrcoef(infomid12,infomid12_real);
title(sprintf('MID12: r=%.3f, p=%.6f', r(1,2), p(1,2)));


suptitle('Comparison b/w Real and Extrapolated Info Values');
print_mfilename(mfilename);
orient('landscape');


return; % end of function plot_information_comparison





function plot_extrap_info_scatter(infosta, infomid1, infomid2, infomid12)


% Plot scatter and histograms of the information comparisons
% ==========================================================
% ==========================================================


figure;
markersize = 2;
% Plot the first mid against the sta
% -----------------------------------------
subplot(2,2,1);
plot(infosta, infomid1, 'ko', 'markersize', markersize);
mxmx = 3.5; %1.05 * max([ max(infosta) max(infomid1) ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.01 0.1 1 3.5];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{STA}');
ylabel('I_{MID1}');
title('MID_{1} versus STA');


% Plot the first mid against the second mid
% -----------------------------------------
subplot(2,2,2);
plot(infomid1, infomid2, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.01 mxmx], [0.01 mxmx], 'k-');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
xlim([0.01 mxmx])
ylim([0.01 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.01 0.1 1 4];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
title('MID_{2} versus MID_{1}');


% Plot the first mid against the first and second mid
% ---------------------------------------------------
subplot(2,2,3);
plot(infomid1, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.02 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID1}');
ylabel('I_{MID12}');
title('MID_{12} versus MID_{1}');


% Plot the synergy
% -----------------------------------------
subplot(2,2,4);
plot(infomid1+infomid2, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.02 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID1} + I_{MID2}');
ylabel('I_{MID12}');
title('MID_{1} + MID_{2} versus MID_{12}');

suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

return; % end of function plot_extrap_info_scatter_hist







function plot_extrap_info_scatter_test_data(infosta, infomid1, infomid2, infomid12)

% Plot scatter and histograms of the information comparisons
% ==========================================================
% ==========================================================

[infosta(:) infomid1(:)];

figure;
markersize = 2;
% Plot the first mid against the sta
% -----------------------------------------
subplot(2,2,1);
plot(infosta, infomid1, 'ko', 'markersize', markersize);
mxmx = 3.5; %1.05 * max([ max(infosta) max(infomid1) ]);
hold on;
plot([0.001 mxmx], [0.001 mxmx], 'k-');
xlim([0.001 mxmx])
ylim([0.001 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.001 0.01 0.1 1 3.5];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
box('off');
xlabel('STA Information');
ylabel('MID1 Information');


% Plot the first mid against the second mid
% -----------------------------------------
subplot(2,2,2);
plot(infomid1, infomid2, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.001 3], [0.001 3], 'k-');
xlabel('I_{MID1}');
ylabel('I_{MID2}');
xlim([0.001 3])
ylim([0.001 3])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.001 0.01 0.1 1 3];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
box('off');
xlabel('MID1 Information');
ylabel('MID2 Information');


% Plot the first mid against the first and second mid
% ---------------------------------------------------
subplot(2,2,3);
plot(infomid1, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.002 10], [0.002 10], 'k-');
xlim([0.002 10])
ylim([0.002 10])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.001 0.01 0.1 1 10];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
box('off');
xlabel('MID1 Information');
ylabel('MID12 Information');


% Plot the synergy
% -----------------------------------------
subplot(2,2,4);
plot(infomid1+infomid2, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.002 10], [0.002 10], 'k-');
xlim([0.002 10])
ylim([0.002 10])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.001 0.01 0.1 1 10];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
box('off');
xlabel('MID1 + MID2 Information');
ylabel('MID12 Information');

suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

return; % end of function plot_extrap_info_scatter_test_data







function plot_extrap_info_hist(infosta, infomid1, infomid2, infomid12)

% Plot scatter and histograms of the information comparisons
% ==========================================================
% ==========================================================

info1_info12 = 100 * infomid1 ./ infomid12;
info2_info12 = 100 * infomid2 ./ infomid12;
info1_info1_info2 = 100 * infomid1 ./ (infomid1 + infomid2);
info12_info1_info2 = 100 * infomid12 ./ (infomid1 + infomid2);
info2_info1 = 100 * infomid2 ./ infomid1;
infosta_info1 = 100 * infosta ./ infomid1;
infosta_info12 = 100 * infosta ./ infomid12;


[infosta(:) infomid1(:)];

[min(infosta_info1) max(infosta_info1)];


figure;
markersize = 2;

% Plot the first mid against the sta
% -----------------------------------------
subplot(2,2,1);
bins = linspace(0,100,11);
count = histc(infosta_info1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:5:25], 'yticklabel', [0:5:25]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * STA / MID_1 Information');
ylabel('Percent of neurons (%)');
title(sprintf('mn=%.2f, sd=%.2f', mean(infosta_info1), std(infosta_info1) ));


% Plot the first mid against the second mid
% -----------------------------------------
subplot(2,2,2);
bins = linspace(0,100,11);
count = histc(info2_info1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:5:35], 'yticklabel', [0:5:35]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * MID_2 / MID_1 Information');
ylabel('Percent of neurons (%)');
title(sprintf('mn=%.2f, sd=%.2f', mean(info2_info1), std(info2_info1) ));


% Plot the first mid against the first and second mid
% ---------------------------------------------------
subplot(2,2,3);
bins = linspace(0,100,11);
count = histc(info1_info12, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
set(gca,'ytick', [0:5:35], 'yticklabel', [0:5:35]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * MID_1 / MID_{12} Information');
ylabel('Percent of neurons (%)');
title(sprintf('mn=%.2f, sd=%.2f', mean(info1_info12), std(info1_info12) ));


% Plot the synergy
% -----------------------------------------
subplot(2,2,4);

% bins = linspace(75,300,10);

bins = linspace(75,200,11);

info12_info1_info2( info12_info1_info2 >= 200 ) = 199;

count = histc(info12_info1_info2, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([65 210 0 1.05*max(count)]);
set(gca,'xtick', [75:25:200], 'xticklabel', [75:25:200]);
set(gca,'ytick', [0:5:50], 'yticklabel', [0:5:50]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * MID_{12} / ( MID_1 + MID_2 ) Information');
ylabel('Percent of neurons(%)');
title(sprintf('mn=%.2f, sd=%.2f', mean(info12_info1_info2), std(info12_info1_info2) ));

[h,p,ci,stat] = ttest(info12_info1_info2, 100, 0.05, 1);

suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

length(find(info12_info1_info2>=125))/length(info12_info1_info2);

return; % end of function plot_extrap_info_scatter_hist




function plot_extrap_info_hist_test_data(infosta, infomid1, infomid2, infomid12)

% Plot scatter and histograms of the information comparisons
% ==========================================================
% ==========================================================

info1_info12 = 100 * infomid1 ./ infomid12;
info2_info12 = 100 * infomid2 ./ infomid12;
info1_info1_info2 = 100 * infomid1 ./ (infomid1 + infomid2);
info12_info1_info2 = 100 * infomid12 ./ (infomid1 + infomid2);
info2_info1 = 100 * infomid2 ./ infomid1;
infosta_info1 = 100 * infosta ./ infomid1;
infosta_info12 = 100 * infosta ./ infomid12;


[infomid1(:) infomid12(:)];


figure;
markersize = 2;

% Plot the first mid against the sta
% -----------------------------------------
subplot(2,2,1);
bins = linspace(0,100,11);
count = histc(infosta_info1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
box('off');
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
% set(gca,'ytick', [0:5:25], 'yticklabel', [0:5:25]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * STA / MID_1 Information');
ylabel('Percent of neurons (%)');
title(sprintf('md=%.2f', median(infosta_info1) ));


% Plot the first mid against the second mid
% -----------------------------------------
subplot(2,2,2);
bins = linspace(0,100,11);
count = histc(info2_info1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
box('off');
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
% set(gca,'ytick', [0:5:35], 'yticklabel', [0:5:35]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * MID_2 / MID_1 Information');
ylabel('Percent of neurons (%)');
title(sprintf('md=%.2f', median(info2_info1) ));


% Plot the first mid against the first and second mid
% ---------------------------------------------------
subplot(2,2,3);
bins = linspace(0,100,11);
count = histc(info1_info12, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
box('off');
% set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
% set(gca,'ytick', [0:5:35], 'yticklabel', [0:5:35]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('100 * MID_1 / MID_{12} Information');
ylabel('Percent of neurons (%)');
title(sprintf('md=%.2f', median(info1_info12) ));


% Plot the synergy
% -----------------------------------------
subplot(2,2,4);


info12_info1_info2(info12_info1_info2>=650) = 400;
info12_info1_info2(info12_info1_info2 < 81) = 81;

% Make sure the histogram catches all the data:

nbins = 14;
bins = 10.^( linspace(log10(80), log10(750), nbins) );
[count] = histc(info12_info1_info2, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
xlim([50 665]);
ylim([0 1.05*max(count)]);
set(hb,'facecolor', 0.75*ones(1,3) ); 
box off;
set(gca,'xtick', [100:100:700], 'xticklabel', [100:100:700]);
set(gca,'ytick', [0:5:50], 'yticklabel', [0:5:50]);
children = get(gca,'children');
set(children(1), 'marker', 'none');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('Synergy');
ylabel('Percent of neurons(%)');
title(sprintf('md=%.2f', median(info12_info1_info2) ));


% bins = linspace(75,500,18)
% count = histc(info12_info1_info2, bins);
% total_count = sum(count);
% count = 100 * count ./ total_count;
% hb = bar(bins, count, 'histc');
% set(hb, 'facecolor', [0.75 0.75 0.75])
% axis([60 510 0 1.05*max(count)]);
% % set(gca,'xtick', [50:25:300], 'xticklabel', [50:25:300]);
% % set(gca,'ytick', [0:5:50], 'yticklabel', [0:5:50]);
% set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% xlabel('100 * MID_{12} / ( MID_1 + MID_2 ) Information');
% ylabel('Percent of neurons(%)');
% title(sprintf('mn=%.2f, sd=%.2f', mean(info12_info1_info2), std(info12_info1_info2) ));

[h,p100,ci,stat] = ttest(info12_info1_info2, 100, 0.05, 1);
[h,p125,ci,stat] = ttest(info12_info1_info2, 125, 0.05, 1);
[h,p150,ci,stat] = ttest(info12_info1_info2, 150, 0.05, 1);

p100 = signrank(info12_info1_info2, 100, 0.05);
p125 = signrank(info12_info1_info2, 100, 0.05);
p150 = signrank(info12_info1_info2, 100, 0.05);

% [100 p100]
% [125 p125]
% [150 p150]


suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

% length(find(info12_info1_info2>=125))/length(info12_info1_info2)
% length(find(info12_info1_info2>=150))/length(info12_info1_info2)

return; % end of function plot_extrap_info_scatter_hist




function plot_extrap_info_synergy_scatter(infosta, infomid1, infomid2, infomid12)


% Plot scatter and histograms of the information comparisons
% ==========================================================
% ==========================================================

info1_info12 = 100 * infomid1 ./ infomid12;
info2_info12 = 100 * infomid2 ./ infomid12;
info1_info1_info2 = 100 * infomid1 ./ (infomid1 + infomid2);
info12_info1_info2 = 100 * infomid12 ./ (infomid1 + infomid2);
info2_info1 = 100 * infomid2 ./ infomid1;
infosta_info1 = 100 * infosta ./ infomid1;
infosta_info12 = 100 * infosta ./ infomid12;


[infosta(:) infomid1(:)];

[min(infosta_info1) max(infosta_info1)];



figure;
markersize = 2;


% Plot the synergy
% -----------------------------------------
subplot(1,2,1);
plot(infomid1+infomid2, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0.02 mxmx])
ylim([0.02 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.02 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID1} + I_{MID2}');
ylabel('I_{MID12}');
title('MID_{1} + MID_{2} versus MID_{12}');


subplot(1,2,2);
plot(infomid1+infomid2, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.02 mxmx], [0.02 mxmx], 'k-');
xlim([0 4])
ylim([0 4])
% set(gca,'xscale', 'log');
% set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0 2 4];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID1} + I_{MID2}');
ylabel('I_{MID12}');
title('MID_{1} + MID_{2} versus MID_{12}');

suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

return; % end of function plot_extrap_info_synergy_scatter




function plot_synergy_versus_mids(infosta, infomid1, infomid2, infomid12)


% Plot scatter and histograms of the information comparisons
% ==========================================================
% ==========================================================

info1_info12 = 100 * infomid1 ./ infomid12;
info2_info12 = 100 * infomid2 ./ infomid12;

info2_info1 = 100 * infomid2 ./ infomid1;

info1_info1_info2 = 100 * infomid1 ./ (infomid1 + infomid2);
info12_info1_info2 = 100 * infomid12 ./ (infomid1 + infomid2);

figure;
markersize = 3;

% Plot the synergy
% -----------------------------------------
subplot(1,2,1);

x = info1_info12;
y = info12_info1_info2;


[beta, s] = polyfit(log10(x), log10(y), 1);

[log10(10) log10(100) 10.^log10(min(x)) 10.^log10(max(x))]

xfit = logspace( log10(10), log10(100), 20 )

yfit = polyval(beta, log10(xfit) );
yfit = 10.^yfit;


plot( x , y , 'ko', 'markersize', markersize, 'markerfacecolor', 'k');
hold on;
plot( xfit, yfit, 'k-');
% plot([min(x) max(x)], [min(y) max(y)], 'k-');
xlim([10 110])
ylim([50 700])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xtick = [10 25 50 75 100];
ytick = [50 100 200 400];
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
% axis('square');
box('off');
xlabel('100 * MID1 / MID12 Information');
ylabel('Synergy');

[r, p] = corrcoef( x, y );
[r(1,2) p(1,2)];

[r, p] = corrcoef( log10(x), log10(y) );
[r(1,2) p(1,2)]

title(sprintf('r=%.2f, p=%.4f', r(1,2), p(1,2)));



% Plot the synergy
% -----------------------------------------
subplot(1,2,2);
x = info2_info1; % second mid compared to first mid
y = info12_info1_info2; % synergy
y=y(x>1);
x=x(x>1);

% y(y>200) = 200;
plot( x , y , 'ko', 'markersize', markersize, 'markerfacecolor', 'k');
hold on;
% plot([min(x) max(x)], [min(y) max(y)], 'k-');
xlim([0 100])
ylim([min(y) 650])
% set(gca,'xscale', 'log');
% set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% tick = [0.02 0.1 1 6];
set(gca,'xtick', 0:20:100, 'xticklabel', 0:20:100);
% set(gca,'ytick', tick, 'yticklabel', tick);
% axis('square');
xlabel('100* MID_{2} / MID_{1} Information');
ylabel('Synergy');
% title('MID_{1} + MID_{2} versus MID_{12}');

[r, p] = corrcoef( x, y );
[r(1,2) p(1,2)]

title(sprintf('r=%.2f, p=%.4f', r(1,2), p(1,2)));

suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

return; % end of function plot_synergy_versus_mids




function plot_training_versus_test_data(infosta, infomid1, infomid2, infomid12, infosta_test, infomid1_test, infomid2_test, infomid12_test)

figure;
markersize = 2;


% Plot the STA data
% -----------------------------------------
subplot(2,2,1);
plot(infosta_test, infosta, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.001 mxmx], [0.001 mxmx], 'k-');
xlim([0.001 mxmx])
ylim([0.001 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.01 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{STA-Test}');
ylabel('I_{STA-Train}');
title('STA Training vs. Test');


% Plot the MID1 data
% -----------------------------------------
subplot(2,2,2);
plot(infomid1_test, infomid1, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.001 mxmx], [0.001 mxmx], 'k-');
xlim([0.001 mxmx])
ylim([0.001 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.01 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID1-Test}');
ylabel('I_{MID1-Train}');
title('MID1 Training vs. Test');


% Plot the MID2 data
% -----------------------------------------
subplot(2,2,3);
plot(infomid2_test, infomid2, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.001 mxmx], [0.001 mxmx], 'k-');
xlim([0.001 mxmx])
ylim([0.001 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.01 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID2-Test}');
ylabel('I_{MID2-Train}');
title('MID2 Training vs. Test');


% Plot the MID12 data
% -----------------------------------------
subplot(2,2,4);
plot(infomid12_test, infomid12, 'ko', 'markersize', markersize);
mxmx = max([ get(gca,'xlim') get(gca,'ylim') ]);
hold on;
plot([0.001 mxmx], [0.001 mxmx], 'k-');
xlim([0.001 mxmx])
ylim([0.001 mxmx])
set(gca,'xscale', 'log');
set(gca,'yscale', 'log');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
tick = [0.01 0.1 1 6];
set(gca,'xtick', tick, 'xticklabel', tick);
set(gca,'ytick', tick, 'yticklabel', tick);
axis('square');
xlabel('I_{MID12-Test}');
ylabel('I_{MID12-Train}');
title('MID12 Training vs. Test');

suptitle('Extrapolated Data Set');
print_mfilename(mfilename);
orient('tall');

return;



function [data, err] = error_check_info_values(sta_info, mid1_info, mid2_info, mid12_info, sta_info_normr, mid1_info_normr, mid2_info_normr, mid12_info_normr)

data = [sta_info mid1_info mid2_info mid12_info];

% normr = [sta_info_normr mid1_info_normr mid2_info_normr mid12_info_normr];

if ( sta_info>mid12_info | sta_info>mid1_info | mid1_info>mid12_info | mid2_info>mid1_info | min(data)<0 | mid12_info<0.01 )
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





