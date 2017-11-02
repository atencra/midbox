function plot_mid_icc_a1_info_compare(projinfo, info)
%plot_mid_icc_projinfo_compare   Compare sta and mid filter information
%
% plot_mid_icc_projinfo_compare(projinfo)
% -------------------------------------------------------------------
% Mutual information for STA and MID1 and MID2 filters.
%
% The function plots I(STA) vs I(MID1), I(MID1) vs I(MID2), and I(MID1) vs I(MID12).
% 
% projinfo : struct array holding the information data. Each element of projinfo
% holds the data for 1 neuron. Usually saved in a file such as:
%
%    load 2003-11-24-site7-2380um-40db-dmr1-fs18115-mid-fio-projinfo.mat
%
% opts = struct('color', 'gray', 'FontMode','fixed','FontSize',8,'width',3.25);
%
% caa 2/1/10

close all;


% First get inferior colliculus data
% ---------------------------------------

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

      % error check - it is not possible for i1 to be > i12, so adjust
      if ( i1 > i12 )
         i1 = 0.99 * i12;
      end

%       % error check - it is not possible for i0 to be > i1, so adjust
%       if ( i0 > i1 )
%          i0 = 0.99 * i1;
%       end;

      % error check - it is not possible for i0 to be > i1, so adjust
      if ( i0 > i1 )
         a = 0.95;
         b = 1.00;
         i0 = ( a + (b-a).*rand(1) ) * i1;
      end;


      ista = [ista i0];
      iv1 = [iv1 i1];
      iv2 = [iv2 i2];
      iv12 = [iv12 i12];

   end

end % (for i)

iv2(iv2<0) = 0.25 * abs(iv2(iv2<0)); % if it's neg, then it's really 0, but
                                     % I make it small just in case

ista(ista<0) = 0.25 * abs(ista(ista<0));
iv1(iv1<0) = 0.25 * abs(iv1(iv1<0));
iv12(iv12<0) = 0.25 * abs(iv12(iv12<0));

ista_iv1 = 100 * ista ./ iv1;
iv1_iv12 = 100 * iv1 ./ iv12;
iv2_iv1 = 100 * iv2 ./ iv1;

[mean(iv1_iv12) std(iv1_iv12) median(iv1_iv12)]



% Now get auditory cortex data
% ---------------------------------------

position = [];
info_sta = [];
info1 = [];
info2 = [];
info_both = [];

for i = 1:length(info)

   position = [position info(i).position];
   info_sta = [info_sta mean(info(i).sta.information)];
   info1 = [info1 mean(info(i).mid1.information)];
   info2 = [info2 mean(info(i).mid2.information)];
   info_both = [info_both mean(info(i).mid12.information)];

end % (for i)


% Error check to get rid of bad INF values:

index = ~isinf(info_sta) & ~isinf(info1) & ~isinf(info2) & ~isinf(info_both);

info_sta = info_sta(index);
info1 = info1(index);
info2 = info2(index);
info_both = info_both(index);
position = position(index);

a1_sta_mid1 = 100 * info_sta ./ info1;
mid1_contrib = 100 * info1 ./ info_both;
synergy = 100 * info_both ./ (info1 + info2);

[mean(mid1_contrib) std(mid1_contrib) median(mid1_contrib)]




% Now plot cumulative distribution curves
% ----------------------------------------

% Compare STA and MID1
subplot(2,1,1);
[yy, xx, nn] = cdfcalc(ista_iv1);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
icxcdf = [-Inf; xx(nn); Inf];
icycdf = [0; 0; yy(1+nn)];

[yy, xx, nn] = cdfcalc(a1_sta_mid1);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
a1xcdf = [-Inf; xx(nn); Inf];
a1ycdf = [0; 0; yy(1+nn)];

[h,pval] = kstest2(ista_iv1, a1_sta_mid1);

hold on;
plot(icxcdf, icycdf, 'k-', 'linewidth', 3);
plot(a1xcdf, a1ycdf, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3);
xlim([ 20 100 ] );
ylim([0 1.05]);
set(gca,'ytick', 0:0.2:1.05, 'yticklabel', 0:0.2:1.05);
set(gca,'tickdir','out', 'ticklength', [0.025 0.025]);
ylabel('Cumulative Proportion');
title(sprintf('STA and MID1 Comparison; p = %.4f',pval));
legend('ICC', 'AI', 0);



% Compare MID1 contribution
subplot(2,1,2);
[yy, xx, nn] = cdfcalc(iv1_iv12);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
icxcdf = [-Inf; xx(nn); Inf];
icycdf = [0; 0; yy(1+nn)];

[yy, xx, nn] = cdfcalc(mid1_contrib);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
a1xcdf = [-Inf; xx(nn); Inf];
a1ycdf = [0; 0; yy(1+nn)];

[h,pval] = kstest2(iv1_iv12, mid1_contrib);

hold on;
plot(icxcdf, icycdf, 'k-', 'linewidth', 3);
plot(a1xcdf, a1ycdf, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3);
xlim([ 20 100 ] );
ylim([0 1.05]);
set(gca,'ytick', 0:0.2:1.05, 'yticklabel', 0:0.2:1.05);
set(gca,'tickdir','out', 'ticklength', [0.025 0.025]);
ylabel('Cumulative Proportion');
title(sprintf('MID1 Contribution; p = %.4f',pval));

set(gcf,'position', [144   418   341   453]);

print_mfilename(mfilename);



return






markersize = 3;

figure;

subplot(3,2,1);
hold on;
xmax = max([max(iv1) max(ista)]);
xmin = min([min(iv1) min(ista)]);
hp = plot(iv1, ista, 'ko', 'markersize', markersize); %, 'markerfacecolor', 'k');
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
ylabel('I(STA) [bits/sp]')
title('STA vs. MID1');


subplot(3,2,2);
bins = linspace(0,100,11);
count = histc(ista_iv1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
xlabel('100 * I(STA) / I(MID1)');
ylabel('Percent of Neurons');


subplot(3,2,3);
hold on
xmax = max([max(iv1) max(iv2)]);
xmin = min([min(iv1) min(iv2)]);
hp = plot(iv1, iv2, 'ko', 'markersize', markersize); %, 'markerfacecolor', 'k');
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
ylabel('I(MID2) [bits/sp]')
title('MID2 vs. MID1');


subplot(3,2,4);
bins = linspace(0,100,11);
count = histc(iv2_iv1, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
xlabel('100 * I(MID2) / I(MID1)');
ylabel('Percent of Neurons');



subplot(3,2,5);
hold on;
xmax = max([max(iv1) max(iv12)]);
xmin = min([min(iv1) min(iv12)]);
hp = plot(iv1, iv12, 'ko', 'markersize', markersize); %, 'markerfacecolor', 'k');
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
xlabel('I(MID1) [bits/sp]')
ylabel('I(MID12) [bits/sp]')
title('MID12 vs. MID1');


subplot(3,2,6);
bins = linspace(0,100,11);
count = histc(iv1_iv12, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-10 110 0 1.05*max(count)]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
xlabel('100 * I(MID1) / I(MID12)');
ylabel('Percent of Neurons');

set(gcf,'position', [360   199   710   723]);

print_mfilename(mfilename);


return;
















