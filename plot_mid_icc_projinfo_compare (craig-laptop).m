function plot_mid_icc_projinfo_compare(projinfo)
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
% opts = struct('color', 'gray', 'FontMode','fixed','FontSize',8,'width',4.5);
%
% caa 2/1/10

close all;


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

% Just remove bad data points later - will skip the following code
% for the final figures
%
%       % error check - it is not possible for i1 to be > i12, so adjust
%       if ( i1 > i12 )
%          i1 = 0.99 * i12;
%       end
% 
%       % error check - it is not possible for i0 to be > i1, so adjust
%       if ( i0 > i1 )
%          i0 = 0.99 * i1;
%       end;

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
index = (ista<=iv1) & (iv1<=iv12);
ista = ista(index);
iv1 = iv1(index);
iv2 = iv2(index);
iv12 = iv12(index);


% iv2(iv2<0) = 0.25 * abs(iv2(iv2<0)); % if it's neg, then it's really 0, but
%                                      % I make it small just in case
% ista(ista<0) = 0.25 * abs(ista(ista<0));
% iv1(iv1<0) = 0.25 * abs(iv1(iv1<0));
% iv12(iv12<0) = 0.25 * abs(iv12(iv12<0));

ista_iv1 = 100 * ista ./ iv1;
iv1_iv12 = 100 * iv1 ./ iv12;
iv2_iv1 = 100 * iv2 ./ iv1;

[r,p] = corrcoef( log10(iv1), log10(iv12) )
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
set(gca,'ytick', [0:20:100], 'yticklabel', [0:20:100]);
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
set(gca,'ytick', [0:30:100], 'yticklabel', [0:30:100]);
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
set(gca,'ytick', [0:30:100], 'yticklabel', [0:30:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
xlabel('100 * I(MID1) / I(MID12)');
ylabel('Percent of Neurons');

set(gcf,'position', [360   199   710   723]);

print_mfilename(mfilename);


return;
