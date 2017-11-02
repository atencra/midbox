function plot_mid_icc_fr_pli_resptype(params, rt, projinfo)
%plot_mid_icc_fr_pli_resptype   Response index, firing rate, phase locking, and filter info 
%
% plot_mid_icc_fr_pli_resptype(params, rt, projinfo)
% -------------------------------------------------------------------
% params : strf params struct array. Obtained from strf_parameters.m
%
% rt : response type struct array. Obtained from response_type_from_fra_raster.m
%
% projinfo : filter information struct array. Obtained from get_fio_info_from_filtstr.m
%
% opts = struct('color', 'gray', 'FontMode','fixed','FontSize',8,'width',4.5);
%
% caa 2/1/10



% % Get response metrics
% fr = [params.w0];
% pli = [params.pli];
% nb_nb_ne = [rt.nb_nb_ne_total];
% ne_nb_ne = [rt.ne_nb_ne_total];
% 
% % Get information metrics
% ista = zeros(1,length(projinfo));
% iv1 = zeros(1,length(projinfo));
% iv2 = zeros(1,length(projinfo));
% iv12 = zeros(1,length(projinfo));
% 
% for i = 1:length(projinfo)
% i
%    ista(i) = projinfo(i).info0_extrap_test(end);
%    iv1(i) = projinfo(i).info1_extrap_test(end);
%    iv2(i) = projinfo(i).info2_extrap_test(end);
%    iv12(i) = projinfo(i).info12_extrap_test(end);
% 
% end % (for i)



% Get response metrics
fr = [];
pli = [];
nb_nb_ne = [];
ne_nb_ne = [];

ista = []; %zeros(1,length(projinfo));
iv1 = []; %zeros(1,length(projinfo));
iv2 = []; %zeros(1,length(projinfo));
iv12 = []; %zeros(1,length(projinfo));

for i = 1:length(projinfo)

   info0 = projinfo(i).info0_extrap_test;
   info1 = projinfo(i).info1_extrap_test;
   info2 = projinfo(i).info2_extrap_test;
   info12 = projinfo(i).info12_extrap_test;

   fr_temp = params(i).w0;
   pli_temp = params(i).pli;
   nb_nb_ne_temp = rt(i).nb_nb_ne_total;
   ne_nb_ne_temp = rt(i).ne_nb_ne_total;

   if ( ~isempty(info0) && ~isempty(info1) && ~isempty(info2) && ~isempty(info12) )

      i0 = info0(end);
      i1 = info1(end);
      i2 = info2(end);
      i12 = info12(end);


% Won't do this for final figures - will do a more conservative check - and will
% remove bad points instead of adjusting points
%
%       % error check - it is not possible for i1 to be > i12, so adjust
%       if ( i1 > i12 )
%          i1 = 0.99 * i12;
%       end
% 
%       % error check - it is not possible for i0 to be > i1, so adjust
%       if ( i0 > i1 )
%          a = 0.95;
%          b = 1.00;
%          i0 = ( a + (b-a).*rand(1) ) * i1;
%       end;

      ista = [ista i0];
      iv1 = [iv1 i1];
      iv2 = [iv2 i2];
      iv12 = [iv12 i12];

   else

      ista = [ista nan];
      iv1 = [iv1 nan];
      iv2 = [iv2 nan];
      iv12 = [iv12 nan];

   end

      fr = [fr fr_temp];
      pli = [pli pli_temp];
      nb_nb_ne = [nb_nb_ne nb_nb_ne_temp];
      ne_nb_ne = [ne_nb_ne ne_nb_ne_temp];

end % (for i)

% iv2(iv2<0) = 0.25 * abs(iv2(iv2<0)); % if it's neg, then it's really 0, but
%                                      % I make it small just in case
% 
% ista(ista<0) = 0.25 * abs(ista(ista<0));
% iv1(iv1<0) = 0.25 * abs(iv1(iv1<0));
% iv12(iv12<0) = 0.25 * abs(iv12(iv12<0));

% Original non-info data

fr_orig = fr;
pli_orig = pli;
nb_nb_ne_orig = nb_nb_ne;
ne_nb_ne_orig = ne_nb_ne;



% Can't have negative information values, so get rid of these points
index = (ista>0) & (iv1>0) & (iv2>0) & (iv12>0);
ista = ista(index);
iv1 = iv1(index);
iv2 = iv2(index);
iv12 = iv12(index);
fr = fr(index);
pli = pli(index);
ne_nb_ne = ne_nb_ne(index);
nb_nb_ne = nb_nb_ne(index);


% Can't have sta > iv1, or iv1 > iv12, so get rid of these points
index = (ista<=iv1) & (iv1<=iv12);
ista = ista(index);
iv1 = iv1(index);
iv2 = iv2(index);
iv12 = iv12(index);
fr = fr(index);
pli = pli(index);
ne_nb_ne = ne_nb_ne(index);
nb_nb_ne = nb_nb_ne(index);


ista_iv1 = 100 * ista ./ iv1;
iv1_iv12 = 100 * iv1 ./ iv12;
iv2_iv1 = 100 * iv2 ./ iv1;


markersize = 3;

close all;

figure;

subplot(3,2,1);
hold on;
xmax = max(fr_orig);
xmin = min(fr_orig);
xrange = xmax - xmin;
ymax = max(pli_orig);
ymin = min(pli_orig);
yrange = ymax - ymin;
hp = plot(fr_orig, pli_orig, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:0.1:1], 'yticklabel', [0:0.1:1]);
% xlim([xmin-0.05*xrange xmax+0.05*xrange]);
% ylim([ymin-0.05*yrange ymax+0.05*yrange]);
xlim([0 xmax+0.05*xrange]);
ylim([0 ymax+0.05*yrange]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('FR [sp/s]')
ylabel('PLI')

% [r,p] = corrcoef(fr_orig, pli_orig);
% fprintf('\nFR vs PLI: r = %.3f, p = %.4f\n', r(1,2), p(1,2));

subplot(3,2,3);
hold on;
xmax = max(fr_orig);
xmin = min(fr_orig);
xrange = xmax - xmin;
ymax = max(ne_nb_ne_orig);
ymin = min(ne_nb_ne_orig);
yrange = ymax - ymin;
plot([xmin xmax], [0.5 0.5], 'k--');
hp = plot(fr_orig, ne_nb_ne_orig, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:0.1:1], 'yticklabel', [0:0.1:1]);
% xlim([xmin-0.05*xrange xmax+0.05*xrange]);
xlim([0 xmax+0.05*xrange]);
% ylim([ymin-0.05*yrange ymax+0.05*yrange]);
ylim([0 0.6]); %1.05*ymax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('FR [sp/s]')
ylabel('Response Index')

% [r,p] = corrcoef(fr_orig, ne_nb_ne_orig);
% fprintf('FR vs Response Index: r = %.3f, p = %.4f\n', r(1,2), p(1,2));




subplot(3,2,5);
hold on;
xmax = max(pli);
xmin = min(pli);
xrange = xmax - xmin;
ymax = max(ne_nb_ne_orig);
ymin = min(ne_nb_ne_orig);
yrange = ymax - ymin;
% plot([xmin xmax], [0.5 0.5], 'k--');
plot([0 xmax], [0.5 0.5], 'k--');
hp = plot(pli_orig, ne_nb_ne_orig, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
set(gca,'xtick', [0:0.1:1], 'xticklabel', [0:0.1:1]);
set(gca,'ytick', [0:0.1:1], 'yticklabel', [0:0.1:1]);
% xlim([xmin-0.05*xrange xmax+0.05*xrange]);
xlim([0 xmax+0.05*xrange]);
ylim([0 0.6]); %1.05*ymax]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('PLI')
ylabel('Response Index')

% [r,p] = corrcoef(pli_orig, ne_nb_ne_orig);
% fprintf('PLI vs Response Index: r = %.3f, p = %.4f\n', r(1,2), p(1,2));



disp('fr_orig');
% median(fr_orig)
% mad(fr_orig,1)

subplot(3,2,2);
bins = linspace(0,100,11);
count = histc(fr_orig, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
box off;
% axis([-10 110 0 1.05*max(count)]);
xlim([min(bins)-0.05*max(bins) max(bins)+0.05*max(bins)]);
ylim([ 0 1.05*max(count)]);
set(gca,'xtick', [0:2*unique(diff(bins)):100], 'xticklabel', [0:2*unique(diff(bins)):100]);
set(gca,'ytick', [0:10:100], 'yticklabel', [0:10:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('FR [sp/s]');
ylabel('Percent of Neurons');

disp('pli_orig');
% median(pli_orig)
% mad(pli_orig,1)

subplot(3,2,4);
bins = linspace(0,0.500,11);
count = histc(pli_orig, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75]);
box off;
xlim([min(bins)-0.05*max(bins) max(bins)+0.05*max(bins)]);
ylim([ 0 1.05*max(count)]);
set(gca,'xtick', [0:2*unique(diff(bins)):1], 'xticklabel', [0:2*unique(diff(bins)):1]);
set(gca,'ytick', [0:10:100], 'yticklabel', [0:10:100]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('PLI');
ylabel('Percent of Neurons');

disp('ne_nb_ne_orig');
% median(ne_nb_ne_orig)
% mad(ne_nb_ne_orig,1)

subplot(3,2,6);
bins = linspace(0,0.60,13);
count = histc(ne_nb_ne_orig, bins);
total_count = sum(count);
count = 100 * count ./ total_count;
max(count)
ytick = linspace(0,25,5);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
box off;
axis([min(bins)-0.05*0.10 max(bins)+0.05*0.10 0 1.05*max(count)]);
set(gca,'xtick', [0:2*unique(diff(bins)):1], 'xticklabel', [0:2*unique(diff(bins)):1]);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
axis([min(bins)-0.05*max(bins) max(bins)+0.05*max(bins) 0 1.05*max(count)]);
xlabel('Response Index');
ylabel('Percent of Neurons');

set(gcf,'position', [360   199   710   723]);

print_mfilename(mfilename);

% close all;


figure;

subplot(3,2,1);
hold on;
xmin = min(fr);
xmax = max(fr);
xrange = xmax - xmin;
ymin = min(ista_iv1);
ymax = max(ista_iv1);
yrange = ymax - ymin;
hp = plot(fr, ista_iv1, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([xmin-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('FR [sp/s]');
ylabel('100 I(STA) / I(MID1)')
title('STA / MID1');

% [r,p] = corrcoef(fr, ista_iv1);
% fprintf('\nFR vs Info STA/MID1: r = %.3f, p = %.4f\n', r(1,2), p(1,2));


subplot(3,2,3);
hold on;
xmin = min(pli);
xmax = max(pli);
xrange = xmax - xmin;
ymin = min(ista_iv1);
ymax = max(ista_iv1);
yrange = ymax - ymin;
hp = plot(pli, ista_iv1, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([0-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', [0:0.1:1], 'xticklabel', [0:0.1:1]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('PLI');
ylabel('100 I(STA)/I(MID1)')

% [r,p] = corrcoef(pli, ista_iv1);
% fprintf('PLI vs Info STA/MID1: r = %.3f, p = %.4f\n', r(1,2), p(1,2));


subplot(3,2,5);
hold on;
xmin = min(ne_nb_ne);
xmax = max(ne_nb_ne);
xrange = xmax - xmin;
ymin = min(ista_iv1);
ymax = max(ista_iv1);
yrange = ymax - ymin;
hp = plot(ne_nb_ne, ista_iv1, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
plot([0.5 0.5],[0 105],'k--');
xlim([0 0.6]); %1.05*ymax]);
% xlim([0-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', [0:0.1:1], 'xticklabel', [0:0.1:1]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('Response Index');
ylabel('100 I(STA) / I(MID1)')

% [r,p] = corrcoef(ne_nb_ne, ista_iv1);
% fprintf('Response Index vs Info STA/MID1: r = %.3f, p = %.4f\n', r(1,2), p(1,2));


subplot(3,2,2);
hold on;
xmin = min(fr);
xmax = max(fr);
xrange = xmax - xmin;
ymin = min(iv1_iv12);
ymax = max(iv1_iv12);
yrange = ymax - ymin;
hp = plot(fr, iv1_iv12, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([xmin-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('FR [sp/s]');
ylabel('100 I(MID1) / I(MID12)');
title('MID1 / MID12');

% [r,p] = corrcoef(fr, iv1_iv12);
% fprintf('FR vs Info MID1/MID12: r = %.3f, p = %.4f\n', r(1,2), p(1,2));


subplot(3,2,4);
hold on;
xmin = min(pli);
xmax = max(pli);
xrange = xmax - xmin;
ymin = min(iv1_iv12);
ymax = max(iv1_iv12);
yrange = ymax - ymin;
hp = plot(pli, iv1_iv12, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([0-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', [0:0.1:1], 'xticklabel', [0:0.1:1]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('PLI');
ylabel('100 I(MID1) / I(MID12)')

% [r,p] = corrcoef(pli, iv1_iv12);
% fprintf('PLI vs Info MID1/MID12: r = %.3f, p = %.4f\n', r(1,2), p(1,2));


subplot(3,2,6);
hold on;
xmin = min(ne_nb_ne);
xmax = max(ne_nb_ne);
xrange = xmax - xmin;
ymin = min(iv1_iv12);
ymax = max(iv1_iv12);
yrange = ymax - ymin;
hp = plot(ne_nb_ne, iv1_iv12, 'ko', 'markersize', markersize);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
plot([0.5 0.5],[0 105],'k--');
xlim([0 0.6]); %1.05*ymax]);
% xlim([0-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca,'xtick', [0:0.1:1], 'xticklabel', [0:0.1:1]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('Response Index');
ylabel('100 I(MID1) / I(MID12)')

% [r,p] = corrcoef(ne_nb_ne, iv1_iv12)
% fprintf('Response Index vs Info MID1/MID12: r = %.3f, p = %.4f\n', r(1,2), p(1,2));


set(gcf,'position', [360   199   710   723]);

print_mfilename(mfilename);

return;





