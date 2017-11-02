function plot_mid_icc_projinfo_fio_curve_fit_params(projinfo, params0, params1, nmse0, nmse1)
%plot_mid_icc_fr_pli_resptype   Response index, firing rate, phase locking, and filter info 
%
% plot_mid_icc_fr_pli_resptype(projinfo)
% -------------------------------------------------------------------
% params : strf params struct array. Obtained from strf_parameters.m
%
% rt : response type struct array. Obtained from response_type_from_fra_raster.m
%
% projinfo : filter information struct array. Obtained from get_fio_info_from_filtstr.m
%
% caa 2/1/10



if ( nargin ~= 5 )
   error('You need 2 or 4 input arguments.');
end



% Get information metrics
ista = zeros(1,length(projinfo));
iv1 = zeros(1,length(projinfo));
iv2 = zeros(1,length(projinfo));
iv12 = zeros(1,length(projinfo));

for i = 1:length(projinfo)

   ista(i) = projinfo(i).info0_extrap_test(end);
   iv1(i) = projinfo(i).info1_extrap_test(end);
   iv2(i) = projinfo(i).info2_extrap_test(end);
   iv12(i) = projinfo(i).info12_extrap_test(end);

end % (for i)

iv2(iv2<0) = 0.25 * abs(iv2(iv2<0)); % if it's neg, then it's really 0, but
                                     % I make it small just in case
ista_iv1 = 100 * ista ./ iv1;
iv1_iv12 = 100 * iv1 ./ iv12;
iv2_iv1 = 100 * iv2 ./ iv1;

i01 = find(nmse0 < 0.25 & nmse1 < 0.25);

params0 = params0(i01,:);
params1 = params1(i01,:);

gr0 = params0(:,2);
tr0 = params0(:,4);

gr1 = params1(:,2);
tr1 = params1(:,4);

ista_iv1 = ista_iv1(i01);
iv1_iv12 = iv1_iv12(i01);
iv2_iv1 = iv2_iv1(i01);

close all;


figure;


% Growth Rate
%------------------------
subplot(2,2,1);
hold on;
xmin = min(gr1);
xmax = max(gr1);
xrange = xmax - xmin;
ymin = min(ista_iv1);
ymax = max(ista_iv1);
yrange = ymax - ymin;
hp = plot(gr1, ista_iv1, 'ko', 'markersize', 4);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([xmin-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('Growth Rate [sp/s/SD]');
ylabel('I(STA) / I(V1)')
title('I(STA) / I(MID1) vs Growth Rate');


subplot(2,2,2);
hold on;
xmin = min(gr1);
xmax = max(gr1);
xrange = xmax - xmin;
ymin = min(iv1_iv12);
ymax = max(iv1_iv12);
yrange = ymax - ymin;
hp = plot(gr1, iv1_iv12, 'ko', 'markersize', 4);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([xmin-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('Growth Rate [sp/s/SD]');
ylabel('I(V1) / I(V12)');
title('MID1 Contrib vs Growth Rate');


% Transition
%------------------------

subplot(2,2,3);
hold on;
xmin = min(tr1);
xmax = max(tr1);
xrange = xmax - xmin;
ymin = min(ista_iv1);
ymax = max(ista_iv1);
yrange = ymax - ymin;
hp = plot(tr1, ista_iv1, 'ko', 'markersize', 4);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([xmin-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('Transition [SD]');
ylabel('I(STA) / I(V1)')
title('I(STA) / I(MID1) vs Transition');


subplot(2,2,4);
hold on;
xmin = min(tr1);
xmax = max(tr1);
xrange = xmax - xmin;
ymin = min(iv1_iv12);
ymax = max(iv1_iv12);
yrange = ymax - ymin;
hp = plot(tr1, iv1_iv12, 'ko', 'markersize', 4);
set(hp, 'markerfacecolor', 'k');
set(hp, 'markeredgecolor', 'k');
xlim([xmin-0.05*xrange xmax+0.05*xrange]);
ylim([0 105]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% set(gca,'xtick', [0:20:100], 'xticklabel', [0:20:100]);
set(gca,'ytick', [0:25:100], 'yticklabel', [0:25:100]);
xlabel('Transition [SD]');
ylabel('I(V1) / I(V12)');
title('MID1 Contrib vs Transition');


print_mfilename(mfilename);

return;

figure;

subplot(2,2,1);

edges = linspace(0, 10, nbins);
np = histc(params1(:,2), edges);
hb = bar(edges, np, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlim([edges(1) edges(end)]);
box off;
title('Growth Rate');

subplot(2,2,2);
edges = linspace(-4, 6, nbins);
np = histc(params1(:,4), edges);
hb = bar(edges, np, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlim([edges(1)-1 edges(end)+1]);
% ylim([]);
box off;
title('Transition (SD)');

subplot(2,2,3);
[r,p] = corrcoef(params1(:,4), params1(:,2));
plot(params1(:,4), params1(:,2), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
box off;
title(sprintf('MID1 Transition Vs. Growth Rate\nr=%.3f, p=%.3f', r(1,2), p(1,2)) );
ylabel('MID1 Growth Rate');
xlabel('MID1 Transition (SD)');

suptitle('MID1 Nonlinearity Curve Fit Params');

print_mfilename(mfilename);


return;






