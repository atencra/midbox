function [params0, params1, nmse0, nmse1] = plot_nonlinearity_curve_fit(fio, sr)
% plot_nonlinearity_curve_fit   Curve fit to 1D nonlinearities
%
% [params0, params1, nmse0, nmse1] = plot_nonlinearity_curve_fit(fio, sr)
% =======================================================================
% fio : struct array holding nonlinearity data. Each element represents
% the data for one neuron.
%
% sr : sampling rate of the spike train. If the spike train was binned
% in 1 ms bins, then sr = 1000. For 5 ms resolution, sr = 200. This is
% needed so that the nonlinearities can be scaled appropriately.
% 
% params0 and params1 are NxM arrays. N is the number of cells, which
% is the same as length(fio), and M is the number of parameters in
% the curve fit.
%
% params(1) : baseline
% params(2) : growth rate
% params(3) : width
% params(4) : transition
%
% nmse0 and nmse1 are the normalized mean squared error of the 
% curve fits. They are Nx1 vectors, with each row corresponding
% to the nmse for the corresponding neuron in fio.
%
% caa 2/1/10

close all;

params0 = [];
params1 = [];

nmse0 = [];
nmse1 = [];


params0_rm = [];
params1_rm = [];

nmse0_rm = [];
nmse1_rm = [];


for i = 1:length(fio)

   xbins = fio(i).x0bins;

   pspk = fio(i).pspk;
   pspkx0 = fio(i).pspkx0;
   pspkx1 = fio(i).pspkx1;

   pspkx0_mtx = [];
   pspkx1_mtx = [];

   for j = 1:length(pspkx0)
      pspkx0_mtx = [pspkx0_mtx pspkx0{j}];
      pspkx1_mtx = [pspkx1_mtx pspkx1{j}];
   end

   pspkx0_mn = nanmean(pspkx0_mtx,2);
   pspkx1_mn = nanmean(pspkx1_mtx,2);

   index0 = ~isnan(pspkx0_mn);
   index1 = ~isnan(pspkx1_mn);

   pspkx0 = sr * double( pspkx0_mn(index0) );
   pspkx1 = sr * double( pspkx1_mn(index1) );

   x0bins = xbins(index0);
   x1bins = xbins(index1);


   % Get rid of non-monotonic effects
   temp_x0bins = x0bins;
   temp_pspkx0 = pspkx0;
   while ( temp_pspkx0(end) < temp_pspkx0(end-1) )
      temp_pspkx0 = temp_pspkx0(1:end-1);
      temp_x0bins = temp_x0bins(1:end-1);
   end


   temp_x1bins = x1bins;
   temp_pspkx1 = pspkx1;
   while ( temp_pspkx1(end) < temp_pspkx1(end-1) )
      temp_pspkx1 = temp_pspkx1(1:end-1);
      temp_x1bins = temp_x1bins(1:end-1);
   end


   % STA curve fit
   % --------------------------------------
   a0 = [1 1 1 1];
   [a0, resnorm] = lsqcurvefit(@sachs_abbas_func, a0, temp_x0bins, temp_pspkx0);
   params0 = [params0; a0(:)'];
   xfit = linspace(min(x0bins), max(x0bins), 100);
   fxfit = sachs_abbas_func(a0, xfit);

   temp = sachs_abbas_func(a0, x0bins);
   [gf0] = gfit2(temp_pspkx0, temp(1:length(temp_pspkx0)), '2'); % gets goodness of fit measure, nmse
   nmse0 = [nmse0; gf0];

   clf;

   subplot(2,1,1);
   hold on;
   plot(xfit, fxfit, 'k-', 'linewidth', 2);
   plot(x0bins, pspkx0, 'ko', 'markerfacecolor', 'k');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   legend('SA Fit', 'Data', 'location', 'northwest');
   title(sprintf('#%.0f : STA : NMSE = %.3f', i, gf0));
   xlim([min(x0bins) max(x0bins)]);
   ylabel('Firing Rate (sp/s)');


   % MID1 curve fit
   % --------------------------------------
   a1 = [1 1 1 1];
   [a1, resnorm] = lsqcurvefit(@sachs_abbas_func, a1, temp_x1bins, temp_pspkx1);
   params1 = [params1; a1(:)'];
   xfit = linspace(min(x1bins), max(x1bins), 100);
   fxfit = sachs_abbas_func(a1, xfit);

   temp = sachs_abbas_func(a1, x1bins);
   [gf1] = gfit2(temp_pspkx1, temp(1:length(temp_pspkx1)), '2'); % gets goodness of fit measure, nmse
   nmse1 = [nmse1; gf1];

   subplot(2,1,2);
   hold on;
   plot(xfit, fxfit, 'k-', 'linewidth', 2);
   plot(x1bins, pspkx1, 'ko', 'markerfacecolor', 'k');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   title(sprintf('MID1 : NMSE = %.3f', gf1));
   xlim([min(x1bins) max(x1bins)]);
   xlabel('Projection (SD)');
   ylabel('Firing Rate (sp/s)');

   print_mfilename(mfilename);

   set(gcf, 'position', [360   340   488   582]);




   %************************* try a different curve fit **********************

   % STA curve fit
   % --------------------------------------
%    beta0 = [1 1 1];
% 
%    [beta] = nlinfit(temp_x0bins, temp_pspkx0, @ringach_malone_func, beta0);
% 
%    [a0, resnorm] = lsqcurvefit(@ringach_malone_func, a0, temp_x0bins, temp_pspkx0);
% 
%    params0_rm = [params0_rm; beta(:)'];
%    xfit = linspace(min(x0bins), max(x0bins), 100);
%    fxfit = ringach_malone_func(a0, xfit);
% 
%    temp = ringach_malone_func(a0, x0bins);
%    [gf0] = gfit2(temp_pspkx0, temp(1:length(temp_pspkx0)), '2'); % gets goodness of fit measure, nmse
%    nmse0_rm = [nmse0_rm; gf0];



%{

   beta0 = [1 1 1];

   [beta] = nlinfit(temp_x0bins, temp_pspkx0, @ringach_malone_func, beta0);

   params0_rm = [params0_rm; beta(:)'];
   xfit = linspace(min(x0bins), max(x0bins), 100);
   fxfit = ringach_malone_func(beta, xfit);

   temp = ringach_malone_func(beta, x0bins);
   [gf0] = gfit2(temp_pspkx0, temp(1:length(temp_pspkx0)), '2'); % gets goodness of fit measure, nmse
   nmse0_rm = [nmse0_rm; gf0];




   figure;

   subplot(2,1,1);
   hold on;
   plot(xfit, fxfit, 'c-', 'linewidth', 2);
   plot(x0bins, pspkx0, 'ko', 'markerfacecolor', 'k');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   legend('RM Fit', 'Data', 'location', 'northwest');
   title(sprintf('#%.0f : STA : NMSE = %.3f', i, gf0));
   xlim([min(x0bins) max(x0bins)]);
   ylabel('Firing Rate (sp/s)');


   % MID1 curve fit
   % --------------------------------------
   a1 = [1 1 1];
   [a1, resnorm] = lsqcurvefit(@ringach_malone_func, a1, temp_x1bins, temp_pspkx1);
   params1_rm = [params1_rm; a1(:)'];
   xfit = linspace(min(x1bins), max(x1bins), 100);
   fxfit = ringach_malone_func(a1, xfit);

   temp = ringach_malone_func(a1, x1bins);
   [gf1] = gfit2(temp_pspkx1, temp(1:length(temp_pspkx1)), '2'); % gets goodness of fit measure, nmse
   nmse1_rm = [nmse1_rm; gf1];

   subplot(2,1,2);
   hold on;
   plot(xfit, fxfit, 'c-', 'linewidth', 2);
   plot(x1bins, pspkx1, 'ko', 'markerfacecolor', 'k');
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   title(sprintf('MID1 : NMSE = %.3f', gf1));
   xlim([min(x1bins) max(x1bins)]);
   xlabel('Projection (SD)');
   ylabel('Firing Rate (sp/s)');

   print_mfilename(mfilename);

   set(gcf, 'position', [360   340   488   582]);
%}

   pause

end % (for)


return;


























index = find(nmse <= 0.15);
params_good = params(index,:);
index = find(params_good(:,4) > -3);
params_good = params_good(index,:);

figure;

subplot(2,3,1);
hist(params_good(:,2), 25);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Growth Rate');

subplot(2,3,2);
hist(params_good(:,3), 25);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Width');

subplot(2,3,3);
hist(params_good(:,4), 25);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition (SD)');

subplot(2,3,4);
plot(params_good(:,2), params_good(:,3), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Width Vs. Growth Rate');
xlabel('Growth Rate');
ylabel('Width');

subplot(2,3,5);
plot(params_good(:,4), params_good(:,2), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition Vs. Growth Rate vs. Transition');
ylabel('Growth Rate');
xlabel('Transition (SD)');

subplot(2,3,6);
plot(params_good(:,3), params_good(:,4), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition Vs. Width');
xlabel('Width');
ylabel('Transition (SD)');



% figure;
% subplot(3,1,1);
% hist(sta_skewness, 25);
% 
% subplot(3,1,2);
% hist(v1_skewness, 25);
% 
% subplot(3,1,3);
% hist(v2_skewness, 25);

return;













