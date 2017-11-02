function fio_curve_fit_single(fio)
%fio_curve_fit - Example curve fit to nonlinearity data
%
% fio_curve_fit(fio)
% ====================================================================
% fio : optional input argument. Holds nonlinearity data. Only the first
%       element of fio will be plotted.
%
% If there are no input arguments, then an example nonlinearity is
% plotted, along with the corresponding curve fit. The curve fit is
% from the paper Ringach and Malone (2007).
%
% Craig Atencio


if ( nargin == 0 ) % use a stored nonlinearity
   x = [ -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7];
   fx =[0 0.0034 0.0059 0.0113 0.0220 0.0390 0.0643 0.0752 0.0868 ...
         0.1077 0.1329 0.1600 0.1913 0.2482 0.3245];
   f0 = 0.06;
else % use the input argument
   x = fio(1).x0bins;
   pspkx0 = fio(1).pspkx0;
   fx = double(nanmean([pspkx0{1} pspkx0{2} pspkx0{3} pspkx0{4}], 2));
   pspk = fio(1).pspk;
   f0 = double(nanmean([pspk{1} pspk{2} pspk{3} pspk{4}], 2));
end 


rmParams = [];
saNmse = [];
rmNmse = [];
rmNmse2 = [];

fx = fx ./ .001; % units of spikes per second
f0 = f0 ./ .001; % units of spikes per second

xfit = linspace(min(x), max(x), 100); % x-axis for curve fit

% Use raw nonlinearity data
b0 = [1 1 1];
[b, resnorm] = lsqcurvefit(@ringach_malone_func, b0, x, fx);
rmParams = [rmParams; b(:)'];
rmFit = ringach_malone_func(b, xfit);

temp = ringach_malone_func(b, x);
[gf] = gfit(fx, temp, '2'); % gets goodness of fit measure, nmse
rmNmse = [rmNmse; gf];


% Set nonlinearity values less than the average rate to be equal to the
% average rate. Make the smallest value equal to zero
b0 = [1 1 1];
fx2 = fx;
fx2(fx<f0) = f0;
fx2 = fx2 - f0;
[b2, resnorm] = lsqcurvefit(@ringach_malone_func, b0, x, fx2);
rmParams = [rmParams; b2(:)'];
rmFit2 = ringach_malone_func(b2, xfit);

temp = ringach_malone_func(b2, x);
[gf] = gfit(fx2, temp, '2'); % gets goodness of fit measure, nmse
rmNmse2 = [rmNmse2; gf];


close all;
figure;

subplot(2,3,1);
hold on;
plot(x, fx, 'ko', 'markerfacecolor', 'k');
plot(xfit, rmFit, 'r-');
plot([-8 8], [f0 f0], 'k-');
plot([b(2) b(2)], [0 max(fx)], 'k-');
xlim([-8 8]);
ylimit = get(gca,'ylim');
set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
legend('Data', 'Curve Fit', 'location', 'northwest');
title(sprintf('NMSE = %.3f', rmNmse));
xlim([-8 8]);
ylabel('Firing Rate (sp/s)');


subplot(2,3,2);
hold on;
thresh = b(2);
rmFit1 = ringach_malone_func([b(1) thresh b(3) ], xfit);
rmFit1m1 = ringach_malone_func([b(1) thresh-1 b(3) ], xfit);
rmFit1p1 = ringach_malone_func([b(1) thresh+1 b(3) ], xfit);
plot(xfit, rmFit1m1, 'k-');
plot(xfit, rmFit1, 'r-');
plot(xfit, rmFit1p1, 'b-');
plot([-8 8], [f0 f0], 'k-');
xlim([-8 8]);
ylimit = get(gca,'ylim');
set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
legend(sprintf('Thr - 1 = %.2f', b(2)-1), sprintf('Thr = %.2f', b(2)), ...
   sprintf('Thr + 1 = %.2f', b(2)+1 ), 'Location', 'NorthWest');


subplot(2,3,3);
hold on;
sigma = b(3);
rmFit1m1 = ringach_malone_func([b(1) b(2) max([sigma-1 0])], xfit);
rmFit1 = ringach_malone_func([b(1) b(2) sigma ], xfit);
rmFit1p1 = ringach_malone_func([b(1) b(2) max([sigma+1 0])], xfit);
plot(xfit, rmFit1m1, 'k-');
plot(xfit, rmFit1, 'r-');
plot(xfit, rmFit1p1, 'b-');
plot([-8 8], [f0 f0], 'k-');
xlim([-8 8]);
ylimit = get(gca,'ylim');
set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
legend(sprintf('Sigma - 1 = %.2f', max([sigma-1 0])), ...
       sprintf('Sigma = %.2f', sigma), ...
       sprintf('Sigma + 1 = %.2f', max([sigma+1 0])), 'Location', 'NorthWest');


subplot(2,3,4);
hold on;
plot(x, fx2, 'ko', 'markerfacecolor', 'k');
plot(xfit, rmFit2, 'r-');
plot([-8 8], [0 0], 'k-');
plot([b2(2) b2(2)], [0 max(fx)], 'k-');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title(sprintf('NMSE = %.3f', rmNmse2));
xlim([-8 8]);
ylimit = get(gca,'ylim');
set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
ylabel('Firing Rate (sp/s)');


subplot(2,3,5);
hold on;
thresh = b2(2);
rmFit2 = ringach_malone_func([b(1) thresh b(3) ], xfit);
rmFit2m1 = ringach_malone_func([b(1) thresh-1 b(3) ], xfit);
rmFit2p1 = ringach_malone_func([b(1) thresh+1 b(3) ], xfit);
plot(xfit, rmFit2m1, 'k-');
plot(xfit, rmFit2, 'r-');
plot(xfit, rmFit2p1, 'b-');
plot([-8 8], [min(fx2) min(fx2)], 'k-');
xlim([-8 8]);
set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
legend(sprintf('Thr - 1 = %.2f', thresh-1), ...
       sprintf('Thr = %.2f', thresh), ...
       sprintf('Thr + 1 = %.2f', thresh+1 ), 'Location', 'NorthWest');


subplot(2,3,6);
hold on;
sigma = b2(3);
rmFit2m1 = ringach_malone_func([b2(1) b2(2) max([sigma-1 0])], xfit);
rmFit2 = ringach_malone_func([b2(1) b2(2) sigma ], xfit);
rmFit2p1 = ringach_malone_func([b2(1) b2(2) max([sigma+1 0])], xfit);
plot(xfit, rmFit2m1, 'k-');
plot(xfit, rmFit2, 'r-');
plot(xfit, rmFit2p1, 'b-');
plot([-8 8], [min(fx2) min(fx2)], 'k-');
xlim([-8 8]);
set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
legend(sprintf('Sigma - 1 = %.2f', max([sigma-1 0])), ...
       sprintf('Sigma = %.2f', sigma), ...
       sprintf('Sigma + 1 = %.2f', max([sigma+1 0])), 'Location', 'NorthWest');
set(gcf,'position', [720 417 1066 534]);
print_mfilename(mfilename);






