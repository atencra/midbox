function [params, nmse] = plot_hansel_vanvreeswijk_examples
%
% plot_hansel_vanvreeswijk_examples - parametric curve fits
%
% ====================================================================
%
%



   x = [ -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7];
   fx =[0 0.0034 0.0059 0.0113 0.0220 0.0390 0.0643 0.0752 0.0868 ...
         0.1077 0.1329 0.1600 0.1913 0.2482 0.3245];
   f0 = 0.06;

rmParams = [];
saNmse = [];
rmNmse = [];
rmNmse2 = [];

   fx = fx ./ .001; % 
   f0 = f0 ./ .001; % 

   xfit = linspace(min(x), max(x), 100);

%{
   % Fit function from Sachs and Abbas (19??)
   a0 = [1 1 1 1];
   [a, resnorm] = lsqcurvefit(@sachs_abbas_func, a0, x, fx);
   saParams = [saParams; a(:)'];
   saFit = sachs_abbas_func(a, xfit);

   temp = sachs_abbas_func(a, x);
   [gf] = gfit2(fx, temp, '2'); % gets goodness of fit measure, nmse
   saNmse = [saNmse; gf];
%}

   % Fit function from Ringach and Malone (2007)
   b0 = [1 1 1];
   [b, resnorm] = lsqcurvefit(@ringach_malone_func, b0, x, fx);
   rmParams = [rmParams; b(:)'];
   rmFit = ringach_malone_func(b, xfit);
   b

   temp = ringach_malone_func(b, x);
   [gf] = gfit2(fx, temp, '2'); % gets goodness of fit measure, nmse
   rmNmse = [rmNmse; gf];

   % Fit function from Ringach and Malone (2007)
   b0 = [1 1 1];
   fx2 = fx;
   fx2(fx<f0) = f0;
   fx2 = fx2 - f0;
   [b2, resnorm] = lsqcurvefit(@ringach_malone_func, b0, x, fx2);
   rmParams = [rmParams; b2(:)'];
   rmFit2 = ringach_malone_func(b2, xfit);

   temp = ringach_malone_func(b2, x);
   [gf] = gfit2(fx2, temp, '2'); % gets goodness of fit measure, nmse
   rmNmse2 = [rmNmse2; gf];
   b2


   close all;
   figure;

   subplot(2,3,1);
   hold on;
   plot(x, fx, 'ko', 'markerfacecolor', 'k');
   %plot(xfit, saFit, 'c-');
   plot(xfit, rmFit, 'r-');
   plot([-8 8], [f0 f0], 'k-');
   plot([b(2) b(2)], [0 max(fx)], 'k-');
   xlim([-8 8]);
   ylimit = get(gca,'ylim');
   set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   legend('Data', 'HVV', 'location', 'northwest');
   %title(sprintf('#%.0f : SA NMSE = %.3f, RM NMSE = %.3f', saNmse, rmNmse));
   title(sprintf('RM NMSE = %.3f', rmNmse));
   xlim([-8 8]);


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
   title(sprintf('RM NMSE = %.3f', rmNmse2));
   xlim([-8 8]);
   ylimit = get(gca,'ylim');
   set(gca,'ylim', [-0.1*max(fx2) max(ylimit)]);

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
   %suptitle('Hansel and van Vreeswijk (2002) Curves');
   print_mfilename(mfilename);




   close all;
   figure;


   subplot(1,2,1);
   hold on;
   thresh1 = 0;
   thresh2 = 2.5;
   thresh3 = 5;
   amp = 10;
   sigma = 0;
   rmFit1 = ringach_malone_func([amp thresh1 sigma], xfit);
   rmFit2 = ringach_malone_func([amp thresh2 sigma], xfit);
   rmFit3 = ringach_malone_func([amp thresh3 sigma], xfit);
   plot(xfit, rmFit1, 'k-');
   plot(xfit, rmFit2, 'r-');
   plot(xfit, rmFit3, 'b-');
%   plot([-8 8], [f0 f0], 'k-');
   xlim([-8 8]);
   ylim([-2.5 65]);
   ylimit = get(gca,'ylim');
%   set(gca,'ylim', [-0.025*max(fx2) max(ylimit)]);
   set(gca,'xtick', [-5 0 2.5 5], 'xticklabel', [-5 0 2.5 5]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   title(sprintf('Sigma = %.2f', sigma)); 
   legend(sprintf('Thr = %.2f', thresh1), sprintf('Thr = %.2f', thresh2), ...
      sprintf('Thr = %.2f', thresh3), 'Location', 'NorthWest');

   subplot(1,2,2);
   hold on;
   thresh = 2.5;
   amp = 10;
   sigma1 = 0;
   sigma2 = 2;
   sigma3 = 4;
   rmFit1 = ringach_malone_func([amp thresh sigma1], xfit);
   rmFit2 = ringach_malone_func([amp thresh sigma2], xfit);
   rmFit3 = ringach_malone_func([amp thresh sigma3], xfit);
   plot(xfit, rmFit1, 'k-');
   plot(xfit, rmFit2, 'r-');
   plot(xfit, rmFit3, 'b-');
%   plot([-8 8], [f0 f0], 'k-');
   xlim([-8 8]);
   ylim([-2.5 65]);
   ylimit = get(gca,'ylim');
%   set(gca,'ylim', [-0.025*max(fx2) max(ylimit)]);
   set(gca,'xtick', [-5 0 2.5 5], 'xticklabel', [-5 0 2.5 5]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   title(sprintf('Threshold = %.2f', thresh)); 
   legend(sprintf('sigma = %.2f', sigma1), sprintf('sigma = %.2f', sigma2), ...
      sprintf('sigma = %.2f', sigma3), 'Location', 'NorthWest');


   set(gcf,'position', [856   345   324   534]);
   %suptitle('Hansel and van Vreeswijk (2002) Curves');
   print_mfilename(mfilename);

   return;
























