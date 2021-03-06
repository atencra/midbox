function plot_fio_curve_fit_examples(fioFit)
%
% plot_fio_curve_fit_examples(fioFit) - data and nonlinearity fits
%
% plot_fio_curve_fit_examples(fioFit)
% ====================================================================
%
% fioFit - struct array holding curve fitting parameters for each 
% nonlinearity. The curve that was fitted was from Hansel and 
% van Vresswijk (2002).
%
% fioFit = get_bandwidth_strf_nonlinearity_hvv_fit(unqbwdb_nonlinearity)
%


for i = 1:2:length(fioFit)

   x = fioFit(i).x;
   fx = fioFit(i).fx;

   xFit = fioFit(i).xFit;
   fxFit = fioFit(i).fxFit;

   fitParams = fioFit(i).fitParams;
   theta = fitParams(2);
   sigma = fitParams(3);

   nmse = fioFit(i).nmse;
   r2 = fioFit(i).r2;

   close all;
   figure;

   subplot(1,2,1);
   x = fioFit(i).x;
   fx = fioFit(i).fx;
   xFit = fioFit(i).xFit;
   fxFit = fioFit(i).fxFit;
   fitParams = fioFit(i).fitParams;
   theta = fitParams(2);
   sigma = fitParams(3);
   nmse = fioFit(i).nmse;
   r2 = fioFit(i).r2;

   hold on;
   plot(xFit, fxFit, '-', 'color', [0.5 0.5 0.5], 'linewidth', 2);
   plot(x, fx, 'ko', 'markerfacecolor', 'k', 'markersize', 3);
   plot([min(x) max(x)], [0 0], 'k-');
   xlim([-7 7]);
   ylimit = get(gca,'ylim');
   set(gca,'ylim', [-0.05*max(fx) 1.05*max(fx)]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   legend('Fit', 'Data', 'location', 'northwest');
   title(sprintf('#%.0f: NMSE=%.2f, R2=%.2f\ntheta=%.1f, sigma=%.1f', ...
   i, nmse, r2, theta, sigma));

   subplot(1,2,2);
   x = fioFit(i+1).x;
   fx = fioFit(i+1).fx;
   xFit = fioFit(i+1).xFit;
   fxFit = fioFit(i+1).fxFit;
   fitParams = fioFit(i+1).fitParams;
   theta = fitParams(2);
   sigma = fitParams(3);
   nmse = fioFit(i+1).nmse;
   r2 = fioFit(i+1).r2;

   hold on;
   plot(xFit, fxFit, '-', 'color', [0.5 0.5 0.5], 'linewidth', 2);
   plot(x, fx, 'ko', 'markerfacecolor', 'k', 'markersize', 3);
   plot([min(x) max(x)], [0 0], 'k-');
   xlim([-7 7]);
   ylimit = get(gca,'ylim');
   set(gca,'ylim', [-0.05*max(fx) 1.05*max(fx)]);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   title(sprintf('#%.0f: NMSE=%.2f, R2=%.2f\ntheta=%.1f, sigma=%.1f', ...
   i+1, nmse, r2, theta, sigma));


   set(gcf,'position', [940 820 321 157]);
   print_mfilename(mfilename);

   pause

end % (for i)

return;
























