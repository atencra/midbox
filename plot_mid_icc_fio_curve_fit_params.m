function plot_mid_icc_fio_curve_fit_params(params0, params1, nmse0, nmse1)
%plot_mid_icc_fio_curve_fit_params Nonlinearity curve fit param histograms
%
% plot_mid_icc_fio_curve_fit_params(params0, params1, nmse0, nmse1)
% --------------------------------------------------------------------
% params0, params1 : parameters from curve fits to nonlinearities for 
% STA (params0) and MID1 (params1).
%
% params0 and params1 are obtained from plot_nonlinearity_curve_fit.m
%
% params0 and params1 are Nx4 arrays. Each row represent the parameters for
% fits to a nonlinearity, or a neuron. Each column is a parameter. The columns 
% represent the parameters as:
%
% params(x,1) : baseline
% params(x,2) : growth rate
% params(x,3) : width
% params(x,4) : transition
%
% caa 2/18/10

if ( nargin ~= 2 && nargin ~= 4 )
   error('You need 2 or 4 input arguments.');
end


if ( nargin == 4 )
   i0 = find(nmse0 < 0.25);
   params0 = params0(i0,:);

   i1 = find(nmse1 < 0.25);
   params1 = params1(i1,:);
end



close all;

nbins = 11;

figure;


subplot(2,2,1);
edges = linspace(0, 10, nbins);
np = histc(params0(:,2), edges);
hb = bar(edges, np, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlim([edges(1) edges(end)]);
box off;
title('Growth Rate');

subplot(2,2,2);
edges = linspace(-4, 6, nbins);
np = histc(params0(:,4), edges);
hb = bar(edges, np, 'histc');
set(hb, 'facecolor', [0.6 0.6 0.6]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlim([edges(1)-1 edges(end)+1]);
box off;
title('Transition (SD)');

subplot(2,2,3);
[r,p] = corrcoef(params0(:,4), params0(:,2));
plot(params0(:,4), params0(:,2), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title(sprintf('STA Transition Vs. Growth Rate\nr=%.3f, p=%.3f', r(1,2), p(1,2)) );
ylabel('STA Growth Rate');
xlabel('STA Transition (SD)');

suptitle('STA Nonlinearity Curve Fit Params');

print_mfilename(mfilename);


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


