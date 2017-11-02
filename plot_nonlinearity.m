function plot_nonlinearity(xbins, pspk, pspkx, titlestr)
%plot_filters_nonlinearities - display filters and nonlinearities.
%
% plot_nonlinearity(xbins, pspk, pspkx)
% -------------------------------------------------------------------
%
% xbins : values at which projection values were binned. This is the
% abscissa for the nonlinearity.
%
% pspk : 1x4 cell array. Probability of a spike for each of the 4 training
% data sets.
%
% pspkx : 1x4 cell array. Each element represents the nonlinearity for one
% of the training sets.
%
% caa 3/13/09

fprintf('\nRunning plot_filters_nonlinearities ...\n');


if ( nargin < 3 | nargin > 4 )
   error('You need 3 or 4 input args.');
end

if ( nargin ~= 4 )
   titlestr = [];
end






% Plot the nonlinearities
% ----------------------------------------------------------

maxmax = max([ max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

subplot(2,2,1);
hold on;
plot(xbins, pspkx{1}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{1} pspk{1}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
ylabel('P(spk|x)');
title('Train Set 1');

subplot(2,2,2);
hold on;
plot(xbins, pspkx{2}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{2} pspk{2}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
ylabel('P(spk|x)');
title('Train Set 2');

subplot(2,2,3);
hold on;
plot(xbins, pspkx{3}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{3} pspk{3}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
ylabel('P(spk|x)');
title('Train Set 3');

subplot(2,2,4);
hold on;
plot(xbins, pspkx{4}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{4} pspk{4}], 'k--');
xlim([-8 8]);
ylim([0 1.05*maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlabel('Projection (SD)');
ylabel('P(spk|x)');
title('Train Set 4');

titlestr

if (  ~isempty(titlestr) )
   suptitle(sprintf('%s : N = %.0f', titlestr, length(xbins)) );
end

set(gcf, 'position', [207   442   758   474]);


return;