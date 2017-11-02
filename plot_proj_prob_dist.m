function [xbins, pspk, px, pxspk, pspkx] = plot_proj_prob_dist(xtrain, xtrain_locator)
%proj_prob_dist - Probability distributions for projection values
%
% For each training data set, we caclulate the probability of a spike, 
% the probability of a projection without respect to a spike, the 
% probability of a projection when a spike occurs, and the probability of
% a spike given a projection value. This last distribution is the 
% nonlinearity for a neuron.
%
% [xbins, pspk, px, pxspk, pspkx] = proj_prob_dist(xtrain, xtrain_locator)
% ------------------------------------------------------------------------
%
% Input arguments:
%
% xtrain : 1x4 cell array. Each element represents the projection values
% onto one of the training set filters.
% 
% xtrain_locator : 1x4 cell array. The locator for each training set. A
% locator is vector whose length is as long as the number of trials in
% a training set. Elements in locator are >= 0, depending on how many 
% spikes occurred during a given trial.
%
% Output arguments:
%
% xbins : the values at which the probability distributions are binned.
% Values of xbins are usually -7:7. These represent normalized units,
% since the projections are normalized relative to the mean and standard
% deviation of the probability of projection without regard to a spike.
%
% pspk : 1x4 cell array. Each element is the probability of a spike 
% during one of the training sets.
%
% px : 1x4 cell array. An element is the probability of a projection 
% value for a training set. 
%
% pxspk : 1x4 cell array. An element is the probability of a projection 
% value given a spike for a training set. 
% 
% pspkx : 1x4 cell array. An element is the probability of a spike
% given a projection for a training set. 
%
%
% This function is set to run for the data in 
%
% file =
%
% 2003-11-24-site8-2400um-50db-dmr1-fs18115-spk-3.isk
% 
% The axes are customized for this data set. Other sets won't come out well
% since the axis limits, tickmarks will be all wrong.
%
% caa 9/11/09


fprintf('\nRunning proj_prob_dist ...\n');

close all;

% Get probability distributions for training data sets

% We need to get: p(spk), p(x), p(x|spk), and p(spk|x)

xbins = -7:7;
xbins = xbins(:);

xedges = xbins + 0.5;
xedges = [xbins(1)-0.5; xedges];

for i = 4:4%length(xtrain)

   locator = xtrain_locator{i};
   xprior = xtrain{i};
   xposterior = xprior( locator > 0 );

   % Probability of a spike in a given trial
   nspikes = sum(locator);
   pspk_temp = nspikes / length(locator); % probability of a spike

   % Normalize projection values to mean, sd of prior
   mean_prior = mean(xprior);
   std_prior = std(xprior);

   x = (xprior - mean_prior) ./ std_prior;
   xspk = x( locator > 0 ); % values corresponding to a spike


   % Probability distributions of prior, posterior projection values

   nx = hist(x, xbins);
   px_temp = nx ./ sum(nx); % p(x)
   px_temp = px_temp(:);

   nxspk = hist(xspk, xbins);
   pxspk_temp = nxspk ./ sum( nxspk ); % p(x|spk)
   pxspk_temp = pxspk_temp(:);

   pspkx_temp = pspk_temp .* pxspk_temp ./ px_temp; % p(spk|x)



   maxmax = max([ max(xprior) max(xposterior) ]);

   [nc_prior,xc_prior] = hist(xprior, 101);
%    nc_prior = nc_prior ./ sum(nc_prior);

   [nc_posterior,xc_posterior] = hist(xposterior, 101);
%    nc_posterior = nc_posterior ./ sum(nc_posterior);

   ymaxmax = max( [ max(nc_prior) max(nc_posterior) ] );
   xmaxmax = max( [ max(xc_prior) max(xc_posterior) ] );

   nrows = 4;
   ncols = 2;

   subplot(nrows,ncols,1);
   hb = bar(xc_prior,nc_prior);
   children = get(gca, 'children');
   set(children, 'facecolor', [0 0 0]);
   xlim(1.05*[-xmaxmax xmaxmax]);
   ylim(1.05*[0 max(nc_prior)]);
   set(gca,'xtick', -6000:3000:6000, 'xticklabel', -6000:3000:6000);
   set(gca,'ytick', 0:80000:160000, 'yticklabel', 0:80000:160000);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   xlabel('X, Projection Value');
   ylabel('Count');
   title('XPrior Raw');

   subplot(nrows,ncols,2);
   hb = bar(xc_posterior,nc_posterior);
   children = get(gca, 'children');
   set(children, 'facecolor', [0 0 0]);
   xlim(1.05*[-xmaxmax xmaxmax]);
   ylim(1.05*[0 max(nc_posterior)]);
   set(gca,'xtick', -6000:3000:6000, 'xticklabel', -6000:3000:6000);
   set(gca,'ytick', 0:500:1000, 'yticklabel', 0:500:1000);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   xlabel('X, Projection Value');
   ylabel('Count');
   title('XPosterior Raw');


   [nc,xc] = hist(x, 101);
%    nc = nc ./ sum(nc);

   [nc_spk,xc_spk] = hist(xspk, 101);
%    nc_spk = nc_spk ./ sum(nc_spk);

   xmaxmax = max( [ max(xc) abs(min(xc)) max(xc_spk) abs(min(xc_spk)) ] );

max(nc)

   subplot(nrows,ncols,3);
   hb = bar(xc,nc);
   children = get(gca, 'children');
   set(children, 'facecolor', [0 0 0]);
%    xlim(1.05*[-xmaxmax xmaxmax]);
   ylim(1.05*[0 max(nc)]);
   xlim([-8 8]);
   set(gca,'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
   set(gca,'ytick', 0:80000:160000, 'yticklabel', 0:80000:160000);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   xlabel('Similarity with STRF (SD)');
   ylabel('Count');
   title('XPrior Norm');


   subplot(nrows,ncols,4);
   hb = bar(xc_spk,nc_spk);
   children = get(gca, 'children');
   set(children, 'facecolor', [0 0 0]);
%    xlim(1.05*[-xmaxmax xmaxmax]);
   ylim(1.05*[0 max(nc_spk)]);
   xlim([-8 8]);
   set(gca,'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
   set(gca,'ytick', 0:500:1000, 'yticklabel', 0:500:1000);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   xlabel('Similarity with STRF (SD)');
   ylabel('Count');
   title('XPost Norm');


   ymaxmax = max([ max(px_temp) max(pxspk_temp) ]);

   subplot(nrows,ncols,5);
   plot(xbins, px_temp, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   xlim([-8 8]);
   ylim(1.05*[0 ymaxmax]);
   set(gca,'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
   set(gca,'ytick', 0:0.4:0.8, 'yticklabel', 0:0.4:0.8);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   xlabel('Similarity with STRF (SD)');
   ylabel('Prob(Projection)');
   title('P(X)');

   subplot(nrows,ncols,6);
   plot(xbins, pxspk_temp, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
   xlim([-8 8]);
   ylim(1.05*[0 ymaxmax]);
   set(gca,'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
   set(gca,'ytick', 0:0.4:0.8, 'yticklabel', 0:0.4:0.8);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   xlabel('Similarity with STRF (SD)');
   ylabel('Prob(Projection|Spike)');
   title('P(X|Spike)');


   ymaxmax = max([ max(px_temp) max(pxspk_temp) ]);

   subplot(nrows,ncols,7);
   hold on;
   plot(xbins, px_temp, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
   plot(xbins, pxspk_temp, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);
   xlim([-8 8]);
   ylim(1.05*[0 ymaxmax]);
   set(gca,'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
   set(gca,'ytick', 0:0.4:0.8, 'yticklabel', 0:0.4:0.8);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box off;
   ylabel('Prob(Projection)');
   legend('P(X)', 'P(X|Spike)');



   subplot(nrows,ncols,8);
   hold on;
   plot(xbins, pspkx_temp, 'ro-', 'markerfacecolor', 'r', 'markersize', 2);

%    [ax, h1, h2]=plotyy(xbins, pspkx_temp, xbins, 1000*pspkx_temp);
%    set(ax(1), 'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
%    set(ax(1), 'xlim', [-8 8], 'ylim', 1.05*[0 max(pspkx_temp)]);
%    set(ax(1), 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    set(h1, 'color', 'r', 'marker', 'o', 'markersize', 2, 'markerfacecolor', 'r', 'linestyle', '-');
% 
%    set(ax(2), 'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
%    set(ax(2), 'xlim', [-8 8], 'ylim', 1.05*[0 max(1000*pspkx_temp)]);
%    set(ax(2), 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    set(h2, 'color', 'r', 'marker', 'o', 'markersize', 2, 'markerfacecolor', 'r', 'linestyle', '-');

%    set(ax, 'color','r', 'markerfacecolor', 'r', 'markersize', 2);
   xlim([-8 8]);
   ylim(1.05*[0 max(pspkx_temp)]);
   set(gca,'xtick', -7.5:3.75:7.5, 'xticklabel', -7.5:3.75:7.5);
   set(gca,'ytick', 0:0.075:0.15, 'yticklabel', 0:0.075:0.15);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   box on;
   xlabel('Similarity with STRF (SD)');
   ylabel('Prob(Spike| Projection)');
   title('Nonlinearity');

   set(gcf,'position', [258    205   806   732]);


   % Assign output data

   pspk{i} = pspk_temp;
   px{i} = px_temp;
   pxspk{i} = pxspk_temp;
   pspkx{i} = pspkx_temp;

end % (for)


return;


% Plot the nonlinearities
% ----------------------------------------------------------

maxmax = max([max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

subplot(2,4,5);
hold on;
plot(xbins, pspkx{1}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{1} pspk{1}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
xlabel('Projection (SD)');
ylabel('P(spk|x)');

subplot(2,4,6);
hold on;
plot(xbins, pspkx{2}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{2} pspk{2}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

subplot(2,4,7);
hold on;
plot(xbins, pspkx{3}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{3} pspk{3}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

subplot(2,4,8);
hold on;
plot(xbins, pspkx{4}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{4} pspk{4}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

if (  ~isempty(titlestr) )
   suptitle(titlestr);
end

set(gcf, 'position', [293 267 1095 571]);

