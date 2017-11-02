function [pspk, px, pxspk, pspkx] = proj_prob_xbins_dist(xbins, xtrain, xtrain_locator)
%proj_prob_dist - Probability distributions for projection values
%
% For each training data set, we caclulate the probability of a spike, 
% the probability of a projection without respect to a spike, the 
% probability of a projection when a spike occurs, and the probability of
% a spike given a projection value. This last distribution is the 
% nonlinearity for a neuron.
%
% [pspk, px, pxspk, pspkx] = proj_prob_xbins_dist(xbins, xtrain, xtrain_locator)
% ------------------------------------------------------------------------
%
% Input arguments:
%
% xbins : the values at which the probability distributions are binned.
% Values of xbins are usually -7:7. These represent normalized units,
% since the projections are normalized relative to the mean and standard
% deviation of the probability of projection without regard to a spike.
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
% caa 9/11/09


fprintf('\nRunning proj_prob_dist ...\n');


% Get probability distributions for training data sets

% We need to get: p(spk), p(x), p(x|spk), and p(spk|x)

% xbins = -7:7;
xbins = xbins(:);

for i = 1:length(xtrain)

   locator = xtrain_locator{i};
   xprior = xtrain{i};

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


   % Assign output data

   pspk{i} = pspk_temp;
   px{i} = px_temp;
   pxspk{i} = pxspk_temp;
   pspkx{i} = pspkx_temp;

end % (for)

return;



