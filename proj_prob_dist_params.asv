function [xspk, xspk_params, pspkx_params] = proj_prob_dist_params(xtrain, xtrain_locator, xtrain_params)
%proj_prob_dist - Probability distributions for projection values
%
% For each training data set, we caclulate the probability of a spike, 
% the probability of a projection without respect to a spike, the 
% probability of a projection when a spike occurs, and the probability of
% a spike given a projection value. This last distribution is the 
% nonlinearity for a neuron.
%
% [xspk, xspk_params, pspkx_params] = proj_prob_dist_params(xtrain, xtrain_locator, xtrain_params)
% ------------------------------------------------------------------------
%
% Input arguments:
%
% xtrain : 1x4 cell array. Each element represents the projection values
% onto one of the training set filters. Includes projections for all
% trials, not just those that lead to a spike.
% 
% xtrain_locator : 1x4 cell array. The locator for each training set. A
% locator is vector whose length is as long as the number of trials in
% a training set. Elements in locator are >= 0, depending on how many 
% spikes occurred during a given trial.
%
% xtrain_params : 1x4 cell array. Each element is an nx2 matrix, with the
% first column the temporal modulation frequency, and the second the
% spectral modulation frequency, of the ripple envelope for each trial.
%
% Output arguments:
%
% xspk : projection values for trials that lead to a spike
% pspkx : 1x4 cell array. An element is the probability of a spike
% given a projection for a training set. 
%
% caa 7/2/09


fprintf('\nRunning proj_prob_dist ...\n');


% Get probability distributions for training data sets

% We need to get: p(spk), p(x), p(x|spk), and p(spk|x)

xbins = -7:7;
xbins = xbins(:);

for i = 1:length(xtrain)

   locator = xtrain_locator{i};
   xprior = xtrain{i};
   xparams = xtrain_params{i}; % mtf params for every projection

% xparams(1:10,:)
% pause

   % Probability of a spike in a given trial
   nspikes = sum(locator);
   pspk_temp = nspikes / length(locator); % probability of a spike

   % Normalize projection values to mean, sd of prior
   mean_prior = mean(xprior);
   std_prior = std(xprior);

   x = (xprior - mean_prior) ./ std_prior;
   xspk = x( locator > 0 ); % projections corresponding to a spike
   xspk_params_temp = xparams(locator>0, :); % mtf params for a spike


   % Probability distributions of prior, posterior projection values

   nx = hist(x, xbins);
   px_temp = nx ./ sum(nx); % p(x)
   px_temp = px_temp(:);

   nxspk = hist(xspk, xbins);
   pxspk_temp = nxspk ./ sum( nxspk ); % p(x|spk)
   pxspk_temp = pxspk_temp(:);

   pspkx_temp = pspk_temp .* pxspk_temp ./ px_temp; % p(spk|x)

% [xbins pspkx_temp]

   % Now map projection values to probabilities, and then map prob's
   % to modulation parameters.

   pspkx_params_temp = zeros(size(xspk));
   for j = 1:length(xspk)
      pspkx_params_temp(j) = interp1( xbins, pspkx_temp, xspk(j) );
% [xspk(j) pspkx_params_temp(j)]
% pause
   end

size(xspk)
size(xspk_params_temp)
size(pspkx_params_temp)
pause

   % pspkx_params_temp should now have the same size as xspk_params

   % Assign output data
   xspk{i} = xspk_params_temp;
   xspk_params{i} = xspk_params_temp;
   pspkx_params{i} = pspkx_params_temp;

   clear pspk_temp px_temp pxspk_temp pspkx_temp xspk_params_temp pspkx_params_temp

end % (for i)

return;

