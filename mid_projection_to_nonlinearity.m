function [xbins, pspk, px, pxspk, pspkx] = mid_projection_to_nonlinearity(xprior, xposterior)
% mid_projection_to_nonlinearity - Probability distributions for projection values
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
% xprior : all projection values.
% 
% xposterior : projection vaulues for spikes.
%
% Output arguments:
%
% xbins : the values at which the probability distributions are binned.
%
% pspk : probability of a spike.
%
% px : normalized probability of a projection value, in units of SD. 
%
% pxspk : probability of a projection value given a spike. 
% 
% pspkx : nonlinearity. Probability of a spike given a projection value. 
%
%


fprintf('%s\n', mfilename);


% Probability of a spike in a given trial
nspikes = length(xposterior);
pspk = nspikes / length(xprior); % probability of a spike

% Normalize projection values to mean, sd of prior
mean_prior = mean(xprior);
std_prior = std(xprior);

x = (xprior - mean_prior) ./ std_prior;
xspk = (xposterior - mean_prior) ./ std_prior;

xprc_low = prctile(x, [0.1]);
xspkprc_low = prctile(xspk, [0.1]);

xprc_hi = prctile(x, [99.9]);
xspkprc_hi = prctile(xspk, [99.9]);

xbins_edges_low = min(min([xprc_low xspkprc_low]));
xbins_edges_hi = max(max([xprc_hi xspkprc_hi]));

xbins_edges = linspace(xbins_edges_low, xbins_edges_hi, 16);
xbins = edge2center(xbins_edges);

nx = hist(x, xbins);
px = nx ./ sum(nx);
px = px(:);

nxspk = hist(xspk, xbins);
pxspk = nxspk ./ sum( nxspk ); % p(x|spk)
pxspk = pxspk(:);
pspkx = pspk .* pxspk ./ px; % p(spk|x)

return;



