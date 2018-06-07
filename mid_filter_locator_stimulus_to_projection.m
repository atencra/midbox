function [xprior, xposterior] = mid_filter_locator_stimulus_to_projection(v, locator, stimulus)
% mid_filter_locator_stimulus_to_projection Project stimulus onto filter
% 
%     [xprior, xposterior] = mid_filter_locator_stimulus_to_projection(v, locator, stimulus)
% 
%     v : NFxNlags receptive field matrix.
% 
%     locator : 1 X Ntrials spike train vector
% 
%     stimulus : NF X Ntrials stimulus matrix
% 
%     xprior : projection values for all stimuli
% 
%     xposterior : projection values corresponding to a spike

fprintf('%s\n', mfilename);


ntrials = length(locator);
[~, nc] = size(v); % # frequencies, # time bins


% Find the prior projection values for all stimuli
xprior = zeros(ntrials, 1);

for i = nc:ntrials

   xprior(i) =  sum( sum ( stimulus( :, i-nc+1:i ) .* v ) ); % inner product

   if ( ~mod(i, 100000) )
      fprintf('i = %.0f\n', i);
   end
end % (for i)


% Find the prior projection values corresponding to a spike
xposterior = xprior( locator > 0 );


return;



