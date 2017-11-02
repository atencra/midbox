function sta = get_sta_from_locator(locator, stimulus, numtbins)
% get_sta_from_locator - calculate STA from spike train vector and
%    a stimulus matrix
%
% sta = get_sta_from_locator(locator, stimulus, numtbins)
% ----------------------------------------------------
%
% locator : vector of integers, where values greater
%           than one imply a spike, and values of 0
%           imply no spike
%
% stimulus : the entire ripple stimulus envelope file
%            as one matrix. Is contained in a .mat
%            file such as:
%
%           dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8-matrix.mat
%
% sta : spike triggered average. The sta has the
%       same number of frequencies as the stimulus
%       matrix while the number of time bins is 20.
%
% numtbins : how much memory to include in the sta. Default is 20.
% 
% caa 12/22/06
%
% sta = get_sta_from_locator(locator, stimulus, numtbins)

if ( nargin == 2 )
   numtbins = 20;
end

sta = zeros(size(stimulus,1), numtbins);


for i = numtbins:length(locator)
   if ( locator(i) )
      sta = sta + locator(i) * stimulus(:,i-numtbins+1:i);
   end
end


return;





