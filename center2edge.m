function edges = center2edge(centers)
% center2edge - convert histogram bin centers to bin edges, with the same bin sizes 
%
% This is useful when histc() is needed though it is
% more natural to list the bin centers.
%
% The vector edges will have the length
%
% length(edges) = length(centers)+1
%
% caa 12/18/06

delta = centers(2)-centers(1);

a = centers(1)-delta/2;

b = centers+delta/2;

edges = [a(:); b(:)];

return;