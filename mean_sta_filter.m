function [v_mean, coeff, projection, mtx] = mean_sta_filter(sta)
%mean_sta_filter - calculate mean training set STA
%
%[v_mean, coeff, projection, mtx] = mean_sta_filter(sta)
%--------------------------------------------------------------
%
% sta : 1x4 cell array of STAs, one for each training set.
%
% v_mean : mean of the 4 training set STAs, in units of SD.
%
% coeff : vector describing if the training set STA share the same
% polarity as the first training set STA. Values are 1 or -1.
%
% projection : vector, length 6. Correlation between all the STAs.
%
% mtx : matrix of the 4 training set STA. Each column is an STA, and 
% corresponds to the STA contained in the STA input cell array.
%
% This function is modeled after get_auditory_filter.m, except for that
% function the inputs were filters from the MID analysis. Here they are
% previously calculated STAs.
%
% caa 3/3/09



mtx = [];

Nparts = length(sta);

temp = sta{1};

numfbins = size(temp,1);
numtbins = size(temp,2);


fsize = numfbins * numtbins;
Nn = fsize * 1; % total number of elements in the STRF

for i = 1:Nparts
   mtx = [mtx reshape(sta{i},Nn,1) ];
end % (for i)

if (isempty(mtx) )
    error('empty mtx in plot_a_vector');
end

mtx = reshape(mtx, Nn, Nparts); % I don't think this does anything

coeff(1) = 1;

for i = 2:Nparts
    coeff(i) = sign(sum(mtx(:,i).*mtx(:,1)));
    if ( coeff(i) == -1 )
        mtx(:,i)=mtx(:,i)*coeff(i);
    end
end

v_mean = mean(mtx,2);

% v_mean = v_mean ./ sqrt(sum(sum(v_mean.*v_mean))); % Now it's a unit vector

% sqrt( sum(var(mtx') ) )
% 
% (Nparts-1)/Nn

v_std = sqrt( sum(var(mtx,0,2) ) * (Nparts-1)/Nn);

v_mean = v_mean ./ v_std;
v_mean = reshape( v_mean, fsize, 1 );

cm = max([abs(min(min(v_mean))),max(max(v_mean))]);

v_mean = reshape(v_mean(:,1), numfbins, numtbins);

k = 1;
for i = 1:Nparts
    for j = i+1:Nparts
        projection(k) = sum( mtx(:,i).*mtx(:,j) );
        k = k+1;
    end
end

return;


