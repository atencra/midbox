function [v_mean, coeff, projection, mtx] = get_auditory_filter(files, Nh, nlags)
%
% [v_mean, coeff, projection, mtx] = get_auditory_filter(files, Nh, nlags)
%
%fname_first : dat filter filename. 
%              Examples: rpsta_707_1_1x25x20_1_1.dat 
%                     or rpdtest2_prelim_v1_707_1_1x25x20_1_4.dat
%
% nlags : number of time bins, or horizontal bins, in filter. This is always
% equal to 20.
%
% Nh : number of vertical bins in filter, i.e. number of frequencies
%
% varargin : number of test reps. Should be 4.
%
% caa 2/16/09


mtx = [];

Nv = 1;
fsize = Nv*Nh;
Nn = fsize*nlags; % total number of elements in the STRF
Nparts = length(files);

for i=1:Nparts

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'rb');

    if ~(fp== -1)
        [v]=fread(fp,'double');
        v=reshape(v,Nn,1);
        mtx=[mtx,v];
        fclose(fp);
    else v=[];
    end
    
end

if (isempty(mtx) )
    files{1}
    display('empty mtx in plot_a_vector');
    return;
end

mtx = reshape(mtx,Nn,Nparts);
coeff(1) = 1;

% This just corrects a possible pos/neg symmetry issue. It's possible some
% of the filters could be similar, though with similar sign. We want to correct
% for that and make sure all the filters overlap correctly
for i=2:Nparts
    coeff(i)=sign(sum(mtx(:,i).*mtx(:,1)));
    if (coeff(i) == -1)
        mtx(:,i)=mtx(:,i)*coeff(i);
    end
end

v_mean = mean(mtx')';
v_mean = v_mean ./ sqrt(sum(sum(v_mean.*v_mean)));
v_std = sqrt(sum(var(mtx'))*(Nparts-1)/Nn);
v_mean = v_mean ./ v_std;
v_mean = reshape(v_mean,fsize,nlags);


% find every possible correlation between the STAs 
k=1;
for i=1:Nparts
    for j=i+1:Nparts
        projection(k) = sum(mtx(:,i).*mtx(:,j)); 
        k=k+1;                                   % 4 estimated STAs?
    end
end


return;


