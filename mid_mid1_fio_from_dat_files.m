function [fio] = mid_mid1_fio_from_dat_files(files, coeff1, Nbins)
%function [fio] = get_dat_mid1_fio(files, coeff1, coeff2, Nbins, Nparts, Nbins_short)
%
%
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
% larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
%
%
% files : cell array of four file names. Each element is a string, and is
% the full path name to a .dat file holding the nonlinearity data.
%
%    Examples:
%
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_1.dat          
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_2.dat          
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_3.dat          
%        rpdx1x2px_pxt_1_707_1_1x25x20_1_4.dat        
%
% coeff1, coeff2 : coefficient vector which tells the polarity of the 
% mid1 and mid2 filters
%
% Npart = 4 : num of test reps in mid code. Should be 4.
%
% Nbins = 15 : num bins used to obtain projection values in mid code. The
% nonlinearity will then have length = 13.
%
% Nbins_short = 14 : num bins used to get nonlinearity values
%
%
% If you call get_data_sta_fio.m as
%
% [fio] = get_dat_mid1_fio(files, coeff1, coeff2);
%
% then Nparts = 4, Nbins = 21, and Nbins_short = 14 by default.
%
% fio is a struct that contains the nonlinearity results. fio has the
% following structure:
%
% fio.nbins = Nbins;
% fio.nbins_short = Nbins_short;
% fio.px1x2_mtx = px1x2_mtx;                 P(x1,x2)
% fio.px1x2spk_mtx = px1x2spk_mtx;           P(x1,x2|spk)
% fio.px1x2spk_pspk_mtx = px1x2spk_pspk_mtx; P(x1,x2|spk) * P(spk)
% fio.pspk_mtx = pspk_mtx;                   P(spk)
% fio.pspkx1x2_mtx = pspkx1x2_mtx;           P(spk|x1,x2)
% fio.x1_mtx = x1_mtx;                       x1
% fio.px1_mtx = px1_mtx;                     P(x1)
% fio.pspkx1_mtx = pspkx1_mtx;               P(spk|x1)
% fio.px1spk_mtx = px1spk_mtx;               P(x1|spk)
% fio.px1spk_pspk_mtx = px1spk_pspk_mtx;     P(x1|spk) * P(spk)
% fio.x1_mean = x1_mean;                     x1
% fio.pspk_mean = pspk_mean;                 P(spk)
% fio.px1_mean = px1_mean;                   P(x1)
% fio.px1_std = px1_std;                     P(x1)
% fio.px1spk_mean = px1spk_mean;             P(x1|spk)
% fio.px1spk_std = px1spk_std;               P(x1|spk)
% fio.pspkx1_mean = pspkx1_mean;             P(spk|x1)
% fio.pspkx1_std = pspkx1_std;               P(spk|x1)
% fio.info;
%
% caa 2/17/09


x1_mtx = [];
pspkx1_mtx = [];
px1_mtx = [];
px1spk_mtx = [];
px1spk_pspk_mtx = [];

pspkx1x2_mtx=[];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];


minpx = 0; % probability can't be lower than 0

% The default for Nparts below was originally 8 in Tanya's original code.
if nargin == 2
    Nbins = 15;
end

if ( length(coeff1)~= length(files) ) 
    error('wrong length of coeff vector');
end

Nparts = length(files);


for i = 1:Nparts

   %fname = sprintf('%s_%u.dat',fname_first,i);
   fp = fopen(files{i},'r');

   if ( fp == -1 )
      display('error opening file');
      display(files{i});
      return
   end

   x1 = fread(fp, Nbins, 'double'); % first filter projection values
   x2 = fread(fp, Nbins, 'double'); % second filter projection values
   px1x2 = fread(fp, Nbins*Nbins, 'double'); % prior distr
   px1x2spk_pspk = fread(fp, Nbins*Nbins, 'double'); % really rbar * pxt
   pspk = fread(fp, 1, 'double');
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike
   Neff = fread(fp, 1, 'double');

   fclose(fp);

   % Get rid of machine error values
   ind0 = find(px1x2+eps<minpx);
   px1x2spk(ind0) = 0;

   % Reshape into a matrix so that marginals can be calculated
   px1x2 = reshape(px1x2, Nbins, Nbins);  %p(x1,x2)
   px1x2spk = reshape(px1x2spk, Nbins, Nbins);
   px1x2spk_pspk = reshape(px1x2spk_pspk, Nbins, Nbins);


   if ( coeff1(i) == 1 )
       px1x2_r = px1x2;
       px1x2spk_r = px1x2spk;
       x1_r = x1;
   else
       x1_r = -x1;
   end


   % assign the complete probability distributions
   px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk,Nbins*Nbins,1)]; % P(x1,x2|spike)
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
   pspk_mtx = [pspk_mtx pspk]; % P(spike)
   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)

   % assign the x1 marginal probability distributions

   % P(x1|spike) * P(spike)
   px1spk_pspk = sum(px1x2spk_pspk,2); % sum across columns to get the first marginal
   px1spk_pspk = px1spk_pspk(:);
   px1spk_pspk_mtx = [px1spk_pspk_mtx px1spk_pspk];

   % P(x1|spike)
   px1spk = sum(px1x2spk,2);
   px1spk = px1spk(:);
   px1spk_mtx = [px1spk_mtx px1spk]; 

   % P(x1)
   px1 = sum(px1x2,2);
   px1 = px1(:);
   px1_mtx = [px1_mtx px1]; 

% %    ior1 = px1sp_psp ./ (px1+eps);
% %    ior1 = ior1(:);
% %    ior1_mtx = [ior1_mtx ior1];


   % P(spike|x1)
   pspkx1 = px1spk_pspk ./ (px1+eps);
   pspkx1 = pspkx1(:);
   pspkx1_mtx = [pspkx1_mtx pspkx1]; 

   xm1 = sum(x1_r .* px1); % expected value: x * p(x) of the projection
   x1_r = x1_r - xm1; % center projections on 0
   x1_r = x1_r(:);
   x1_mtx = [x1_mtx x1_r];

end % (for i)


if ( isempty(pspkx1x2_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end



fio.nbins = Nbins;
fio.pspk_mtx = pspk_mtx;
fio.x1_mtx = x1_mtx;
fio.px1_mtx = px1_mtx;
fio.px1spk_mtx = px1spk_mtx;
fio.pspkx1_mtx = pspkx1_mtx;

return;




