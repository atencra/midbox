function [fio] = mid_sta_fio_from_dat_files(files, coeff, Nbins)
%function [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] =
%get_dat_sta_fio(folder, location, cellnum, nv, nh, coeff, varargin)
%
%[mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = ...
%     get_dat_sta_fio(folder, location, cellnum, coeff, Nparts, Nbins, Nbins_short, minpx)
%
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
% larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
%
%
% The function assumes that you have files of the form:
%
%     rpx1pxpxt_sta_707_1_1x25x20_1_1.dat            
%     rpx1pxpxt_sta_707_1_1x25x20_1_2.dat            
%     rpx1pxpxt_sta_707_1_1x25x20_1_3.dat            
%     rpx1pxpxt_sta_707_1_1x25x20_1_4.dat     
%
% Parts of the file name are variable, such as 707, 25 , 20. These are specified by
% input arguments.
%
% folder : location of data. Example: 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\';
%
% location : Tanya's site identifier. Examples: 707, 708, 709, etc.
%
% cellnum : Tanya's global cell identifier. Examples: 1, 2, 3, ...
%
% coeffs : coefficient vector which tells the polarity of the filters
%
% Npart = 4 : num of test reps in mid code. Should be 4.
%
% Nbins = 21 : num bins used to obtain projection values in mid code
%
% Nbins_short = 14 : num bins used to get nonlinearity values
%
%
% If you call get_data_sta_fio.m as
%
% [fio] = get_dat_sta_fio(folder, location, cellnum, nv, nh, coeff)
%
% then Nparts = 4, Nbins = 21, and Nbins_short = 14 by default.
%
% caa 2/17/09

% Parameters that Tatyana uses in her programs
% Nbins = 21; % used for input/output function but not RF
% Nbins_short = 14; % used for input/output function but not RF
% Nparts = 4; % used to identify files, not sure what the significance
%             % is, though



%     files{1} = rpx1pxpxt_sta_707_1_1x25x20_1_1.dat            
%     files{2} = rpx1pxpxt_sta_707_1_1x25x20_1_2.dat            
%     files{3} = rpx1pxpxt_sta_707_1_1x25x20_1_3.dat            
%     files{4} = rpx1pxpxt_sta_707_1_1x25x20_1_4.dat     


x_mtx = [];
px_mtx = [];
ior_mtx = [];
mean_firing = 0;
Nparts = length(files);
minpx = 0;

% The default for Nparts below was originally 8 in Tanya's original code.
if nargin == 2
    Nbins = 15;
end


if (length(coeff)~=Nparts) 
    length(coeff)
    coeff
    Nparts
    display('wrong length of coeff vector');
    return
end

if ( Nparts ~= 4 )
   error('For time being, Nparts must equal 4.');
end


pspk_mtx = [];
pxspk_mtx = [];
pspkx_mtx = [];

for i=1:Nparts

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'r');

    if (fp==-1)
        display('error opening file');
        display(files{i});
        pause
        return
    end


    x = fread(fp,Nbins,'double'); % projection value
    px = fread(fp,Nbins,'double'); % prob of a projection, p(x)

    ind0 = find( px < minpx );
    xm = sum( x.*px );
    x = x - xm; % center the projection values to 0

    pxt = fread(fp,Nbins,'double');
    pxt(ind0) = 0;

    pxspk_pspk = pxt; % p(x|spk) * p(spk)


    rbar = fread(fp, 1, 'double');
    pspk = rbar; % p(spk)
    pxspk = pxspk_pspk ./ pspk; % p(x|spk)

    mean_firing = mean_firing + rbar/Nparts;
    mean_pspk = mean_firing;

    Nrep_eff = fread(fp,1,'double');
    fclose(fp);

   x_r = x;

   if ( coeff(i)==1 )
      ior = pxt ./ (px+eps);
      pspkx = ior;
   else 
      for j = 1:Nbins
         ior(Nbins+1-j) = pxt(j) / (px(j)+eps);
      end
      pspkx = ior;
      x_r = -x;
   end

   x_r = reshape(x_r,Nbins,1); % projection value
   px = reshape(px,Nbins,1);
   ior = reshape(ior,Nbins,1); % p(spk|x)

   maxior = max(ior);
   minior = min(ior);

   ior_mtx = [ior_mtx ior]; % matrix of p(spk|x)
   x_mtx = [x_mtx x_r]; % projection matrix
   px_mtx = [px_mtx px];

   pxspk = reshape(pxspk, Nbins,1);
   pspkx = reshape(pspkx, Nbins,1);

   pspk_mtx = [pspk_mtx pspk];
   pxspk_mtx = [pxspk_mtx pxspk];
   pspkx_mtx = [pspkx_mtx pspkx];
    
end


if (isempty(ior_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end


ior_mtx = reshape(ior_mtx,Nbins,Nparts);


fio.nbins = Nbins;
fio.pspk_mtx = pspk_mtx;
fio.pspk_mean = mean_pspk;
fio.x_mtx = x_mtx;
fio.px_mtx = px_mtx;
fio.pxspk_mtx = pxspk_mtx;
fio.pspkx_mtx = pspkx_mtx;
fio.ior_mtx = ior_mtx;

return;











