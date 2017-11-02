function [mid, info] = get_sgi_dmr_rn_sta_mid_fio_info(filestruct)
% get_rippledir_detailed_sta_mid_fio_info - put MID data files into a struct
%    array so data can be analysed
%
% mid = get_sta_mid_fio_data(filestruct)
%
% The folder C:\MATLAB65\work\tatyana\Filters contains all the data from
% Tatyana's maximally informative dimensions analysis. This function goes
% through that folder, extracts all the files for one location, and then
% creates a struct array of the file information, where each element of the
% struct array corresponds to one single unit in the location.
%
% location : a scalar. May be 500, 515, 516, 517, 519, 532, 534, 537, 602,
% 603, 604, 608, 609, 611, 614, or 616.
%
%
% caa 4/11/06



if ( nargin~=1 )
   error('You need one input arg.');
end

% current_directory = pwd;
% 
% if ( isempty(findstr(current_directory, 'C:\MATLAB65\work\tatyana\Filters')) )
%    error('You need to run this function in C:\MATLAB65\work\tatyana\Filters');
% end


% Parameters that Tatyana uses in her programs
nlag_start = 1; % helps give size of RF estimate
nlag_end = 1;   % if 1 then the RF is one image, not multiple images
% nlags = 1;      % multiple images as in vision work

Nbins = 21; % used for input/output function but not RF

Nbins = 15;

Nbins_medium = 15; % used for input/output function but not RF

Nbins_short = 16; %14; % used for nonlinearity to downsample it; we won't do this
                  % The length of the new noninearity is Nbins_short -1;
                  % since we want it equal to 15, we make it 16

Nparts = 4; % used to identify files, not sure what the significance
            % is, though

% Number of time and frequency bins in the receptive fields
Nh = filestruct(1).nh;
Nv = filestruct(1).nv;
nlags = filestruct(1).nlags;

% Assign the contents of filestruct to mid and info
mid = filestruct;
info = filestruct;

for i = 1:length(filestruct)

%    % Assign basic data structure parameters
%    mid(i) = filestruct(i);
%    info(i) = filestruct(i);


   % Create additional fields if they don't exist
   if isfield(filestruct, 'spl')
      mid(i).spl = filestruct(i).spl;
      info(i).spl = filestruct(i).spl;
   else
      mid(i).spl = 9999;
      info(i).spl = 9999;
   end

   if isfield(filestruct, 'sm')
      mid(i).sm = filestruct(i).sm;
      info(i).sm = filestruct(i).sm;
   else
      mid(i).sm = 9999;
      info(i).sm = 9999;
   end

   if isfield(filestruct, 'tm')
      mid(i).tm = filestruct(i).tm;
      info(i).tm = filestruct(i).tm;
   else
      mid(i).tm = 9999;
      info(i).tm = 9999;
   end

   if isfield(filestruct, 'mdb')
      mid(i).mdb = filestruct(i).mdb;
      info(i).mdb = filestruct(i).mdb;
   else
      mid(i).mdb = 9999;
      info(i).mdb = 9999;
   end

%    mid(i).location = filestruct(i).location;
%    mid(i).unit = filestruct(i).unit;
%    mid(i).x0 = filestruct(i).x0;
%    mid(i).nh = filestruct(i).nh;
%    mid(i).nv = filestruct(i).nv;
%    mid(i).nlags = filestruct(i).nlags;
%    mid(i).numtbins = filestruct(i).tbins;
%    mid(i).numfbins = filestruct(i).fbins;
% 
% 
%    info(i).location = filestruct(i).location;
%    info(i).unit = filestruct(i).unit;
%    info(i).x0 = filestruct(i).x0;
%    info(i).nh = filestruct(i).nh;
%    info(i).nv = filestruct(i).nv;
%    info(i).nlags = filestruct(i).nlags;
%    info(i).numtbins = filestruct(i).tbins;
%    info(i).numfbins = filestruct(i).fbins;


   % Process the sta files
   %====================================================================

   sta_filter_files = filestruct(i).rpsta;
   sta_fio_files = filestruct(i).rpx1pxpxt_sta;

   [sta, rpsta, rpx1pxpxt_sta] = ...
      get_sta_filters_fio(sta_filter_files, sta_fio_files, Nh, Nv, nlags, Nbins, Nbins_short, Nparts);

   info(i).sta = sta;

   mid(i).rpsta = rpsta;
   mid(i).rpx1pxpxt_sta = rpx1pxpxt_sta;


   % Process the rpdbest1 files
   %====================================================================

   best1_v1_files = filestruct(i).rpdbest1_v1;
   best1_v2_files = filestruct(i).rpdbest1_v2;

   [rpdbest1_v1, rpdbest1_v2] = ...
      get_best_filters(best1_v1_files, best1_v2_files, Nh, Nv, nlags);

   mid(i).rpdbest1_v1 = rpdbest1_v1;
   mid(i).rpdbest1_v2 = rpdbest1_v2;



   % Process the rpdbest2 files
   %====================================================================

   best2_v1_files = filestruct(i).rpdbest2_v1;
   best2_v2_files = filestruct(i).rpdbest2_v2;

   [rpdbest2_v1, rpdbest2_v2] = ...
      get_best_filters(best2_v1_files, best2_v2_files, Nh, Nv, nlags);

   mid(i).rpdbest2_v1 = rpdbest2_v1;
   mid(i).rpdbest2_v2 = rpdbest2_v2;



   % Process the rpdtest1 files
   %====================================================================

   test1_v1_files = filestruct(i).rpdtest1_v1;
   test1_v2_files = filestruct(i).rpdtest1_v2;
   test1_fio_files = filestruct(i).rpdx1x2px_pxt_1;

   [dummy1, dummy2, dummy3, rpdtest1_v1, rpdtest1_v2, rpdx1x2px_pxt_1] = ...
      get_test_filters_fio(test1_v1_files, test1_v2_files, test1_fio_files, Nh, Nv, nlags, Nbins, Nbins_short, Nparts);

   % assign the data to the struct array
   mid(i).rpdtest1_v1 = rpdtest1_v1;
   mid(i).rpdtest1_v2 = rpdtest1_v2;
   mid(i).rpdx1x2px_pxt_1 = rpdx1x2px_pxt_1;



   % Process the rpdtest2 files
   %====================================================================

   test2_v1_files = filestruct(i).rpdtest2_v1;
   test2_v2_files = filestruct(i).rpdtest2_v2;
   test2_fio_files = filestruct(i).rpdx1x2px_pxt_2;

   [mid1, mid2, mid12, rpdtest2_v1, rpdtest2_v2, rpdx1x2px_pxt_2] = ...
      get_test_filters_fio(test2_v1_files, test2_v2_files, test2_fio_files, Nh, Nv, nlags, Nbins, Nbins_short, Nparts);

   % assign the data to the struct array
   info(i).mid1 = mid1;
   info(i).mid2 = mid2;
   info(i).mid12 = mid12;

   mid(i).rpdtest2_v1 = rpdtest2_v1;
   mid(i).rpdtest2_v2 = rpdtest2_v2;
   mid(i).rpdx1x2px_pxt_2 = rpdx1x2px_pxt_2;

end % (for)

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v_mean, coeff, projection, mtx] = get_auditory_strf(files, Nh, Nv, nlags)

mtx = [];


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
% of the filters could be similar, though with opposite sign. We want to correct
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

k=1;

for i=1:Nparts
    for j=i+1:Nparts
        projection(k) = sum(mtx(:,i).*mtx(:,j)); % find every possible correlation between the 
        k=k+1;                                   % 4 estimated STAs?
    end
end

return;



function [sta, rpsta, rpx1pxpxt_sta] = ...
      get_sta_filters_fio(sta_filter_files, sta_fio_files, Nh, Nv, nlags, Nbins, Nbins_short, Nparts )


   if ( ~isempty(sta_filter_files) )

      % get the spike triggered average
      %====================================================================
      [v_sta, coeff_sta, projection_sta, mtx_sta] = ...
         get_auditory_strf(sta_filter_files, Nh, Nv, nlags);

      rpsta.filter = v_sta;
      rpsta.coeff = coeff_sta;
      rpsta.projection = projection_sta;
      rpsta.mtx = mtx_sta;


      % get the sta information 
      %====================================================================
      [x_mtx, px_mtx, pxspk_mtx, pspk_mtx, pspkx_mtx, info_sta] = ...
         get_sta_information(sta_fio_files, coeff_sta, Nbins, Nparts);

      sta.x = x_mtx;
      sta.px = px_mtx;
      sta.pxspk = pxspk_mtx;
      sta.pspk = pspk_mtx;
      sta.pspkx = pspkx_mtx;
      sta.information = info_sta;



      % the sta nonlinearity
      %====================================================================
      [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = ...
         get_sta_fior(sta_fio_files, coeff_sta, Nbins, Nbins_short, Nparts);

      rpx1pxpxt_sta.mean_firing = mean_firing;
      rpx1pxpxt_sta.x = x_mean;
      rpx1pxpxt_sta.ior_mean = ior_mean;
      rpx1pxpxt_sta.ior_std = ior_std;
      rpx1pxpxt_sta.px_mean = px_mean;
      rpx1pxpxt_sta.px_std = px_std;


   else

      sta.x = [];
      sta.px = [];
      sta.pxspk = [];
      sta.pspk = [];
      sta.pspkx = [];
      sta.information = [];

      rpsta.filter = [];
      rpsta.coeff = [];
      rpsta.projection = [];
      rpsta.mtx = [];

      rpx1pxpxt_sta.mean_firing = [];
      rpx1pxpxt_sta.x = [];
      rpx1pxpxt_sta.ior_mean = [];
      rpx1pxpxt_sta.ior_std = [];
      rpx1pxpxt_sta.px_mean = [];
      rpx1pxpxt_sta.px_std = [];

   end

return;



function [rpdbest_v1, rpdbest_v2] = ...
      get_best_filters(best_v1_files, best_v2_files, Nh, Nv, nlags)


   if ( ~isempty(best_v1_files) )

      % Filter for projection 1 and projection 2
      [v1, coeff1, projection1, mtx1] = ...
         get_auditory_strf(best_v1_files, Nh, Nv, nlags);

      % Filter for projection 2
      [v2, coeff2, projection2, mtx2] = ...
         get_auditory_strf(best_v2_files, Nh, Nv, nlags);

      % Assign the data
      rpdbest_v1.filter = v1;
      rpdbest_v1.coeff = coeff1;
      rpdbest_v1.projection = projection1;
      rpdbest_v1.mtx = mtx1;

      rpdbest_v2.filter = v2;
      rpdbest_v2.coeff = coeff2;
      rpdbest_v2.projection = projection2;
      rpdbest_v2.mtx = mtx2;

   else

      % Assign the data
      rpdbest_v1.filter = [];
      rpdbest_v1.coeff = [];
      rpdbest_v1.projection = [];
      rpdbest_v1.mtx = [];

      rpdbest_v2.filter = [];
      rpdbest_v2.coeff = [];
      rpdbest_v2.projection = [];
      rpdbest_v2.mtx = [];

   end

return;



function [mid1, mid2, mid12, rpdtest_v1, rpdtest_v2, rpdx1x2px_pxt] = ...
      get_test_filters_fio(test_v1_files, test_v2_files, test_fio_files, Nh, Nv, nlags, Nbins, Nbins_short, Nparts )


   if ( ~isempty(test_fio_files) ) % make sure the 2d nonlinearity exists

      % mid1 filter
      [v1, coeff1, projection1, mtx1] = get_auditory_strf(test_v1_files, Nh, Nv, nlags);

      % mid2 filter
      [v2, coeff2, projection2, mtx2] = get_auditory_strf(test_v2_files, Nh, Nv, nlags);


      % mid1 information
      [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, info_mid1] = ...
         get_mid1_information(test_fio_files, coeff1, coeff2, Nbins, Nparts);

      mid1.x = x1_mtx;
      mid1.px = px1_mtx;
      mid1.pxspk = px1spk_mtx;
      mid1.pspk = pspk_mtx;
      mid1.pspkx = pspkx1_mtx;
      mid1.information = info_mid1;


      % mid2 information
      [x2_mtx, px2_mtx, px2spk_mtx, pspk_mtx, pspkx2_mtx, info_mid2] = ...
         get_mid2_information(test_fio_files, coeff1, coeff2, ...
         Nbins, Nparts);

      mid2.x = x2_mtx;
      mid2.px = px2_mtx;
      mid2.pxspk = px2spk_mtx;
      mid2.pspk = pspk_mtx;
      mid2.pspkx = pspkx2_mtx;
      mid2.information = info_mid2;


      % mid1 and mid2 information
      [x1_mtx, x2_mtx, px1x2_mtx, px1x2spk_mtx, pspk_mtx, pspkx1x2_mtx, info_mid12] = ...
         get_mid1_mid2_information(test_fio_files, coeff1, ...
         coeff2, Nbins, Nparts);

      mid12.x1 = x1_mtx;
      mid12.x2 = x2_mtx;
      mid12.px1x2 = px1x2_mtx;
      mid12.px1x2spk = px1x2spk_mtx;
      mid12.pspk = pspk_mtx;
      mid12.pspkx1x2 = pspkx1x2_mtx;
      mid12.information = info_mid12;


      % Nonlinearity For projection v1, v2, and v1&v2 simultaneously
      [x1_mean, ior1_mean, ior1_std, px1_mean, px1_std] = ...
         get_mid_fior_for_filter1(test_fio_files, ...
         coeff1, coeff2, Nbins, Nbins_short, Nparts);

      [x2_mean, ior2_mean, ior2_std, px2_mean, px2_std] = ...
         get_mid_fior_for_filter2(test_fio_files, ...
         coeff1, coeff2, Nbins, Nbins_short, Nparts);

      [x12, y12, ior12, ior12_mean, ior12_std] = ...
         get_mid_fior_for_filter1_and_filter2(test_fio_files, ...
         coeff1, coeff2, Nbins, Nparts);

% imagesc(x12, y12, ior12_mean);
% xlabel('MID1 Proj');
% ylabel('MID2 Proj');
% axis square;
% axis xy;


      % Assign the data
      rpdtest_v1.filter = v1;
      rpdtest_v1.coeff = coeff1;
      rpdtest_v1.projection = projection1;
      rpdtest_v1.mtx = mtx1;

      rpdtest_v2.filter = v2;
      rpdtest_v2.coeff = coeff2;
      rpdtest_v2.projection = projection2;
      rpdtest_v2.mtx = mtx2;

      rpdx1x2px_pxt.x1 = x1_mean;
      rpdx1x2px_pxt.ior1_mean = ior1_mean;
      rpdx1x2px_pxt.ior1_std = ior1_std;
      rpdx1x2px_pxt.px1_mean = px1_mean;
      rpdx1x2px_pxt.px1_std = px1_std;

      rpdx1x2px_pxt.x2 = x2_mean;
      rpdx1x2px_pxt.ior2_mean = ior2_mean;
      rpdx1x2px_pxt.ior2_std = ior2_std;
      rpdx1x2px_pxt.px2_mean = px2_mean;
      rpdx1x2px_pxt.px2_std = px2_std;

      rpdx1x2px_pxt.x12 = x12;
      rpdx1x2px_pxt.y12 = y12;
      rpdx1x2px_pxt.ior12 = ior12;
      rpdx1x2px_pxt.ior12_mean = ior12_mean;
      rpdx1x2px_pxt.ior12_std = ior12_std;

   else

      % assign the information data
      mid1.x = [];
      mid1.px = [];
      mid1.pxspk = [];
      mid1.pspk = [];
      mid1.pspkx = [];
      mid1.information = [];

      mid2.x = [];
      mid2.px = [];
      mid2.pxspk = [];
      mid2.pspk = [];
      mid2.pspkx = [];
      mid2.information = [];

      mid12.x = [];
      mid12.px = [];
      mid12.pxspk = [];
      mid12.pspk = [];
      mid12.pspkx = [];
      mid12.information = [];


      % Assign the filter/nonlinearity data
      rpdtest_v1.filter = [];
      rpdtest_v1.coeff = [];
      rpdtest_v1.projection = [];
      rpdtest_v1.mtx = [];

      rpdtest_v2.filter = [];
      rpdtest_v2.coeff = [];
      rpdtest_v2.projection = [];
      rpdtest_v2.mtx = [];

      rpdx1x2px_pxt.x1 = [];
      rpdx1x2px_pxt.ior1_mean = [];
      rpdx1x2px_pxt.ior1_std = [];
      rpdx1x2px_pxt.px1_mean = [];
      rpdx1x2px_pxt.px1_std = [];

      rpdx1x2px_pxt.x2 = [];
      rpdx1x2px_pxt.ior2_mean = [];
      rpdx1x2px_pxt.ior2_std = [];
      rpdx1x2px_pxt.px2_mean = [];
      rpdx1x2px_pxt.px2_std = [];

      rpdx1x2px_pxt.x12 = [];
      rpdx1x2px_pxt.y12 = [];
      rpdx1x2px_pxt.ior12 = [];
      rpdx1x2px_pxt.ior12_mean = [];
      rpdx1x2px_pxt.ior12_std = [];

   end

return;





function [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = ...
   get_sta_fior(files, coeff, Nbins, Nbins_short, Nparts)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

if ( nargin < 2 | nargin > 5 )
   error('You need 2 to 5 input args.');
end

if ( nargin == 2 )
   Nbins = 21;
   Nbins_short = 14;
   Nparts = 4;
end

if ( nargin == 3 )
   Nbins_short = 14;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx = 0;


x_mtx = [];    
px_mtx = [];    
ior_mtx = [];
% Nbins_short = 14;
mean_firing = 0;



if ( length(coeff) ~= Nparts ) 
    length(coeff)
    coeff
    Nparts
    display('wrong length of coeff vector');
    return
end

for i = 1:Nparts

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'r');

    if (fp==-1)
        display('error opening file');
        display(files{i});
        pause
        return
    end

    x = fread(fp,Nbins,'double');
    px = fread(fp,Nbins,'double');
    ind0 = find(px<minpx);
    xm = sum(x.*px);
    x = x - xm;
    pxt = fread(fp,Nbins,'double');
    pxt(ind0)=0;

    if (abs(sum(px)-1)>0.001) 
        display(sprintf('sum of px=%f',sum(px))); 
        return;
    end

    rbar = fread(fp,1,'double');
    mean_firing = mean_firing + rbar / Nparts;
    Nrep_eff = fread(fp,1,'double');
    fclose(fp);

    x_r = x;
 
    if (coeff(i)==1)
        ior = pxt ./ (px+eps);
    else 
        for j=1:Nbins
            ior(Nbins+1-j) = pxt(j) / (px(j)+eps);
        end
        x_r = -x;
    end

    ior = reshape(ior,Nbins,1);
    px = reshape(px,Nbins,1);
    maxior = max(ior);
    minior = min(ior);

    x_r = reshape(x_r,Nbins,1); %
    ior = reshape(ior,Nbins,1); %
    ior_mtx = [ior_mtx,ior];
    px_mtx = [px_mtx,px];
    x_mtx = [x_mtx,x_r];
    
end % (for i)


if (isempty(ior_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end


ior_mtx = reshape(ior_mtx,Nbins,Nparts);

x0 = min(min(x_mtx));
x1 = max(max(x_mtx));
maxmax = max([abs(x0) x1]);
edges = linspace(-maxmax, maxmax, Nbins_short);

%edges=linspace(x0,x1,Nbins_short);

ior_rescaled = zeros(Nbins_short-1,Nparts);
px_rescaled = zeros(Nbins_short-1,Nparts);
npoints = zeros(Nbins_short-1,Nparts);

for i=1:Nbins_short-1

    for j=1:Nparts

        ind = find( (x_mtx(:,j)>=edges(i)) & (x_mtx(:,j)<edges(i+1)) );

        if ( ~isempty(ind) )
            npoints(i,j) = npoints(i,j) + length(ind);
            px_rescaled(i,j) = px_rescaled(i,j) + sum( px_mtx(ind,j) );
            ior_rescaled(i,j) = ior_rescaled(i,j) + sum( ior_mtx(ind,j) .* px_mtx(ind,j) );
        end

    end

    x_mean(i) = 0.5 * ( edges(i) + edges(i+1) );

end

ior_rescaled = ior_rescaled ./ (px_rescaled+eps);
px_mean = mean(px_rescaled');

ior_mean = mean(ior_rescaled');

ior_std = sqrt(var(ior_rescaled') / Nparts);
px_std = sqrt(var(px_rescaled') / Nparts);    

dx = (max(x_mean)-min(x_mean))/(Nbins_short-1);

px_mean = px_mean ./ dx;
px_std = px_std ./ dx;

return;




function [x1_mean, ior1_mean, ior1_std, px1_mean, px1_std] = ...
   get_mid_fior_for_filter1(files, coeff1, coeff2, Nbins, Nbins_short, Nparts)
%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density


if ( nargin < 3 | nargin > 6 )
   error('You need 3 to 6 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nbins_short = 14;
   Nparts = 4;
end

if ( nargin == 4 )
   Nbins_short = 14;
   Nparts = 4;
end

if ( nargin == 5 )
   Nparts = 4;
end

minpx = 0;


x1_mtx=[];
ior1_mtx=[];
px1_mtx=[];
x2_mtx=[];
ior2_mtx=[];
ior_mtx=[];


if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   display('wrong length of coeff vector');
   return
end


Nok = 0;

for i = 1:Nparts

   %fname = sprintf('%s_%u.dat',fname_first,i);
   fp = fopen(files{i},'r');

   if (fp==-1)
       display('error opening file');
       display(files{i});
       return
   end

   x1 = fread(fp,Nbins,'double');
   x2 = fread(fp,Nbins,'double');
   px = fread(fp,Nbins*Nbins,'double');
   pxt = fread(fp,Nbins*Nbins,'double');
   rbar = fread(fp,1,'double');
   Neff = fread(fp,1,'double');

   fclose(fp);

   ind0 = find(px+eps<minpx);
   pxt(ind0) = 0;
   px = reshape(px,Nbins,Nbins);  %px(x1,x2)
   pxt = reshape(pxt,Nbins,Nbins);
   if (coeff1(i)==1)
       px_r = px;
       pxt_r = pxt;
       x1_r = x1;
   else
       x1_r = -x1;
   end

   if (coeff2(i)==-1)
       x2_r = -x2;
   end

   if (coeff2(i)==1) 
      x2_r = x2; 
   end

   ior_mtx = [ior_mtx,reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
   pxt1 = sum(pxt_r');
   px1 = sum(px_r');
   pxt2 = sum(pxt_r);
   px2 = sum(px_r);
   xm1 = sum(x1_r.*px1');
   xm2 = sum(x2_r.*px2');
   x1_r = x1_r - xm1;
   x2_r = x2_r - xm2;
   ior2 = pxt2 ./ (px2+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
   ior1 = pxt1 ./ (px1+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
   ior1 = reshape(ior1,Nbins,1); % stolbets
   px1 = reshape(px1,Nbins,1); % stolbets
   ior2 = reshape(ior2,Nbins,1); % stolbets

   x1_r = reshape(x1_r,Nbins,1); % stolbets
   x2_r = reshape(x2_r,Nbins,1); % stolbets
   ior1_mtx = [ior1_mtx,ior1];
   ior2_mtx = [ior2_mtx,ior2];
   px1_mtx = [px1_mtx,px1];
   x1_mtx = [x1_mtx,x1_r];
   x2_mtx = [x2_mtx,x2_r];
end

if (isempty(ior_mtx) )
   display('empty ior_mtx in plot_an_ior');
   return;
end

x0 = min(min(x1_mtx));
x1 = max(max(x1_mtx));
maxmax = max([abs(x0) x1]);
edges = linspace(-maxmax, maxmax, Nbins_short);

ior1_rescaled = zeros(Nbins_short-1,Nparts);
px1_rescaled = zeros(Nbins_short-1,Nparts);
npoints = zeros(Nbins_short-1,Nparts);

for i=1:Nbins_short-1
   for j=1:Nparts
       ind = find( (x1_mtx(:,j)>=edges(i)) & (x1_mtx(:,j)<edges(i+1)) );
       if ~isempty(ind)
           npoints(i,j) = 1;%npoints(i,j)+length(ior1_mtx(ind,j));
           ior1_rescaled(i,j) = ior1_rescaled(i,j) + sum(ior1_mtx(ind,j) .* px1_mtx(ind,j));
           px1_rescaled(i,j) = px1_rescaled(i,j) + sum(px1_mtx(ind,j));
       end
   end
   x1_mean(i) = 0.5 * (edges(i) + edges(i+1)); %mean(x1_mtx');
end

nsamples = sum(npoints');
ior1_rescaled = ior1_rescaled./(px1_rescaled+eps);
px1_mean = sum(px1_rescaled')./(nsamples+eps);
px1_std = sqrt((sum(px1_rescaled'.^2) ./ (nsamples+eps)-px1_mean.^2) ./ (nsamples-1));
sumpx1 = sum(px1_mean);

% if abs(sumpx1-1)>0.1
%     sumpx1
%     fname
%     px1_mean
%     px1_rescaled
%     px1_mtx
%     pause
% end

px1_mean = px1_mean ./ sum(px1_mean);
dx = ( max(x1_mean) - min(x1_mean) ) / ( Nbins_short - 1 );

ior1_mean = sum(ior1_rescaled') ./ (nsamples+eps);
ior1_std = sqrt((sum(ior1_rescaled'.^2)./(nsamples+eps)-ior1_mean.^2)./(nsamples-1).*(nsamples+eps));

px1_mean = px1_mean ./ dx;
px1_std = px1_std ./ dx;

return;








function [x2_mean, ior2_mean, ior2_std, px2_mean, px2_std] = ...
   get_mid_fior_for_filter2(files, coeff1, coeff2, Nbins, Nbins_short, Nparts)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

if ( nargin < 3 | nargin > 6 )
   error('You need 3 to 6 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nbins_short = 14;
   Nparts = 4;
end

if ( nargin == 4 )
   Nbins_short = 14;
   Nparts = 4;
end

if ( nargin == 5 )
   Nparts = 4;
end

minpx = 0;


x1_mtx = [];    
ior1_mtx = [];
x2_mtx = [];    
ior2_mtx = [];  
px2_mtx = [];  
ior_mtx = [];


if (length(coeff1)~=Nparts) 
    length(coeff1)
    Nparts
    Nbins
    Nbins_short
    pause
    display('wrong length of coeff vector');
    return
end

for i=1:Nparts

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'r');

    if (fp==-1)
        display('error opening file');
        display(files{i});
        return
    end
    x1=fread(fp,Nbins,'double');
    x2=fread(fp,Nbins,'double');
    px=fread(fp,Nbins*Nbins,'double');
    pxt=fread(fp,Nbins*Nbins,'double');
    px=reshape(px,Nbins,Nbins);  %px(x1,x2)  
    pxt=reshape(pxt,Nbins,Nbins);
    rbar=fread(fp,1,'double');
    Nrep_eff=fread(fp,1,'double');
    fclose(fp);
    ind0=find(px<minpx);
    pxt(ind0)=0;
    if (coeff1(i)==1)
        px_r=px;
        pxt_r=pxt;
        x1_r=x1;
    else 
        x1_r=-x1;
    end
    if (coeff2(i)==-1)
        x2_r=-x2;
    end
    
    if (coeff2(i)==1) x2_r=x2; end
    ior_mtx=[ior_mtx,reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
    pxt1=sum(pxt_r');
    px1=sum(px_r');
    pxt2=sum(pxt_r);
    px2=sum(px_r);
    ior2=pxt2./(px2+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);   
    ior1=pxt1./(px1+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1=reshape(ior1,Nbins,1); % stolbets
    ior2=reshape(ior2,Nbins,1); % stolbets
    xm1=sum(x1_r.*px1');
    xm2=sum(x2_r.*px2');
    x1_r=x1_r-xm1;
    x2_r=x2_r-xm2;
    x1_r=reshape(x1_r,Nbins,1); % stolbets
    x2_r=reshape(x2_r,Nbins,1); % stolbets
    px2=reshape(px2,Nbins,1); % stolbets
    ior1_mtx=[ior1_mtx,ior1];
    px2_mtx=[px2_mtx,px2];
    ior2_mtx=[ior2_mtx,ior2];
    x1_mtx=[x1_mtx,x1_r];
    x2_mtx=[x2_mtx,x2_r];
end
if (isempty(ior_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end

% x0=min(min(x2_mtx));
% x1=max(max(x2_mtx));
% edges=linspace(x0,x1,Nbins_short);

x0 = min(min(x2_mtx));
x1 = max(max(x2_mtx));
maxmax = max([abs(x0) x1]);
edges = linspace(-maxmax, maxmax, Nbins_short);

ior2_rescaled=zeros(Nbins_short-1,Nparts);
px2_rescaled=zeros(Nbins_short-1,Nparts);
npoints=zeros(Nbins_short-1,Nparts);

for i=1:Nbins_short-1
    for j=1:Nparts
        ind=find((x2_mtx(:,j)>=edges(i))&(x2_mtx(:,j)<edges(i+1)));
        if ~isempty(ind)
            npoints(i,j)=1;%npoints(i,j)+length(ior2_mtx(ind,j));
            ior2_rescaled(i,j)=ior2_rescaled(i,j)+sum(ior2_mtx(ind,j).*px2_mtx(ind,j));
            px2_rescaled(i,j)=px2_rescaled(i,j)+sum(px2_mtx(ind,j));
            %        ior2_rescaled(i,j)=ior2_rescaled(i,j)+sum(ior2_mtx(ind,j));
        end
    end
    x2_mean(i)=0.5*(edges(i)+edges(i+1));
end

ior2_rescaled=ior2_rescaled./(px2_rescaled+eps);

nsamples = sum(npoints');

px2_mean = sum(px2_rescaled')./(nsamples+eps);
px2_std = sqrt(var(px2_rescaled')*(Nparts-1)./(nsamples-1)./(nsamples+eps));

sumpx2 = sum(px2_mean);

% if abs(sumpx2-1)>0.1 
%     sumpx2
%     fname
%     px2_mean
%     px2_rescaled
%     px2_mtx
%     pause
% end

px2_mean = px2_mean./sum(px2_mean);
dx = (max(x2_mean)-min(x2_mean))/(Nbins_short-1);

ior2_mean = sum(ior2_rescaled')./(nsamples+eps);
ior2_std = sqrt((sum(ior2_rescaled'.^2)./(nsamples+eps)-ior2_mean.^2)./(nsamples-1));

px2_mean = px2_mean ./ dx;
px2_std = px2_std ./ dx;

% subplot(Nwiny,Nwinx,(row-1)*Nwinx+column)
% errorbar(x2_mean, ior2_mean, ior2_std,'g')
% hold on
% errorbar(x2_mean, px2_mean, px2_std,'k')
% xmax = max(x2_mean)*1.1;
% xmin = min(x2_mean)*1.1;
% xlim([xmin xmax])
% ylim([0 max([0.05,ior2_mean+ior2_std])*1.2])

return;





function [x12, y12, ior_rescaled, ior12_mean, ior12_std] = get_mid_fior_for_filter1_and_filter2(files, coeff1, coeff2, Nbins, Nparts)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density1_mtx=[];
% also properly rotated the two-dimensional matrix, and its labels

if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx = 0;

ior1_mtx = [];
x2_mtx = [];
x1_mtx = [];
ior2_mtx = [];
ior_mtx = [];


% if isempty(varargin)
%     alpha1=0;
%     alpha2=0;
%     beta1=0;
%     beat2=0;
% else
%     alpha1=varargin{1};
%     alpha2=varargin{2};
%     beta1=varargin{3};
%     beta2=varargin{4};
%     non_ortho=varargin{5}
% end


if ( length(coeff1) ~= Nparts )
    length(coeff1)
    display('wrong length of coeff vector');
    return
end

for i = 1:Nparts
    %fname = sprintf('%s_%u.dat',fname_first,i);

    fp = fopen(files{i},'r');

    if ( fp == -1 )
        display('error opening file');
        display(files{i});
        return
    end

    x1 = fread(fp,Nbins,'double');
    x2 = fread(fp,Nbins,'double');
    px = fread(fp,Nbins*Nbins,'double');
    pxt = fread(fp,Nbins*Nbins,'double');

    px = reshape(px,Nbins,Nbins);
    pxt = reshape(pxt,Nbins,Nbins);

    if (coeff1(i)==1)
        px_r = px;
        pxt_r = pxt;
        x1_r = x1;
    else
        for j=1:Nbins
            px_r(Nbins+1-j,:) = px(j,:);
            pxt_r(Nbins+1-j,:) = pxt(j,:);
            x1_r = -x1;
        end
    end

    if (coeff2(i)==-1)
        for j=1:Nbins
            px_r(:,Nbins+1-j) = px(:,j);
            pxt_r(:,Nbins+1-j) = pxt(:,j);
            x2_r = -x2;
        end
    end

    if (coeff2(i)==1) 
      x2_r = x2; 
   end

    ior_mtx = [ior_mtx,reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
    px_mtx = [ior_mtx,reshape(px_r,Nbins*Nbins,1)];
    pxt1 = sum(pxt_r,2);
    px1 = sum(px_r,2);
    pxt2 = sum(pxt_r,1);
    px2 = sum(px_r,1);
    ior2 = pxt2 ./ (px2+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1 = pxt1 ./ (px1+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1 = reshape(ior1,Nbins,1); % stolbets
    ior2 = reshape(ior2,Nbins,1); % stolbets
    x1_r = reshape(x1_r,Nbins,1); % stolbets
    x2_r = reshape(x2_r,Nbins,1); % stolbets
    ior1_mtx = [ior1_mtx,ior1];
    ior2_mtx = [ior2_mtx,ior2];
    x1_mtx = [x1_mtx,x1_r];
    x2_mtx = [x2_mtx,x2_r];

    if ( i==1 )
        rbar=fread(fp,1,'double');
        Nrep_eff=fread(fp,1,'double');
    end

    fclose(fp);

end

if (isempty(ior_mtx) )
   display('empty ior_mtx in plot_an_ior');
   return;
end

x0 = min(min(x1_mtx));
x1 = max(max(x1_mtx));
maxmax = max([abs(x0) x1]);
xedges = linspace(-maxmax, maxmax, Nbins);

y0 = min(min(x2_mtx));
y1 = max(max(x2_mtx));
maxmax = max([abs(y0) y1]);
yedges = linspace(-maxmax, maxmax, Nbins);

ior_rescaled=zeros(Nbins*Nbins,Nparts);

for i=1:Nbins-1
    for k=1:Nbins-1
        for j=1:Nparts
            xind = find((x1_mtx(:,j)>=xedges(i))&(x1_mtx(:,j)<xedges(i+1)));
            yind = find((x2_mtx(:,j)>=yedges(k))&(x2_mtx(:,j)<yedges(k+1)));
            temp = 0;
            sum_px = 0;
            for ii=1:length(xind)
                for jj=1:length(yind)
                    temp = temp+ior_mtx((yind(jj)-1)*Nbins+xind(ii),j).*px_mtx((yind(jj)-1)*Nbins+xind(ii),j);
                    sum_px = sum_px+sum(px_mtx((yind(jj)-1)*Nbins+xind(ii),j));
                end
            end
            if (temp>0)
                temp = temp/sum_px;
            end
            ior_rescaled((k-1)*Nbins+i,j) = ior_rescaled((k-1)*Nbins+i,j)+temp;
        end
    end
    x1_mean(i) = 0.5*(xedges(i)+xedges(i+1));
    x2_mean(i) = 0.5*(yedges(i)+yedges(i+1));
end


ior_mean = mean(ior_rescaled');
ior12_mean = reshape(ior_mean, Nbins, Nbins)';

ior_std = sqrt(var(ior_rescaled') / Nparts);
ior12_std = reshape(ior_std, Nbins, Nbins)';

x12 = x1_mean;
y12 = x2_mean;

% subplot(Nwiny,Nwinx,(row-1)*Nwinx+column);
% imagesc(x12, y12, ior12_mean);
% xlabel('MID1 Proj');
% ylabel('MID2 Proj');
% axis square;
% colorbar;
% ind=find(ior_mean==max(ior_mean));
% indy=mod(ind,Nbins);
% indx=(ind-indy)/Nbins+1;

return;








function [x_mtx, px_mtx, pxspk_mtx, pspk_mtx, pspkx_mtx, information] = ...
   get_sta_information(files, coeff, Nbins, Nparts)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density


if ( nargin < 2 | nargin > 4 )
   error('You need 2 to 4 input args.');
end

if ( nargin == 2 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 3 )
   Nparts = 4;
end

minpx = 0;


x_mtx = [];    
px_mtx = [];
pxspk_mtx = [];
pspkx_mtx = [];
pspk_mtx = [];


if (length(coeff)~=Nparts) 
    length(coeff)
    coeff
    Nparts
    display('wrong length of coeff vector');
    return
end

for i = 1:Nparts

   fp = fopen(files{i},'r');

   if (fp==-1)
       display('error opening file');
       display(files{i});
       pause
       return
   end

   x = fread(fp,Nbins,'double');
   px = fread(fp,Nbins,'double');
   pxspk_pspk = fread(fp,Nbins,'double');
   pspk = fread(fp,1,'double');
   Nrep_eff = fread(fp,1,'double');
   fclose(fp);

   ind0 = find(px < minpx);
   pxspk_pspk(ind0)=0;
   pxspk = pxspk_pspk ./ pspk; % prob of projection given a spike

    xm = sum(x.*px);
    x = x - xm;

    if (abs(sum(px)-1)>0.001) 
        display(sprintf('sum of px=%f',sum(px))); 
        return;
    end

% sum(px)
% sum(pxspk)
% plot(x,px,'r-', x, pxspk, 'k-');
% pause

    x_r = x;
 
   if (coeff(i)==1)
      pspkx = pxspk_pspk ./ (px+eps);
   else 
      for j=1:Nbins
         pspkx(Nbins+1-j) = pxspk_pspk(j) ./ (px(j)+eps);
      end
      x_r = -x;
   end

   x_r = reshape(x_r,Nbins,1);
   px = reshape(px,Nbins,1);
   pxspk = reshape(pxspk,Nbins,1);
   pspkx = reshape(pspkx,Nbins,1);

   x_mtx = [x_mtx x_r];
   px_mtx = [px_mtx px];
   pxspk_mtx = [pxspk_mtx pxspk];
   pspk_mtx = [pspk_mtx pspk];
   pspkx_mtx = [pspkx_mtx pspkx];

end

if ( size(pspkx_mtx,2) ~= Nparts )
    error('empty pxspk_pspk');
end

% size(x_mtx)
% size(px_mtx)
% size(pxspk_mtx)
% size(pspkx_mtx)
% size(pspk_mtx)
% pause


% Calculate information values using a different technique:
pxs1 = pxspk_mtx(:,1);
pxs2 = pxspk_mtx(:,2);
pxs3 = pxspk_mtx(:,3);
pxs4 = pxspk_mtx(:,4);

px1 = px_mtx(:,1);
px2 = px_mtx(:,2);
px3 = px_mtx(:,3);
px4 = px_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

return;


function [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, information] = ...
   get_mid1_information(files, coeff1, coeff2, Nbins, Nparts)

%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density


%       % mid1 information
%       [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, info_mid1] = ...
%          get_mid1_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, ...
%          Nbins_medium, Nparts, Nbins_short);




% Variable name definitions and what they mean:
% 
% pspx1 = P(spike|x1)
% psp = P(spike)
% px1sp = P(x1|spike)
% px1 = P(x1)
% px1x2 = P(x1,x2)
% px1x2sp = P(x1,x2|spike)
% px1sp_psp = P(x1|sp) * P(spike)


if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx = 0;



% initialize variables
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



%Nparts=8;
%Nbins=21;




% Don't really this stuff
% if isempty(varargin)
%     Nparts=8;
%     Nbins=21;
%     Nbins_short=19;
% elseif length(varargin)==1
%     Nparts=8;
%     Nbins_short=19;
%     Nbins=varargin{1};
% elseif length(varargin)==2
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=19;
% elseif length(varargin)==3
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
% elseif length(varargin)==4
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
%     minpx=varargin{4};
% end


if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   display('wrong length of coeff vector');
   return
end

Nok=0;

for i=1:Nparts

   %fname = sprintf('%s_%u.dat',fname_first,i);
   fp = fopen(files{i},'r');

   if (fp==-1)
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

   if (coeff1(i)==1)
       px1x2_r = px1x2;
       px1x2spk_r = px1x2spk;
       x1_r = x1;
   else
       x1_r = -x1;
   end

   if (coeff2(i) == -1)
       x2_r = -x2;
   end

   if (coeff2(i)==1) x2_r=x2; end

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

% sum(px1sp_psp)

   % P(x1|spike)
   px1spk = sum(px1x2spk,2);
   px1spk = px1spk(:);
   px1spk_mtx = [px1spk_mtx px1spk]; 

   % P(x1)
   px1 = sum(px1x2,2);
   px1 = px1(:);
   px1_mtx = [px1_mtx px1]; 

   % P(spike|x1)
   pspkx1 = px1spk_pspk ./ (px1+eps);
   pspkx1 = pspkx1(:);
   pspkx1_mtx = [pspkx1_mtx pspkx1]; 

   xm1 = sum(x1_r .* px1); % expected value: x * p(x) of the projection
   x1_r = x1_r - xm1;
   x1_r = x1_r(:);
   x1_mtx = [x1_mtx x1_r];

end

if (isempty(pspkx1x2_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end

% sum(pspkx1_mtx)
% sum(px1spk_mtx)
% sum(px1_mtx)
% pause

% Calculate information values using a different technique:
pxs1 = px1spk_mtx(:,1);
pxs2 = px1spk_mtx(:,2);
pxs3 = px1spk_mtx(:,3);
pxs4 = px1spk_mtx(:,4);

px1 = px1_mtx(:,1);
px2 = px1_mtx(:,2);
px3 = px1_mtx(:,3);
px4 = px1_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

return;




function [x2_mtx, px2_mtx, px2spk_mtx, pspk_mtx, pspkx2_mtx, information] = ...
   get_mid2_information(files, coeff1, coeff2, Nbins, Nparts)
%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

% Variable name definitions and what they mean:
% 
% pspx1 = P(spike|x1)
% psp = P(spike)
% px1sp = P(x1|spike)
% px1 = P(x1)
% px1x2 = P(x1,x2)
% px1x2sp = P(x1,x2|spike)
% px1sp_psp = P(x1|sp) * P(spike)

if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx=0;



x2_mtx = [];
pspkx2_mtx = [];
px2_mtx = [];
px2spk_mtx = [];
px2spk_pspk_mtx = [];

pspkx1x2_mtx=[];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];


% if isempty(varargin)
%     Nparts=8;
%     Nbins=21;
%     Nbins_short=19;
% elseif length(varargin)==1
%     Nparts=8;
%     Nbins_short=19;
%     Nbins=varargin{1};
% elseif length(varargin)==2
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=19;
% elseif length(varargin)==3
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
% elseif length(varargin)==4
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
%     minpx=varargin{4};
% end



if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   error('wrong length of coeff vector');
end

Nok = 0;

for i=1:Nparts

   fp = fopen(files{i},'r');

   if (fp==-1)
      error(sprintf('error opening file %s', files{i}));
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

   if ( coeff1(i) == -1 )
       x1_r = -x1;
   else
       x1_r = x1;
   end

   if ( coeff2(i) == -1 )
       x2_r = -x2;
   else
       x2_r = x2;
   end

   if (coeff2(i)==1) 
      x2_r = x2; 
   end

   % assign the complete probability distributions
   px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk,Nbins*Nbins,1)]; % P(x1,x2|spike)
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
   pspk_mtx = [pspk_mtx pspk]; % P(spike)
   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)

   % assign the x2 marginal probability distributions

   % P(x2|spike) * P(spike)
   px2spk_pspk = sum(px1x2spk_pspk,1); % sum across columns to get the first marginal
   px2spk_pspk = px2spk_pspk(:);
   px2spk_pspk_mtx = [px2spk_pspk_mtx px2spk_pspk];


   % P(x2|spike)
   px2spk = sum(px1x2spk,1);
   px2spk = px2spk(:);
   px2spk_mtx = [px2spk_mtx px2spk];

   % P(x2)
   px2 = sum(px1x2,1);
   px2 = px2(:);
   px2_mtx = [px2_mtx px2];

   % P(spike|x2)
   pspkx2 = px2spk_pspk ./ (px2+eps);
   pspkx2 = pspkx2(:);
   pspkx2_mtx = [pspkx2_mtx pspkx2];

   xm2 = sum( x2_r .* px2 ); % expected value: x * p(x) of the projection
   x2_r = x2_r - xm2;
   x2_r = x2_r(:);
   x2_mtx = [x2_mtx x2_r];

end


% Calculate information values:
pxs1 = px2spk_mtx(:,1);
pxs2 = px2spk_mtx(:,2);
pxs3 = px2spk_mtx(:,3);
pxs4 = px2spk_mtx(:,4);

px1 = px2_mtx(:,1);
px2 = px2_mtx(:,2);
px3 = px2_mtx(:,3);
px4 = px2_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

return;




function [x1_mtx, x2_mtx, px1x2_mtx, px1x2spk_mtx, pspk_mtx, pspkx1x2_mtx, information] = ...
      get_mid1_mid2_information(files, coeff1, coeff2, Nbins, Nparts)
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within 
% the larger bin, and dividing the two numbers, also introduced dx, and probability distribution 
% is now shown as a probability density1_mtx=[];
% also properly rotated the two-dimensional matrix, and its labels


if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx = 0;

x1_mtx = [];
x2_mtx = [];

pspkx1x2_mtx = [];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];

% if isempty(varargin)
%     alpha1=0;
%     alpha2=0;
%     beta1=0;
%     beat2=0;
% else
%     alpha1=varargin{1};
%     alpha2=varargin{2};
%     beta1=varargin{3};
%     beta2=varargin{4};
%     non_ortho=varargin{5}
% end


if ( length(coeff1) ~= Nparts )
    length(coeff1)
    display('wrong length of coeff vector');
    return
end

for i=1:Nparts

    fp = fopen(files{i},'r');

    if (fp==-1)
        error('error opening file %s', files{i});
        return
    end

% files{i}
% temp = fread(fp,inf,'double');
% size(temp)
% temp
% pause

   x1 = fread(fp,Nbins,'double');
   x2 = fread(fp,Nbins,'double');
% pause
   px1x2 = fread(fp,Nbins*Nbins,'double');
% px1x2(1:10)
% pause
   px1x2spk_pspk = fread(fp,Nbins*Nbins,'double');
   pspk = fread(fp, 1, 'double');

   px1x2 = reshape(px1x2,Nbins,Nbins); % make it a matrix
   px1x2spk_pspk = reshape(px1x2spk_pspk,Nbins,Nbins); % make it a matrix
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike

   fclose(fp);

    if ( coeff1(i)==1 )
        px1x2_r = px1x2;
        px1x2spk_r = px1x2spk;
        px1x2spk_pspk_r = px1x2spk_pspk;
        x1_r = x1;
    else
        for j=1:Nbins
            px1x2_r(Nbins+1-j,:) = px1x2(j,:);
            px1x2spk_pspk_r(Nbins+1-j,:) = px1x2spk_pspk(j,:);
            px1x2spk_r(Nbins+1-j,:) = px1x2spk(j,:);
            x1_r = -x1;
        end
    end

   if ( coeff2(i) == -1 )
      for j=1:Nbins
         px1x2_r(:,Nbins+1-j) = px1x2(:,j);
         px1x2spk_pspk_r(:,Nbins+1-j) = px1x2spk_pspk(:,j);
         px1x2spk_r(:,Nbins+1-j) = px1x2spk(:,j);
         x2_r = -x2;
      end
   end

   if ( coeff2(i)==1 ) 
      x2_r = x2; 
   end

   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk_r./(px1x2_r+eps),Nbins*Nbins,1)];

   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk_r,Nbins*Nbins,1)]; % make it a column vector
   px1x2_mtx = [px1x2_mtx reshape(px1x2_r,Nbins*Nbins,1)]; % make it a column vector
   pspk_mtx = [pspk_mtx pspk];
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk_r,Nbins*Nbins,1)]; % make it a column vector

   x1_r = reshape(x1_r,Nbins,1); % stolbets
   x2_r = reshape(x2_r,Nbins,1); % stolbets

   x1_mtx = [x1_mtx x1_r];
   x2_mtx = [x2_mtx x2_r];

end


% Calculate information values:
pxs1 = px1x2spk_mtx(:,1);
pxs2 = px1x2spk_mtx(:,2);
pxs3 = px1x2spk_mtx(:,3);
pxs4 = px1x2spk_mtx(:,4);

px1 = px1x2_mtx(:,1);
px2 = px1x2_mtx(:,2);
px3 = px1x2_mtx(:,3);
px4 = px1x2_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

return;



















