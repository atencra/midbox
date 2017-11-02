function [proj, fio, projinfo] = get_projection_from_filtstr(filtstr, stimulus, locator_mtx)
% get_fio_info_from_filtstr   Projections, nonlinearities, and info from stimulus, spikes, and filters
%
% [proj, fio, projinfo] = get_fio_info_from_filtstr(filtstr, stimulus, locator_mtx)
% ------------------------------------------------------------------------
%
% filtstr : struct holding filters from mid calculations. It is obtained from
% get_filters_from_dat_files.m Filtstr must have length one. Longer lengths
% are not allowed, since they take up too much memory and cause matlab to 
% crash.
%
% stimulus : the entire ripple stimulus envelope file
%            as one matrix. 
%
% locator_mtx : matrix of integers, where values greater
%           than one imply a spike, and values of 0
%           imply no spike. Each row corresponds to an element 
%           in the struct array filtstr.
%
%           You get the locator from a .isk file and using the call
%           locator = locator_from_isk_file(iskfile);
%
%           If this argument is not supplied, or is empty, then the 
%           function will search the current directory for the *-*.isk
%           file corresponding to the element in filtstr which is being 
%           analyzed.
%
% Same as get_fio_info_filtstr except the projection values are returned. This 
% is very memory intensive - so be careful or matlab may not be able to hold
% the output data.
%
% caa 1/25/10

if ( nargin < 2 | nargin > 3 )
   error('You need 2 or 3 input args.');
elseif ( nargin == 2 )
end


if ( nargin < 1 | nargin > 3 )
   error('You need between 1 and 3 input args.');
end

if ( nargin == 1 )
   stimfile = 'D:\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat';
   s = load(stimfile, 'stimulus');
   stimulus = s.stimulus;
   stimulus = single(stimulus);
   locator_mtx = [];
end

if ( nargin == 2 )
   locator_mtx = [];
end

if ( length(filtstr) ~= 1 )
   error('filtstr must be a single struct, not a struct array.');
end


close all;

fraction = [90 92.5 95 97.5 100]; % previously had 80 and 85 - led to bad fits 

for i = 1:length(filtstr)

   location = filtstr.location;
   unit = filtstr.unit;

   x0 = filtstr.x0; %20;
   nh = filtstr.nh;
   nv = filtstr.nv;
   nlags = filtstr.nlags;
   numtbins = filtstr.tbins;
   numfbins = filtstr.fbins;

   index_freq = (x0):(numfbins-1+x0);

   if ( ~isempty(locator_mtx) )
      locator = uint8( locator_mtx(i,:) );
   else
      d = dir(sprintf('*-%.0f.isk', unit));
      if ( length(d) > 1 )
         error('More than file with that unit number.');
      elseif ( length(d) == 0 )
         error('No isk files for unit %.0f in current folder.', unit);
      else
         locator = uint8( locator_from_isk_file(d.name) );
         
      end
   end % (if)

   if ( length(locator) ~= size(stimulus,2) )
      error('Spike train and envelope file have different number of trials.');
   end


   %--------------------------------------------------------------------
   %   STA Filter from MID Code
   %--------------------------------------------------------------------

   fprintf('\ni = %.0f   STA\n', i);

   sta{1} = reshape( filtstr.mtx_sta(:,1), numfbins, numtbins );
   sta{2} = reshape( filtstr.mtx_sta(:,2), numfbins, numtbins );
   sta{3} = reshape( filtstr.mtx_sta(:,3), numfbins, numtbins );
   sta{4} = reshape( filtstr.mtx_sta(:,4), numfbins, numtbins );


   [x0train, x0test, x0train_locator, x0test_locator] = ...
      train_test_projection(sta, locator, stimulus, index_freq);

   [x0bins, pspk, px0, px0spk, pspkx0] = ...
      proj_prob_dist(x0train, x0train_locator); % for the nonlinearity

   [ifrac0_train] = train_test_info_fraction(x0bins, x0train_locator, x0train, fraction); % training information

   [ifrac0_test] = train_test_info_fraction(x0bins, x0test_locator, x0test, fraction); % test information

   [ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_train] = ...
      train_test_info_fraction_mean_std(fraction, ifrac0_train);

   [ifrac0_mn_test, ifrac0_std_test, ifrac0_mtx_test] = ...
      train_test_info_fraction_mean_std(fraction, ifrac0_test);

   [info0_extrap_train, info0_extrap_test] = ...
      train_test_info_extrapolate(fraction, ifrac0_mtx_train, ifrac0_mn_train, ...
      ifrac0_mtx_test, ifrac0_mn_test);

%    plot_filters_nonlinearities(sta, x0bins, pspk, pspkx0, 'STA Analysis');
% 
%    plot_train_test_information(fraction, ifrac0_mtx_train, ...
%       ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_test, ifrac0_mn_test, ...
%       ifrac0_std_test, 'STA Info Analysis');


   % Save STA data
   %---------------------------------------------------------------
   proj.location = location;
   proj.unit = unit;
   proj.x0 = x0;
   proj.nh = nh;
   proj.nv = nv;
   proj.nlags = nlags;
   proj.numtbins = numtbins;
   proj.numfbins = numfbins;
   proj.index_freq = index_freq;
   proj.locator = uint8(locator);

   proj.x0train = x0train;
   proj.x0test = x0test;
   proj.x0train_locator = x0train_locator;
   proj.x0test_locator = x0test_locator;
   proj.x0bins = x0bins;


   projinfo.location = location;
   projinfo.unit = unit;
   projinfo.x0 = x0;
   projinfo.nh = nh;
   projinfo.nv = nv;
   projinfo.nlags = nlags;
   projinfo.numtbins = numtbins;
   projinfo.numfbins = numfbins;
   projinfo.index_freq = index_freq;

   projinfo.fraction = fraction;
   projinfo.ifrac0_train = ifrac0_train;
   projinfo.ifrac0_test = ifrac0_test;
   projinfo.ifrac0_mn_train = ifrac0_mn_train;
   projinfo.ifrac0_std_train = ifrac0_std_train;
   projinfo.ifrac0_mtx_train = ifrac0_mtx_train;
   projinfo.ifrac0_mn_test = ifrac0_mn_test;
   projinfo.ifrac0_std_test = ifrac0_std_test;
   projinfo.ifrac0_mtx_test = ifrac0_mtx_test;
   projinfo.info0_extrap_train = info0_extrap_train;
   projinfo.info0_extrap_test = info0_extrap_test;


   fio.location = location;
   fio.unit = unit;
   fio.x0 = x0;
   fio.nh = nh;
   fio.nv = nv;
   fio.nlags = nlags;
   fio.numtbins = numtbins;
   fio.numfbins = numfbins;
   fio.index_freq = index_freq;

   fio.x0bins = x0bins;
   fio.pspk = pspk;
   fio.px0 = px0;
   fio.px0spk = px0spk;
   fio.pspkx0 = pspkx0;



   %--------------------------------------------------------------------
   %   MID1 Filter from MID Code
   %--------------------------------------------------------------------

   fprintf('\ni = %.0f   MID1\n', i);

   mid1{1} = reshape( filtstr.mtx_v1(:,1), numfbins, numtbins );
   mid1{2} = reshape( filtstr.mtx_v1(:,2), numfbins, numtbins );
   mid1{3} = reshape( filtstr.mtx_v1(:,3), numfbins, numtbins );
   mid1{4} = reshape( filtstr.mtx_v1(:,4), numfbins, numtbins );

   [x1train, x1test, x1train_locator, x1test_locator] = ...
      train_test_projection(mid1, locator, stimulus, index_freq);

   [x1bins, pspk, px1, px1spk, pspkx1] = proj_prob_dist(x1train, x1train_locator); % for the nonlinearity

   [ifrac1_train] = train_test_info_fraction(x1bins, x1train_locator, x1train, fraction); % training information

   [ifrac1_test] = train_test_info_fraction(x1bins, x1test_locator, x1test, fraction); % test information

   [ifrac1_mn_train, ifrac1_std_train, ifrac1_mtx_train] = ...
      train_test_info_fraction_mean_std(fraction, ifrac1_train);

   [ifrac1_mn_test, ifrac1_std_test, ifrac1_mtx_test] = ...
      train_test_info_fraction_mean_std(fraction, ifrac1_test);

   [info1_extrap_train, info1_extrap_test] = ...
      train_test_info_extrapolate(fraction, ifrac1_mtx_train, ifrac1_mn_train, ...
      ifrac1_mtx_test, ifrac1_mn_test);


%    plot_filters_nonlinearities(mid1, x1bins, pspk, pspkx1, 'MID1 Analysis');
% 
%    plot_train_test_information(fraction, ifrac1_mtx_train, ...
%       ifrac1_mn_train, ifrac1_std_train, ifrac1_mtx_test, ifrac1_mn_test, ...
%       ifrac1_std_test, 'MID1 Info Analysis');


   % Save MID1 data
   %---------------------------------------------------------------
   proj.x1train = x1train;
   proj.x1test = x1test;
   proj.x1train_locator = x1train_locator;
   proj.x1test_locator = x1test_locator;
   proj.x1bins = x1bins;


   projinfo.ifrac1_train = ifrac1_train;
   projinfo.ifrac1_test = ifrac1_test;
   projinfo.ifrac1_mn_train = ifrac1_mn_train;
   projinfo.ifrac1_std_train = ifrac1_std_train;
   projinfo.ifrac1_mtx_train = ifrac1_mtx_train;
   projinfo.ifrac1_mn_test = ifrac1_mn_test;
   projinfo.ifrac1_std_test = ifrac1_std_test;
   projinfo.ifrac1_mtx_test = ifrac1_mtx_test;
   projinfo.info1_extrap_train = info1_extrap_train;
   projinfo.info1_extrap_test = info1_extrap_test;


   fio.x1bins = x1bins;
   fio.px1 = px1;
   fio.px1spk = px1spk;
   fio.pspkx1 = pspkx1;


   %--------------------------------------------------------------------
   %   MID2 Filter from MID Code
   %--------------------------------------------------------------------

   fprintf('\ni = %.0f   MID2\n', i);

   mid2{1} = reshape( filtstr.mtx_v2(:,1), numfbins, numtbins );
   mid2{2} = reshape( filtstr.mtx_v2(:,2), numfbins, numtbins );
   mid2{3} = reshape( filtstr.mtx_v2(:,3), numfbins, numtbins );
   mid2{4} = reshape( filtstr.mtx_v2(:,4), numfbins, numtbins );

   [x2train, x2test, x2train_locator, x2test_locator] = ...
      train_test_projection(mid2, locator, stimulus, index_freq);

   [x2bins, pspk, px2, px2spk, pspkx2] = proj_prob_dist(x2train, x2train_locator); % for the nonlinearity

   [ifrac2_train] = train_test_info_fraction(x2bins, x2train_locator, x2train, fraction); % training information

   [ifrac2_test] = train_test_info_fraction(x2bins, x2test_locator, x2test, fraction); % test information

   [ifrac2_mn_train, ifrac2_std_train, ifrac2_mtx_train] = ...
      train_test_info_fraction_mean_std(fraction, ifrac2_train);

   [ifrac2_mn_test, ifrac2_std_test, ifrac2_mtx_test] = ...
      train_test_info_fraction_mean_std(fraction, ifrac2_test);

   [info2_extrap_train, info2_extrap_test] = ...
      train_test_info_extrapolate(fraction, ifrac2_mtx_train, ifrac2_mn_train, ...
      ifrac2_mtx_test, ifrac2_mn_test);


%    plot_filters_nonlinearities(mid2, x2bins, pspk, pspkx2, 'MID2 Analysis');
% 
%    plot_train_test_information(fraction, ifrac2_mtx_train, ...
%       ifrac2_mn_train, ifrac2_std_train, ifrac2_mtx_test, ifrac2_mn_test, ...
%       ifrac2_std_test, 'MID2 Info Analysis');


   % Save MID2 data
   %---------------------------------------------------------------
   proj.x2train = x2train;
   proj.x2test = x2test;
   proj.x2train_locator = x2train_locator;
   proj.x2test_locator = x2test_locator;
   proj.x2bins = x2bins;


   projinfo.ifrac2_train = ifrac2_train;
   projinfo.ifrac2_test = ifrac2_test;
   projinfo.ifrac2_mn_train = ifrac2_mn_train;
   projinfo.ifrac2_std_train = ifrac2_std_train;
   projinfo.ifrac2_mtx_train = ifrac2_mtx_train;
   projinfo.ifrac2_mn_test = ifrac2_mn_test;
   projinfo.ifrac2_std_test = ifrac2_std_test;
   projinfo.ifrac2_mtx_test = ifrac2_mtx_test;
   projinfo.info2_extrap_train = info2_extrap_train;
   projinfo.info2_extrap_test = info2_extrap_test;


   fio.x2bins = x2bins;
   fio.px2 = px2;
   fio.px2spk = px2spk;
   fio.pspkx2 = pspkx2;


   %--------------------------------------------------------------------
   %   MID1 and MID2 Information
   %--------------------------------------------------------------------

   fprintf('\ni = %.0f   MID1 and MID2\n', i);

   [xbins, pspk, px1x2, px1x2spk, pspkx1x2] = proj_prob_dist_2d(x1train_locator, x1train, x2train); % for 2D nonlinearity

   [ifrac12_train] = train_test_info_fraction_2d(x1bins, x1train_locator, x1train, x2train, fraction);

   [ifrac12_test] = train_test_info_fraction_2d(x1bins, x1test_locator, x1test, x2test, fraction);

   [ifrac12_mn_train, ifrac12_std_train, ifrac12_mtx_train] = ...
      train_test_info_fraction_mean_std(fraction, ifrac12_train);

   [ifrac12_mn_test, ifrac12_std_test, ifrac12_mtx_test] = ...
      train_test_info_fraction_mean_std(fraction, ifrac12_test);

   [info12_extrap_train, info12_extrap_test] = ...
      train_test_info_extrapolate(fraction, ifrac12_mtx_train, ifrac12_mn_train, ...
      ifrac12_mtx_test, ifrac12_mn_test);

%    plot_train_test_information(fraction, ifrac12_mtx_train, ...
%       ifrac12_mn_train, ifrac12_std_train, ifrac12_mtx_test, ifrac12_mn_test, ...
%       ifrac12_std_test, 'MID1 and MID2 Info Analysis');


   % Save MID1 and MID2 data
   %---------------------------------------------------------------
   projinfo.ifrac12_train = ifrac12_train;
   projinfo.ifrac12_test = ifrac12_test;
   projinfo.ifrac12_mn_train = ifrac12_mn_train;
   projinfo.ifrac12_std_train = ifrac12_std_train;
   projinfo.ifrac12_mtx_train = ifrac12_mtx_train;
   projinfo.ifrac12_mn_test = ifrac12_mn_test;
   projinfo.ifrac12_std_test = ifrac12_std_test;
   projinfo.ifrac12_mtx_test = ifrac12_mtx_test;
   projinfo.info12_extrap_train = info12_extrap_train;
   projinfo.info12_extrap_test = info12_extrap_test;

   fio.x12bins = xbins;
   fio.px1x2 = px1x2;
   fio.px1x2spk = px1x2spk;
   fio.pspkx1x2 = pspkx1x2;


   % Assign variables for easy display of information results
   %---------------------------------------------------------------
   infoextrap.train0 = info0_extrap_train;
   infoextrap.train1 = info1_extrap_train;
   infoextrap.train2 = info2_extrap_train;
   infoextrap.train12 = info12_extrap_train;

   infoextrap.test0 = info0_extrap_test;
   infoextrap.test1 = info1_extrap_test;
   infoextrap.test2 = info2_extrap_test;
   infoextrap.test12 = info12_extrap_test;

   % fprintf('STA Train Info : %.4f %.4f %.4f %.4f %.4f\n', infoextrap.train0)
   % fprintf('MID1 Train Info : %.4f %.4f %.4f %.4f %.4f\n', infoextrap.train1)
   % fprintf('MID2 Train Info : %.4f %.4f %.4f %.4f %.4f\n', infoextrap.train1)
   % fprintf('MID12 Train Info : %.4f %.4f %.4f %.4f %.4f\n', infoextrap.train12)

   fprintf('\n');
   fprintf('STA Test Info   : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test0);
   fprintf('MID1 Test Info  : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test1);
   fprintf('MID2 Test Info  : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test2);
   fprintf('MID12 Test Info : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test12);
   fprintf('\n');
   fprintf('STA / MID1      : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test0 ./ infoextrap.test1)
   fprintf('MID1 / MID12    : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test1 ./ infoextrap.test12)
   fprintf('\n');

end % (for i)

return;







