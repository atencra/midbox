function [fio1] = get_fio_info_from_projection(proj)
% get_fio_info_from_projection   Nonlinearities and info from projections
%
% Calculates values for one unit, for multiple number of bins in nonlinearities
%
% [fio1] = get_fio_info_from_projection(proj)
% ------------------------------------------------------------------------
%
%
%
% caa 1/25/10

if ( nargin ~= 1 )
   error('You need 1 input arg.');
end


if ( length(proj) ~= 1 )
   error('proj must be a single struct, not a struct array.');
end


close all;

fraction = [90 92.5 95 97.5 100]; % previously had 80 and 85 - led to bad fits 

x0train = proj.x0train;
x0test = proj.x0test;
x0train_locator = proj.x0train_locator;
x0test_locator = proj.x0test_locator;

x1train = proj.x1train;
x1test = proj.x1test;
x1train_locator = proj.x1train_locator;
x1test_locator = proj.x1test_locator;

x2train = proj.x2train;
x2test = proj.x2test;
x2train_locator = proj.x2train_locator;
x2test_locator = proj.x2test_locator;


nbins = [5 7 9 11 13 15 17 19 21 23 25 27 29 31];

close all;

fio1 = struct(...
'location', [], ...
'unit',     [], ...
'xbins',    [], ...
'pspk',     [], ...
'px0',      [], ...
'px0spk',   [], ...
'pspkx0',   [], ...
'px1',      [], ...
'px1spk',   [], ...
'pspkx1',   [], ...
'px2',      [], ...
'px2spk',   [], ...
'pspkx2',   []);

location = proj.location;
unit = proj.unit;

for i = 1:length(nbins)

   %--------------------------------------------------------------------
   %   STA Filter from MID Code
   %--------------------------------------------------------------------

   fprintf('\nSTA\n');

   x0bins = linspace( -7, 7, nbins(i) );

   [pspk, px0, px0spk, pspkx0] = proj_prob_xbins_dist(x0bins, x0train, x0train_locator);

%    [ifrac0_train] = train_test_info_fraction(x0bins, x0train_locator, x0train, fraction); % training information
% 
%    [ifrac0_test] = train_test_info_fraction(x0bins, x0test_locator, x0test, fraction); % test information
% 
%    [ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_train] = ...
%       train_test_info_fraction_mean_std(fraction, ifrac0_train);
% 
%    [ifrac0_mn_test, ifrac0_std_test, ifrac0_mtx_test] = ...
%       train_test_info_fraction_mean_std(fraction, ifrac0_test);
% 
%    [info0_extrap_train, info0_extrap_test] = ...
%       train_test_info_extrapolate(fraction, ifrac0_mtx_train, ifrac0_mn_train, ...
%       ifrac0_mtx_test, ifrac0_mn_test);

%    plot_nonlinearity(x0bins, pspk, pspkx0, 'STA Analysis');


%    plot_train_test_information(fraction, ifrac0_mtx_train, ...
%       ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_test, ifrac0_mn_test, ...
%       ifrac0_std_test, 'STA Info Analysis');


%    % Save STA data
%    %---------------------------------------------------------------
   xbins{i} = x0bins;
   px_sta{i} = px0;
   pxspk_sta{i} = px0spk;
   pspkx_sta{i} = pspkx0;



   %--------------------------------------------------------------------
   %   MID1 Filter from MID Code
   %--------------------------------------------------------------------

   fprintf('\nMID1\n');

   [pspk, px1, px1spk, pspkx1] = proj_prob_xbins_dist(x0bins, x1train, x1train_locator);

%    [ifrac1_train] = train_test_info_fraction(x1bins, x1train_locator, x1train, fraction); % training information
% 
%    [ifrac1_test] = train_test_info_fraction(x1bins, x1test_locator, x1test, fraction); % test information
% 
%    [ifrac1_mn_train, ifrac1_std_train, ifrac1_mtx_train] = ...
%       train_test_info_fraction_mean_std(fraction, ifrac1_train);
% 
%    [ifrac1_mn_test, ifrac1_std_test, ifrac1_mtx_test] = ...
%       train_test_info_fraction_mean_std(fraction, ifrac1_test);
% 
%    [info1_extrap_train, info1_extrap_test] = ...
%       train_test_info_extrapolate(fraction, ifrac1_mtx_train, ifrac1_mn_train, ...
%       ifrac1_mtx_test, ifrac1_mn_test);



%    plot_nonlinearity(x0bins, pspk, pspkx1, 'MID1 Analysis');


%    plot_train_test_information(fraction, ifrac1_mtx_train, ...
%       ifrac1_mn_train, ifrac1_std_train, ifrac1_mtx_test, ifrac1_mn_test, ...
%       ifrac1_std_test, 'MID1 Info Analysis');


   % Save MID1 data
   %---------------------------------------------------------------
   px_mid1{i} = px1;
   pxspk_mid1{i} = px1spk;
   pspkx_mid1{i} = pspkx1;


   %--------------------------------------------------------------------
   %   MID2 Filter from MID Code
   %--------------------------------------------------------------------

   fprintf('\nMID2\n');

   [pspk, px2, px2spk, pspkx2] = proj_prob_xbins_dist(x0bins, x2train, x2train_locator);

%    [ifrac2_train] = train_test_info_fraction(x2bins, x2train_locator, x2train, fraction); % training information
% 
%    [ifrac2_test] = train_test_info_fraction(x2bins, x2test_locator, x2test, fraction); % test information
% 
%    [ifrac2_mn_train, ifrac2_std_train, ifrac2_mtx_train] = ...
%       train_test_info_fraction_mean_std(fraction, ifrac2_train);
% 
%    [ifrac2_mn_test, ifrac2_std_test, ifrac2_mtx_test] = ...
%       train_test_info_fraction_mean_std(fraction, ifrac2_test);
% 
%    [info2_extrap_train, info2_extrap_test] = ...
%       train_test_info_extrapolate(fraction, ifrac2_mtx_train, ifrac2_mn_train, ...
%       ifrac2_mtx_test, ifrac2_mn_test);


   plot_nonlinearity(x0bins, pspk, pspkx2, 'MID2 Analysis');

pause
 
%    plot_train_test_information(fraction, ifrac2_mtx_train, ...
%       ifrac2_mn_train, ifrac2_std_train, ifrac2_mtx_test, ifrac2_mn_test, ...
%       ifrac2_std_test, 'MID2 Info Analysis');

   % Save MID2 data
   %---------------------------------------------------------------
   px_mid2{i} = px2;
   pxspk_mid2{i} = px2spk;
   pspkx_mid2{i} = pspkx2;

end


fio1.location = location;
fio1.unit = unit;
fio1.xbins = xbins;
fio1.pspk = pspk;

fio1.px0 = px_sta; % this and the following are cell arrays of cell arrays
fio1.px0spk = pxspk_sta;
fio1.pspkx0 = pspkx_sta;

fio1.px1 = px_mid1;
fio1.px1spk = pxspk_mid1;
fio1.pspkx1 = pspkx_mid1;

fio1.px2 = px_mid2;
fio1.px2spk = pxspk_mid2;
fio1.pspkx2 = pspkx_mid2;

return;



   %--------------------------------------------------------------------
   %   MID1 and MID2 Information
   %--------------------------------------------------------------------

   fprintf('\nMID1 and MID2\n');

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


return;







