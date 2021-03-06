function [x0train, x0train_locator, x0train_params, x0spktrain, x0spktrain_params] = dmr_params_from_mid_projections(locator, stimulus, params, spr_paramfile, x0train, x0train_locator, x0train_params)
% dmr_params_from_mid_projections - modulation parameters for each evoked spike
%
% [x0train, x0train_locator, x0train_params] = ...
%    dmr_params_from_mid_projections(locator, stimulus, params, spr_paramfile)
% ----------------------------------------------------
% Address simple question: Given a spike, what were the temporal and spectral
% modulation parameters?
%
% That's all we care about with this function!
%
% locator : vector of integers, where values greater
%           than one imply a spike, and values of 0
%           imply no spike.
%
%           You get the locator from a .isk file and using the call
%           locator = locator_from_isk_file(iskfile);
%
% stimulus : the entire ripple stimulus envelope file
%            as one matrix.
%
% params : struct array holding data specs for the cell to be analyzed.
%
%     params.prefix = where the MID analysis data are located
%     params.location = location, or site number, of cell
%     params.cell = cell number, which cell at the location
%     params.optlevel = optimization level; 1 if MID code was run for first
%                       MID only; 2 if the code was run for both MIDs
%     params.cy = how many time bins to average; should usually be 1
%     params.nt = number of time bins in filter; always 20
%     params.nf = num of freq bins; always 25
%     params.x0 = what was the first stimulus frequency used for analysis.
%                 Usually something like 20.
%
% caa 12/22/06


if ( nargin < 2 | nargin > 7 )

   error('You need 2 or 4 input args.');

elseif ( nargin == 2 | nargin == 7)

   params.prefix = 'C:\MATLABR2007b\work\mid_icc\2003-11-24-site8\dat_files\';
   params.location = 708;
   params.cell = 3;
   params.optlevel = 2;
   params.cy = 1;
   params.nt = 20;
   params.nf = 25;
   params.x0 = 20;

   prefix = params.prefix; %prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';
   location = params.location; %707;
   cell = params.cell; % 1;
   optlevel = params.optlevel; %2;
   cy = params.cy; %1;
   numtbins = params.nt; %20;
   numfbins = params.nf; %25;
   x0 = params.x0; %20;
   index_freq = (x0):(numfbins-1+x0);

%    spr_paramfile = 'C:\MATLABR2007b\work\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat';
   spr_paramfile = 'C:\Users\Craig\stimuli\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat';

elseif ( nargin == 3 )

   prefix = params.prefix; %prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';
   location = params.location; %707;
   cell = params.cell; % 1;
   optlevel = params.optlevel; %2;
   cy = params.cy; %1;
   numtbins = params.nt; %20;
   numfbins = params.nf; %25;
   x0 = params.x0; %20;
   index_freq = (x0):(numfbins-1+x0);

   spr_paramfile = 'C:\MATLABR2007b\work\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat';

else

   prefix = params.prefix; %prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';
   location = params.location; %707;
   cell = params.cell; % 1;
   optlevel = params.optlevel; %2;
   cy = params.cy; %1;
   numtbins = params.nt; %20;
   numfbins = params.nf; %25;
   x0 = params.x0; %20;
   index_freq = (x0):(numfbins-1+x0);

end

load(spr_paramfile, 'faxis', 'taxis', 'MaxRD', 'MaxFM');

if ( length(locator) ~= size(stimulus,2) )
   error('Spike train and envelope file have different number of trials.');
end


fraction = [80 85 90 92.5 95 97.5 100];


%--------------------------------------------------------------------
%   STA Filter from MID Code
%--------------------------------------------------------------------

   file_sta = sprintf('%srpsta_%u_%u_1x%ux%u_1', prefix, location, cell, numfbins, numtbins);

   sta = sta_testrep(locator, stimulus, index_freq, numfbins, numtbins);

   close('all');

if ( nargin < 7 )

%    subplot(2,2,1);
%    imagesc(sta{1});
%    subplot(2,2,2);
%    imagesc(sta{2});
%    subplot(2,2,3);
%    imagesc(sta{3});
%    subplot(2,2,4);
%    imagesc(sta{4});

   [tmf, xmf, rtf] = strf2rtf(taxis, faxis, sta{1}, MaxFM, MaxRD, 'y');
% want to use the function strf2rtf.m here. Or copy and rename it to sta2rtf.m
% [Fm, RD, RTF] = strf2rtf(taxis, faxis, STRF, MaxFm, MaxRD, Display)

   pause

   close all;

   % return;

tic
%    [x0train, x0train_locator, x0train_params] = ...
   [x0train, x0train_locator, x0train_params, x0spktrain, x0spktrain_params] = ...
   train_test_projection_params(sta, locator, stimulus, index_freq, taxis, faxis, MaxFM, MaxRD);
toc

   [x0bins, pspk, px0, px0spk, pspkx0] = proj_prob_dist(x0train, x0train_locator); % for the nonlinearity

else

   [x0bins, pspk, px0, px0spk, pspkx0] = proj_prob_dist_params(x0train, x0train_locator, x0train_params); % for the nonlinearity

end


% return;

% [ifrac0_train] = train_test_info_fraction(x0bins, x0train_locator, x0train, fraction); % training information
% 
% [ifrac0_test] = train_test_info_fraction(x0bins, x0test_locator, x0test, fraction); % test information
% 
% [ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_train] = ...
% train_test_info_fraction_mean_std(fraction, ifrac0_train);
% 
% [ifrac0_mn_test, ifrac0_std_test, ifrac0_mtx_test] = ...
% train_test_info_fraction_mean_std(fraction, ifrac0_test);
% 
% [info0_extrap_train, info0_extrap_test] = ...
% train_test_info_extrapolate(fraction, ifrac0_mtx_train, ifrac0_mn_train, ...
% ifrac0_mtx_test, ifrac0_mn_test);

plot_filters_nonlinearities(sta, x0bins, pspk, pspkx0, 'STA Analysis');

return;




% [v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_filter(file_sta, numtbins, numfbins);
% 
% sta{1} = reshape(mtx_sta(:,1), numfbins, numtbins);
% sta{2} = reshape(mtx_sta(:,2), numfbins, numtbins);
% sta{3} = reshape(mtx_sta(:,3), numfbins, numtbins);
% sta{4} = reshape(mtx_sta(:,4), numfbins, numtbins);
% 
% [x0train, x0test, x0train_locator, x0test_locator] = ...
% train_test_projection(sta, locator, stimulus, index_freq);
% 
% [x0bins, pspk, px0, px0spk, pspkx0] = proj_prob_dist(x0train, x0train_locator); % for the nonlinearity
% 
% [ifrac0_train] = train_test_info_fraction(x0bins, x0train_locator, x0train, fraction); % training information
% 
% [ifrac0_test] = train_test_info_fraction(x0bins, x0test_locator, x0test, fraction); % test information
% 
% [ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_train] = ...
% train_test_info_fraction_mean_std(fraction, ifrac0_train);
% 
% [ifrac0_mn_test, ifrac0_std_test, ifrac0_mtx_test] = ...
% train_test_info_fraction_mean_std(fraction, ifrac0_test);
% 
% [info0_extrap_train, info0_extrap_test] = ...
% train_test_info_extrapolate(fraction, ifrac0_mtx_train, ifrac0_mn_train, ...
% ifrac0_mtx_test, ifrac0_mn_test);
% 
% plot_filters_nonlinearities(sta, x0bins, pspk, pspkx0, 'STA Analysis');
% 
% plot_train_test_information(fraction, ifrac0_mtx_train, ...
% ifrac0_mn_train, ifrac0_std_train, ifrac0_mtx_test, ifrac0_mn_test, ...
% ifrac0_std_test, 'STA Info Analysis');
% 


return;

%--------------------------------------------------------------------
%   MID1 Filter from MID Code
%--------------------------------------------------------------------

% prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';
file_mid1 = sprintf('%srpdtest%u_prelim_v1_%u_%u_1x%ux%u_%u', prefix, optlevel, location, cell, numfbins, numtbins, cy);

[v_mid1, coeff_mid1, projection_mid1, mtx_mid1] = get_auditory_filter(file_mid1, numtbins, numfbins);

mid1{1} = reshape(mtx_mid1(:,1), numfbins, numtbins);
mid1{2} = reshape(mtx_mid1(:,2), numfbins, numtbins);
mid1{3} = reshape(mtx_mid1(:,3), numfbins, numtbins);
mid1{4} = reshape(mtx_mid1(:,4), numfbins, numtbins);

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

plot_filters_nonlinearities(mid1, x1bins, pspk, pspkx1, 'MID1 Analysis');

plot_train_test_information(fraction, ifrac1_mtx_train, ...
ifrac1_mn_train, ifrac1_std_train, ifrac1_mtx_test, ifrac1_mn_test, ...
ifrac1_std_test, 'MID1 Info Analysis');



%--------------------------------------------------------------------
%   MID2 Filter from MID Code
%--------------------------------------------------------------------

% prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';
file_mid2 = sprintf('%srpdtest%u_prelim_v2_%u_%u_1x%ux%u_%u', prefix, optlevel, location, cell, numfbins, numtbins, cy);

[v_mid2, coeff_mid2, projection_mid2, mtx_mid2] = get_auditory_filter(file_mid2, numtbins, numfbins);

mid2{1} = reshape(mtx_mid2(:,1), numfbins, numtbins);
mid2{2} = reshape(mtx_mid2(:,2), numfbins, numtbins);
mid2{3} = reshape(mtx_mid2(:,3), numfbins, numtbins);
mid2{4} = reshape(mtx_mid2(:,4), numfbins, numtbins);

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


plot_filters_nonlinearities(mid2, x2bins, pspk, pspkx2, 'MID2 Analysis');

plot_train_test_information(fraction, ifrac2_mtx_train, ...
ifrac2_mn_train, ifrac2_std_train, ifrac2_mtx_test, ifrac2_mn_test, ...
ifrac2_std_test, 'MID2 Info Analysis');


%--------------------------------------------------------------------
%   MID1 and MID2 Information
%--------------------------------------------------------------------


% [xbins, pspk, px1x2, px1x2spk, pspkx1x2] = proj_prob_dist_2d(x1train_locator, x1train, x2train); % for 2D nonlinearity


[ifrac12_train] = train_test_info_fraction_2d(x1bins, x1train_locator, x1train, x2train, fraction);

[ifrac12_test] = train_test_info_fraction_2d(x1bins, x1test_locator, x1test, x2test, fraction);


[ifrac12_mn_train, ifrac12_std_train, ifrac12_mtx_train] = ...
train_test_info_fraction_mean_std(fraction, ifrac12_train);

[ifrac12_mn_test, ifrac12_std_test, ifrac12_mtx_test] = ...
train_test_info_fraction_mean_std(fraction, ifrac12_test);


[info12_extrap_train, info12_extrap_test] = ...
train_test_info_extrapolate(fraction, ifrac12_mtx_train, ifrac12_mn_train, ...
ifrac12_mtx_test, ifrac12_mn_test);

plot_train_test_information(fraction, ifrac12_mtx_train, ...
ifrac12_mn_train, ifrac12_std_train, ifrac12_mtx_test, ifrac12_mn_test, ...
ifrac12_std_test, 'MID1 and MID2 Info Analysis');


% infoextrap.train0 = info0_extrap_train;
infoextrap.train1 = info1_extrap_train;
infoextrap.train2 = info2_extrap_train;
infoextrap.train12 = info12_extrap_train;

% infoextrap.test0 = info0_extrap_test;
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
fprintf('MID2 Test Info  : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test1);
fprintf('MID12 Test Info : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test12);
fprintf('\n');
fprintf('STA / MID1      : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test0 ./ infoextrap.test1)
fprintf('MID1 / MID12    : %.4f  %.4f  %.4f  %.4f  %.4f\n\n', infoextrap.test1 ./ infoextrap.test12)
fprintf('\n');

return;









