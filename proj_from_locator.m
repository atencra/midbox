function proj_from_locator(locator, stimulus)
% get_sta_from_locator - calculate STA from spike train vector and
%    a stimulus matrix
%
% sta = proj_from_locator(locator, stimulus)
% ----------------------------------------------------
%
% locator : vector of integers, where values greater
%           than one imply a spike, and values of 0
%           imply no spike
%
% stimulus : the entire ripple stimulus envelope file
%            as one matrix.
%
% sta : spike triggered average. Will always have dimensions 25x20.
% 
% caa 12/22/06


if ( length(locator) ~= size(stimulus,2) )
   error('Spike train and envelope file have different number of trials.');
end


% % From Tanya's MID code:
% % ----------------------------------------------------
% if ( signl == 701) {     
% dimx = 61;
% dimy = 1;
% Nh = 25;     
% x0 = 20; // dimx-Nh+1;     
% cy = 1;
% Movie_length = 1115226;
% }
% 
% 
% for(i=1;i<=Ntrials*Nh;i++)  stimuli[i]=0;     
% 
% for(i=1;i<=Ntrials;i++){       
% 	for(k=1;k<=cy;k++){ 	
% 		fread(bfloat,sizeof(float),dimx*dimy,inf); 	
% 		if (ferror(inf)) myerror("error reading from file dmr-500flo-40000fhi-4SM-40TM-    		40db-96khz-48DF-21min_6_carriers_per_octave.spr"); 	
% 
% 		for(j=1;j<=Nh;j++){ 	  
% 		// for x0 such that last Nh points are taken, x0=dimx-Nh+1 	  
% 		// stimuli[(i-1)*Nh+j] = (double)bfloat[j-1+(dimx-Nh)]; previous version 	  
% 		// stimuli[(i-1)*Nh+j] = (double)bfloat[j+x0-1];  	  
% 		// x0 = dimx - Nh + 1;
% 
% 		stimuli[(i-1)*Nh+j] +=  (double) bfloat[j-1+x0-1];  	
% 		} 
% 	}
% }







numtbins = 20;
numfbins = 25;
x0 = 20;
% Index to get correct frequency bins. I think the code says use:
% (x0-1):(numfbins-1+x0-1), but (x0):(numfbins-1+x0) gives the correct 
% result.
index_freq = (x0):(numfbins-1+x0);


%--------------------------------------------------------------------
%   STA Filter from MID Code
%--------------------------------------------------------------------
location = 707;
cell = 1;
Nh = 20;
Nv = 25;
nlags = 1;
Nparts = 4;

prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';
file_sta = sprintf('%srpsta_%u_%u_1x%ux%u_1', prefix, location, cell, numfbins, numtbins);
[v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_filter(file_sta, Nh, Nv, nlags, Nparts);

sta = sta_testrep(locator, stimulus, index_freq, numfbins, numtbins);

[sta_mean, coeff, projection, mtx] = mean_sta_filter(sta);


figure;
subplot(1,2,1);
imagesc( sta_mean );
minmin = min(min(sta_mean));
maxmax = max(max(sta_mean));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
h = colorbar;
cmap = color_brewer_colormap('rdbu');
colormap(jet);

subplot(1,2,2);
imagesc( reshape(v_sta, numfbins, numtbins) );
minmin = min(min(v_sta));
maxmax = max(max(v_sta));
boundary = max([abs(minmin) abs(maxmax)]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
h = colorbar;
cmap = color_brewer_colormap('rdbu');
colormap(jet);

return;


% Prior distribution of projection values for training set spikes 
% onto each training set STA

tic
[xtrain, xtest, xtrain_locator, xtest_locator] = train_test_projection(sta, locator, stimulus, index_freq);
toc

tic
[xbins, pspk, px, pxspk, pspkx] = proj_prob_dist(xtrain, xtrain_locator); % for the nonlinearity
toc

tic
[train_fraction, info_train_fraction] = train_test_info_fraction(xbins, xtrain_locator, xtrain); % training information
toc

tic
[test_fraction, info_test_fraction] = train_test_info_fraction(xbins, xtest_locator, xtest); % test information
toc

tic
[info_frac_test_mn, info_frac_test_std, info_frac_test_mtx] = train_test_info_fraction_mean_std(test_fraction, info_test_fraction);
toc

figure;

subplot(2,4,1);
imagesc(sta{1});
axis xy;
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
colorbar;

subplot(2,4,2);
imagesc(sta{2});
axis xy;
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
colorbar;

subplot(2,4,3);
imagesc(sta{3});
axis xy;
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
colorbar;

subplot(2,4,4);
imagesc(sta{4});
axis xy;
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
colorbar;

maxmax = max([max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

subplot(2,4,5);
hold on;
plot(xbins, pspkx{1}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{1} pspk{1}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

subplot(2,4,6);
hold on;
plot(xbins, pspkx{2}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{2} pspk{2}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

subplot(2,4,7);
hold on;
plot(xbins, pspkx{3}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{3} pspk{3}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

subplot(2,4,8);
hold on;
plot(xbins, pspkx{4}, 'ko-', 'markerfacecolor', 'k');
plot([-7 7], [pspk{4} pspk{4}], 'k--');
xlim([-8 8]);
ylim([0 maxmax]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);


figure; 

subplot(2,2,1);
plot( 1./test_fraction, info_frac_test_mtx(1,:), 'ko');
ylim([0 1.1*max( info_frac_test_mtx(1,:) )]);

subplot(2,2,2);
plot( 1./test_fraction, info_frac_test_mtx(2,:), 'ko');
ylim([0 1.1*max( info_frac_test_mtx(2,:) )]);

subplot(2,2,3);
plot( 1./test_fraction, info_frac_test_mtx(3,:), 'ko');
ylim([0 1.1*max( info_frac_test_mtx(3,:) )]);

subplot(2,2,4);
plot( 1./test_fraction, info_frac_test_mtx(4,:), 'ko');
ylim([0 1.1*max( info_frac_test_mtx(4,:) )]);


figure; 

errorbar( 1./test_fraction, info_frac_test_mn, info_frac_test_std, 'ko');
ylim( [0 1.1*max(info_frac_test_mn)+ max(info_frac_test_std)] );


return;



% Prior distribution of projection values for test set spikes the 
% onto training set STA


% First we get the prior distribution of projection values
xprior = zeros(size(locator));
[nr, nc] = size(sta); % # frequencies, # time bins

for i = nc:length(locator)

   xprior(i) =  sum( sum ( stimulus( index_freq, i-nc+1:i ) .* sta ) ); % inner product

   if ( ~mod(i, 100000) )
      fprintf('i = %.0f\n', i);
   end
end % (for i)




% Probability of a spike in a given trial
nspikes = sum(locator);
pspk = nspikes / length(locator); % probability of a spike

% Normalize projection values to mean, sd of prior
mean_prior = mean(xprior);
std_prior = std(xprior);

x = (xprior - mean_prior) ./ std_prior;
xspk = x( locator > 0 ); % values corresponding to a spike



% Probability distributions of prior, posterior projection values

xbins = -7:7;
xbins = xbins(:);

nx = hist(x, xbins);
px = nx ./ sum(nx); % p(x)
px = px(:);

nxspk = hist(xspk, xbins);
pxspk = nxspk ./ sum( nxspk ); % p(x|spk)
pxspk = pxspk(:);

pspkx = pspk .* pxspk ./ px;

% [xbins px pxspk pspkx]

iplugin = info_px_pxspk(px, pxspk);

frac = 75;
% [imn, isd] = info_fraction(locator, x, xbins, frac);

[i75] = info_fraction(locator, x, xbins, 75);
[i80] = info_fraction(locator, x, xbins, 80);
[i85] = info_fraction(locator, x, xbins, 85);
[i90] = info_fraction(locator, x, xbins, 90);
[i95] = info_fraction(locator, x, xbins, 95);
[i100] = info_fraction(locator, x, xbins, 100);

i100


close all;

% figure;
% imagesc(sta);
% colorbar;

figure;

subplot(2,2,1);
hist(x, 500);
xlim([-10 10]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

subplot(2,2,2);
hold on;
plot(xbins, px, 'k-');
plot([-7 7], [pspk pspk], 'k--');
xlim([-7 7]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);


subplot(2,2,3);
hist(xspk, 500);
xlim([-10 10]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);


subplot(2,2,4);
hold on;
plot(xbins, px, 'k-');
plot(xbins, pxspk, 'r-');
plot(xbins, pspkx, 'g-');
plot([-7 7], [pspk pspk], 'k--');
legend('P(x)', 'P(x|spk)', 'P(spk|x)', 'P(spk)');
xlim([-7 7]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

set(gcf, 'position', [172 114 968 772]);

return;


function iplugin = info_px_pxspk(px, pxspk)

index = find( px>0 & pxspk>0 );

iplugin = sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) );

return;


function [ifrac] = info_fraction(locator, x, xbins, frac)

frac = frac / 100;

nreps = 5;

ntrials = length(locator);

nfrac = round(frac * ntrials); % number of reduced trials

ifrac = zeros(1,nreps);

rand('state', 0);

for i = 1:nreps

   first = ceil( rand(1) * (ntrials-nfrac+1) );
   last = first + nfrac - 1;

   index_trials = first:last;

   locator_temp = locator( index_trials );
   xtemp = x( index_trials );
   xspktemp = xtemp( locator_temp > 0);


   nx = hist(xtemp, xbins);
   px = nx ./ sum(nx); % p(x)
   px = px(:);

   nxspk = hist(xspktemp, xbins);
   pxspk = nxspk ./ sum( nxspk ); % p(x|spk)
   pxspk = pxspk(:);

   iplugin = info_px_pxspk(px, pxspk);

   ifrac(i) = iplugin;

end % (for i)

return;


function [info_frac_mn, info_frac_std, info_mtx] = train_test_info_fraction_mean_std(fraction, info_fraction)
%
% xbins : vector of bins for probability histograms
%
% locator : Tells if there was a spike at a given trial. Cell array 
% with 4 elements. Each element is a vector.
%
% x : Contains the projection values. Cell array with 4 elements. Each 
% element is a vector.
%
% caa 3/13/09



% We have four data sets. length(x) = 4, with x a cell array.

% For one of the data sets, we need to get the information fractions
% for 75, 80, 85, 90, 95, and 100 percent of the data.

% So we need two loops: one for data set number, and one for data fraction

% As output we'll use a cell array that is 4xlength(fraction),
% where

% fraction = [75 80 85 90 95 100];


if ( length(fraction) ~= size(info_fraction,2) )
   error('fraction and size(info_fraction,2) do not match.');
end


info_mtx = zeros( size(info_fraction,1), size(info_fraction,2) );

for i = 1:size(info_fraction,1) % for each data set

   for j = 1:size(info_fraction,2) % for each fraction of the data

      temp = info_fraction{i}{j};

      info_mtx(i, j) = mean(temp);

   end % (for j)

end % (for i)

info_frac_mn = mean( info_mtx, 1 );
info_frac_std = std( info_mtx, 0, 1 );

return;



function [fraction, ifraction] = train_test_info_fraction(xbins, locator, x)
%
% xbins : vector of bins for probability histograms
%
% locator : Tells if there was a spike at a given trial. Cell array 
% with 4 elements. Each element is a vector.
%
% x : Contains the projection values. Cell array with 4 elements. Each 
% element is a vector.
%
% caa 3/13/09



% We have four data sets. length(x) = 4, with x a cell array.

% For one of the data sets, we need to get the information fractions
% for 75, 80, 85, 90, 95, and 100 percent of the data.

% So we need two loops: one for data set number, and one for data fraction

% As output we'll use a cell array that is 4xlength(fraction),
% where

fraction = [75 80 85 90 92.5 95 97.5 100];

ifraction = cell( length(x), length(fraction) );

for i = 1:length(x)

   for j = 1:length(fraction)

      [ifrac] = info_fraction( locator{i}, x{i}, xbins, fraction(j) );

      ifraction{i}{j} = ifrac;

   end % (for j)

end % (for i)

return;



function sta = sta_testrep(locator, stimulus, index_freq, numfbins, numtbins)

ntrials = length(locator);

[irep1, irep2, irep3, irep4] = training_set_indices(ntrials);


% STA for testrep #1
% ---------------------------------------------------

sta1 = zeros(numfbins, numtbins);

for i = 1:length(irep1)

   index = irep1(i);

   if ( locator( index ) & index >= 20 )

      index_stim = index-numtbins+1:index;

      sta1 = sta1 + locator( index ) * stimulus(index_freq, index_stim);

   end

end


% STA for testrep #2
% ---------------------------------------------------

sta2 = zeros(numfbins, numtbins);

for i = 1:length(irep2)

   index = irep2(i);

   if ( locator( index ) & index >= 20 )

      index_stim = index-numtbins+1:index;

      sta2 = sta2 + locator( index ) * stimulus(index_freq, index_stim);

   end

end



% STA for testrep #3
% ---------------------------------------------------

sta3 = zeros(numfbins, numtbins);

for i = 1:length(irep3)

   index = irep3(i);

   if ( locator( index ) & index >= 20 )

      index_stim = index-numtbins+1:index;

      sta3 = sta3 + locator( index ) * stimulus(index_freq, index_stim);

   end

end



% STA for testrep #4
% ---------------------------------------------------

sta4 = zeros(numfbins, numtbins);

for i = 1:length(irep4)

   index = irep4(i);

   if ( locator( index ) & index >= 20 )

      index_stim = index-numtbins+1:index;

      sta4 = sta4 + locator( index ) * stimulus(index_freq, index_stim);

   end

end


sta{1} = sta1;
sta{2} = sta2;
sta{3} = sta3;
sta{4} = sta4;


return;


function [irep1, irep2, irep3, irep4] = training_set_indices(ntrials)
% testrep_indices - indices for test repetition 
%
% ntrials : total number of trials in data set. Usually a number like
% 1115226, or something similar and of that magnitude.
%
% We separate data into training and test sets. Training sets are made
% up from 3/4 of the data. We can divide the data four ways:
%
% training set 1 = [2 3 4]
% training set 2 = [1 3 4]
% training set 3 = [1 2 4]
% training set 4 = [1 2 3]
%
% where 1, 2, 3, 4 represent the first 1/4, the second 1/4, the third 1/4,
% and the fourth 1/4 of the data set.
%
% The test set then is comprised of 1/4 of the data. We can notate this as
%
% test set 1 = [1]
% test set 2 = [2]
% test set 3 = [3]
% test set 4 = [4]
%
% So this function returns the indices into the 4 different training sets.
%
% caa 3/12/09

ntrials_testrep = floor( ntrials / 4 );

index1 = [1:ntrials_testrep];
index2 = [ntrials_testrep+1:2*ntrials_testrep];
index3 = [2*ntrials_testrep+1:3*ntrials_testrep];
index4 = [3*ntrials_testrep+1:ntrials];

irep1 = [index2 index3 index4];
irep2 = [index1 index3 index4];
irep3 = [index1 index2 index4];
irep4 = [index1 index2 index3];

return;



function [irep1, irep2, irep3, irep4] = test_set_indices(ntrials)
% testrep_indices - indices for test repetition 
%
% ntrials : total number of trials in data set. Usually a number like
% 1115226, or something similar.
%
% We separate data into training and test sets. Training sets are made
% up from 3/4 of the data. We can divide the data four ways:
%
% training set 1 = [2 3 4]
% training set 2 = [1 3 4]
% training set 3 = [1 2 4]
% training set 4 = [1 2 3]
%
% where 1, 2, 3, 4 represent the first 1/4, the second 1/4, the third 1/4,
% and the fourth 1/4 of the data set.
%
% The test set then is comprised of 1/4 of the data. We can notate this as
%
% test set 1 = [1]
% test set 2 = [2]
% test set 3 = [3]
% test set 4 = [4]
%
% So this function returns the indices into the 4 different test sets.
%
% caa 3/12/09
 

ntrials_testrep = floor( ntrials / 4 );

irep1 = [1:ntrials_testrep];
irep2 = [ntrials_testrep+1:2*ntrials_testrep];
irep3 = [2*ntrials_testrep+1:3*ntrials_testrep];
irep4 = [3*ntrials_testrep+1:ntrials];

return;



function [xtrain, xtest, xtrain_locator, xtest_locator] = train_test_projection(filters, locator, stimulus, index_freq)

ntrials = length(locator);

[itrain1, itrain2, itrain3, itrain4] = training_set_indices(ntrials);
[itest1, itest2, itest3, itest4] = test_set_indices(ntrials);

% Prior distribution of projection values for test set spikes the 
% onto training set STA


% Project all trials onto filter to get the prior distribution.
% This distribution will include the training set projections and the
% test set projections. All trials will have a projection value
% associated with them.


% Find the prior projection values for all filters at the same time:
x1prior = zeros(ntrials, 1);
x2prior = zeros(ntrials, 1);
x3prior = zeros(ntrials, 1);
x4prior = zeros(ntrials, 1);

f1 = filters{1};
f2 = filters{2};
f3 = filters{3};
f4 = filters{4};

[nr, nc] = size(f1); % # frequencies, # time bins

for i = nc:ntrials

   x1prior(i) =  sum( sum ( stimulus( index_freq, i-nc+1:i ) .* f1 ) ); % inner product
   x2prior(i) =  sum( sum ( stimulus( index_freq, i-nc+1:i ) .* f2 ) ); % inner product
   x3prior(i) =  sum( sum ( stimulus( index_freq, i-nc+1:i ) .* f3 ) ); % inner product
   x4prior(i) =  sum( sum ( stimulus( index_freq, i-nc+1:i ) .* f4 ) ); % inner product

   if ( ~mod(i, 100000) )
      fprintf('i = %.0f\n', i);
   end
end % (for i)


% For the training set:
x1train = x1prior(itrain1); % include prior and posterior
x1train_locator = locator(itrain1);
x1spktrain = x1train( x1train_locator > 0 ); % includes only posterior

x2train = x2prior(itrain2); % include prior and posterior
x2train_locator = locator(itrain2);
x2spktrain = x2train( x2train_locator > 0 ); % includes only posterior

x3train = x3prior(itrain3); % include prior and posterior
x3train_locator = locator(itrain3);
x3spktrain = x3train( x3train_locator > 0 ); % includes only posterior

x4train = x4prior(itrain4); % include prior and posterior
x4train_locator = locator(itrain4);
x4spktrain = x4train( x4train_locator > 0 ); % includes only posterior


% For the test set:
x1test = x1prior(itest1); % includes prior and posterior
x1test_locator = locator(itest1);
x1spktest = x1test( x1test_locator > 0 ); % includes only posterior

x2test = x2prior(itest2); % includes prior and posterior
x2test_locator = locator(itest2);
x2spktest = x2test( x2test_locator > 0 ); % includes only posterior

x3test = x3prior(itest3); % includes prior and posterior
x3test_locator = locator(itest3);
x3spktest = x3test( x3test_locator > 0 ); % includes only posterior

x4test = x4prior(itest4); % includes prior and posterior
x4test_locator = locator(itest4);
x4spktest = x4test( x4test_locator > 0 ); % includes only posterior



% Save training data to output cell arrays
xtrain{1} = x1train;
xtrain{2} = x2train;
xtrain{3} = x3train;
xtrain{4} = x4train;

xtrain_locator{1} = x1train_locator;
xtrain_locator{2} = x2train_locator;
xtrain_locator{3} = x3train_locator;
xtrain_locator{4} = x4train_locator;

xspktrain{1} = x1spktrain;
xspktrain{2} = x2spktrain;
xspktrain{3} = x3spktrain;
xspktrain{4} = x4spktrain;



% Save test data to output cell arrays
xtest{1} = x1test;
xtest{2} = x2test;
xtest{3} = x3test;
xtest{4} = x4test;

xtest_locator{1} = x1test_locator;
xtest_locator{2} = x2test_locator;
xtest_locator{3} = x3test_locator;
xtest_locator{4} = x4test_locator;

xspktest{1} = x1spktest;
xspktest{2} = x2spktest;
xspktest{3} = x3spktest;
xspktest{4} = x4spktest;


% [length(xtrain) length(xtest) length(xtrain_locator) length(xtest_locator)]


return;




function [xbins, pspk, px, pxspk, pspkx] = proj_prob_dist(xtrain, xtrain_locator)


% Get probability distributions for training data sets

% We need to get: p(spk), p(x), p(x|spk), and p(spk|x)

xbins = -7:7;
xbins = xbins(:);

for i = 1:length(xtrain)

   locator = xtrain_locator{i};
   xprior = xtrain{i};

   % Probability of a spike in a given trial
   nspikes = sum(locator);
   pspk_temp = nspikes / length(locator); % probability of a spike

   % Normalize projection values to mean, sd of prior
   mean_prior = mean(xprior);
   std_prior = std(xprior);

   x = (xprior - mean_prior) ./ std_prior;
   xspk = x( locator > 0 ); % values corresponding to a spike


   % Probability distributions of prior, posterior projection values

   nx = hist(x, xbins);
   px_temp = nx ./ sum(nx); % p(x)
   px_temp = px_temp(:);

   nxspk = hist(xspk, xbins);
   pxspk_temp = nxspk ./ sum( nxspk ); % p(x|spk)
   pxspk_temp = pxspk_temp(:);

   pspkx_temp = pspk_temp .* pxspk_temp ./ px_temp; % p(spk|x)


   % Assign output data

   pspk{i} = pspk_temp;
   px{i} = px_temp;
   pxspk{i} = pxspk_temp;
   pspkx{i} = pspkx_temp;

end % (for)

return;




function [v_mean, coeff, projection, mtx] = mean_sta_filter(sta)
%
%[v_mean, coeff, projection, mtx] = mean_sta_filter(sta)
%
%
% caa 3/13/09



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
    coeff(i)=sign(sum(mtx(:,i).*mtx(:,1)));
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





