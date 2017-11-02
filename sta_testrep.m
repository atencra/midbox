function sta = sta_testrep(locator, stimulus, index_freq, numfbins, numtbins)
% sta_testrep - Calculate STA for training set data
% 
% sta = sta_testrep(locator, stimulus, index_freq, numfbins, numtbins)
% --------------------------------------------------------------------
% locator : vector, with each element >= 0, which tells if a spike occurred
% at a given trial
% 
% stimulus : matrix, rows represent frequencies, columns represent time, or
% trial. Number of columns in stimulus equals the length of locator.
% 
% index_freq : index to frequencies that were used in the MID analysis. The
% MID analysis used 25 frequencies, so these have to be selected from the
% stimulus
% 
% numfbins : Number of frequency bins used to make the STAs. Always 25.
% 
% numtbins : Number of time bins used to make the STAs. Always 20.
% 
% sta : 1x4 cell array. Each element is the STA estimated from one of the 
% four training data sets.
% 
% caa 3/1/09

fprintf('\nRunning sta_testrep ...\n');


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







