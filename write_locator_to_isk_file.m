function write_locator_to_isk_file(sprfile, spk, trigger)
% write_locator_to_isk_file  Make file w/ spike train
%
% write_locator_to_isk_file(sprfile, spk, trigger)
% -----------------------------------------------------------------
%
% Takes the stimulus file sprfile and the spike train struct array, with
% corresponding trigger variable, and makes a binned spike train. The 
% spike trainis binned at the resolutio of the sprfile. The spike train
% is an integer vector, having values, 0, 1, 2, etc., though usually
% the values are 0 or 1.
%
% The function uses the stimulus and spike times to create a vector that 
% has as many elements as there are time bins in the stimulus. At any time 
% bin, a spike may have occurred - the vector keeps track of this number. 
% The vector is then save to a text file that end in *.isk
%
% The *.isk nomenclature follows the format from Tanya Sharpee's MID code.
%
% write_locator_to_isk_file(sprfile, spk, trigger);
%
% sprfile : full path to the ripple envelope file.
%
% For experiment 20040217:
%
%    C:\MATLABR2007b\work\stimuli\20040217\dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8.spr
%
% If you input [] for sprfile, the spr file for 20040217 will be used.
%
% spk : struct array of spike times for different neurons.
%
% trigger : vector of sample numbers, each of which corresponds to a 
% pulse during experimental stimulation. The triggers allow us to
% synchronize the stimulus and the spike train.
%
% No output variable are returned. A *.isk is created in the current
% directory from which the function is called.
%
% caa 1/1/09


% Old spr file that was used in an IC experiment
% sprfile = 'dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8.spr';


if ( isempty(sprfile) )
   sprfile = 'C:\MATLABR2007b\work\stimuli\20040217\dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min_DFt2_DFf8.spr';
end

% Don't really need the following two variables. The function
% get_locator_for_mid_analysis.m requires them, but we can change this
% since that function doesn't do anything with them.

t1 = 0;
t2 = 0.1;


for i = 1:length(spk) %Ncells

   % Get spike times from ms to sample number
	fsspk = spk(i).fs;
	spet = spk(i).spiketimes / 1000 * fsspk;

   % Get the locator variable, which contains number of spikes in a trial
	[t, f, locator, ns, ar] = get_locator_for_mid_analysis(sprfile, t1, t2, spet, trigger, fsspk);


   % Make the file name for locator
	exp = spk(i).exp;
	site = spk(i).site;
	depth = spk(i).depth;
	atten = spk(i).atten;
	stim = spk(i).stim;


	spkfile = sprintf('%s-site%.0f-%.0fum-%.0fdb-%s-fs%.0f-%u.isk', ...
		exp, site, depth, atten, stim, fsspk, i);

   fprintf('spkfile = %s\n', spkfile);

   % Write the locator to the file
	fid = fopen( spkfile, 'w' );
	count = fprintf( fid, '%u \n', locator );
	fclose( fid );

   fprintf( 'Ntrials = %.0f\n\n', length(locator) );

end

return;


