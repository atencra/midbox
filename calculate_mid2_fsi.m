function calculate_mid2_fsi(mid, spk, trigger, sprfile)
% calculate_filter_info - compute projection values for STA, MID1, and
%    MID2. 
%
% calculate_filter_info(mid, spk, trigger, sprfile)
%
% For each filter, either the STA or the MIDs, the stimuli that
% elicit a spike are projected onto each filter, and the projection
% values are saved. Also, a spike train of 50,000 random spikes is
% is created and the stimuli associated with the spikes are also
% projected onto each filter.
%
% The projections values may be later used to calculate information values
% for each filter, or both MID filters.
%
% mid holds the filters, while spk and trigger contain the spike trains
% and stimulus trigger times. The struct array is saved in files of the
% form nonlinearity_params_*.mat
%
% proj is a struct array holding all the projection data.
%
% The sprfile for the data is either  
%
% sprfile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min_DFt22_DFf7.spr';
%
% or
%
% sprfile =
% 'dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min_DFt10_DFf8.spr'
%
% caa 12/18/06

if ( length(mid) ~= length(spk) )
   error('mid and spk need to match.');
end

if ( ~strcmp(mid(1).exp,spk(1).exp) | mid(1).site~=spk(1).site | ...
     mid(1).chan~=spk(1).chan | mid(1).model~=spk(1).model | ...
     ~strcmp(mid(end).exp,spk(end).exp) | mid(end).site~=spk(end).site | ...
     mid(end).chan~=spk(end).chan | mid(end).model~=spk(end).model )
   error('The elements of mid and spk do not match.');
end


T1 = 0.005; % 5 ms 
T2 = 0.095; % 95 ms

i1 = findstr(sprfile,'-');
i2 = findstr(sprfile,'SM');
i3 = findstr(sprfile,'TM');
i4 = findstr(sprfile,'db');

sm = str2num(sprfile(i2(end)-1));
tm = str2num(sprfile(i3(end)-2:i3(end)-1));
mdb = str2num(sprfile(i4(end)-2:i4(end)-1));

atten = spk(1).atten;
fs = spk(1).fs;
spl = 105 - atten;

stype = spk(1).stim;

if ( ~isempty(findstr(stype,'dmr')) )
   snd = 'MR';
elseif ( ~isempty(findstr(stype,'rn')) )
   snd = 'RN';
else
   error('stim in struct array "spk" must be "dmr" or "rn".');
end

modtype = 'dB';
nblocks = 500;  % update display every 50 blocks


mid2fsi = struct(...
'exp',        [], ...
'site',       [], ...
'chan',       [], ...
'model',      [], ...
'depth',      [], ...
'position',   [], ...
'stim',       [], ...
'atten',      [], ...
'spl',        [], ...
'sm',         [], ...
'tm',         [], ...
'mdb',        [], ...
'spetrand',   [], ...
'ccrand',   [], ...
'ccspk',    [], ...
'fsi_cdf',   [], ...
'fsi_mean',    [], ...
'fsi_med', []);


% For experiments up to 2003-4-8
tbins = mid(1).tbins;
fbins = mid(1).fbins;


% create random spike train with 50,000 spikes
mintrig = trigger(3);
maxtrig = trigger(end-3);
spetrand = ceil( (maxtrig-mintrig) * rand(1,50000) + mintrig );
[spetrand, temp] = sort(spetrand);

for i = 1:length(mid)

   sp = spk(i).spiketimes;
   spet = round(sp / 1000 * fs); % convert ms to sample number

   mid2vec = mid(i).rpdtest2_v2.filter; % use the mean mid2
   mid2 = reshape(mid2vec, fbins, tbins);

   % First we compute the prior distribution of correlation coefficient
	% values. For this we use the random spike train. 
   [ccrand, spindexrand] = get_filter_projection(mid2, sprfile, 0.005, 0.095, spetrand, trigger, fs, spl, mdb);


   % Now we want to find the corr coefficient value given a spike, 
	% so use the original spike train and the original filter
   [ccspk, spindex] = get_filter_projection(mid2, sprfile, 0.005, 0.095, spet, trigger, fs, spl, mdb);


	% Now compute the feature selectivity indices

	[fsi_cdf, fsi_mean, fsi_med] = fsi(ccspk, ccrand);

	% Now we're done and we can save the data

   % Save the data in the struct array proj
   mid2fsi(i).exp = mid(i).exp;
   mid2fsi(i).site = mid(i).site;
   mid2fsi(i).chan = mid(i).chan;
   mid2fsi(i).model = mid(i).model;
   mid2fsi(i).depth = mid(i).depth;
   mid2fsi(i).position = mid(i).position;
   mid2fsi(i).stim = mid(i).stim;
   mid2fsi(i).atten = mid(i).atten;
   mid2fsi(i).spl = mid(i).spl;
   mid2fsi(i).sm = mid(i).sm;
   mid2fsi(i).tm = mid(i).tm;
   mid2fsi(i).mdb = mid(i).mdb;
   mid2fsi(i).spetrand = spetrand(:);
   mid2fsi(i).ccrand = ccrand;
   mid2fsi(i).ccspk = ccspk;
   mid2fsi(i).fsi_cdf = fsi_cdf;
   mid2fsi(i).fsi_mean = fsi_mean;
   mid2fsi(i).fsi_med = fsi_med;

   fprintf('%s - site%.0f - chan%.0f - model%.0f  -  %.0f of %.0f completed.\n', ...
      spk(i).exp, spk(i).site, spk(i).chan, spk(i).model, i, length(spk));
   pause(0.5);

end % (for)


return;


