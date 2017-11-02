function info = calculate_filter_info(proj)
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


proj = struct(...
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
'xstarand',   [], ...
'xstaspk',    [], ...
'x1x2rand',   [], ...
'x1x2spk',    [], ...
'x1xrandspk', []);


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


   %=======================================================================
   % For the STA:
   %=======================================================================
   stavec = mid(i).rpsta.filter; % use the mean sta
   sta = reshape(stavec,fbins,tbins);

   % First we compute the prior distribution of projection values. 
   % For this we use the random spike train. 
   [xstarand, spindexrand] = get_filter_projection(sta, sprfile, 0.005, 0.095, spetrand, trigger, fs, spl, mdb);


   % Now we want to find the p(proj|spike), so use the original spike
   % train and the original filter
   [xstaspk, spindex] = get_filter_projection(sta, sprfile, 0.005, 0.095, spet, trigger, fs, spl, mdb);


   %=======================================================================
   % For the first and second MID:
   %=======================================================================
   mid1vec = mid(i).rpdtest2_v1.filter; % use the mean mid1
   mid1 = reshape(mid1vec, fbins, tbins);

   mid2vec = mid(i).rpdtest2_v2.filter; % use the mean mid2
   mid2 = reshape(mid2vec, fbins, tbins);


   % Next we want p(projection onto mid1 & mid2 | random spikes), so use the original filter but this time
   % with 50,000 random spikes 
   [x1x2rand, spindexrand] = get_mid1_mid2_projection(mid1, mid2, sprfile, 0.005, 0.095, spetrand, trigger, fs, spl, mdb);


   % First we want to find the p(proj onto mid1 & mid2 | spike), so use the original spike
   % train and the original filter
   [x1x2spk, spindex] = get_mid1_mid2_projection(mid1, mid2, sprfile, 0.005, 0.095, spet, trigger, fs, spl, mdb);


   % Finally, as a control, we compute projections of the stimulus onto the
   % first MID and onto a random matrix. For the random matrix we will use
   % a randomized MID2 so that the power in the matrix stays the same.
   r = rand(1,length(mid2vec));
   [temp, index] = sort(r);
   mid2vec_rand = mid2vec(index);
   mid2_rand = reshape(mid2vec_rand, fbins, tbins);
   [x1xrandspk, spindex] = get_mid1_mid2_projection(mid1, mid2_rand, sprfile, 0.005, 0.095, spet, trigger, fs, spl, mdb);


   % Save the data in the struct array proj
   proj(i).exp = mid(i).exp;
   proj(i).site = mid(i).site;
   proj(i).chan = mid(i).chan;
   proj(i).model = mid(i).model;
   proj(i).depth = mid(i).depth;
   proj(i).position = mid(i).position;
   proj(i).stim = mid(i).stim;
   proj(i).atten = mid(i).atten;
   proj(i).spl = mid(i).spl;
   proj(i).sm = mid(i).sm;
   proj(i).tm = mid(i).tm;
   proj(i).mdb = mid(i).mdb;
   proj(i).spetrand = spetrand(:);
   proj(i).xstarand = xstarand;
   proj(i).xstaspk = xstaspk;
   proj(i).x1x2rand = x1x2rand;
   proj(i).x1x2spk = x1x2spk;
   proj(i).x1xrandspk = x1xrandspk;

   fprintf('%s - site%.0f - chan%.0f - model%.0f  -  %.0f of %.0f completed.\n', ...
      spk(i).exp, spk(i).site, spk(i).chan, spk(i).model, i, length(spk));
   pause(0.5);

end % (for)





%    % Now we have raw projection values. Next we want to normalize the
%    % values with respect to the standard deviation of p(proj|no spike)
%    mn = mean(projrand);
%    sd = std(projrand);
%    projspk_scaled = (projspk - mn) ./ sd;
%    projrand_scaled = (projrand - mn) ./ sd;
% 
% 
%    % Finally, we calculate probability distributions for different
%    % bin sizes
%    [binsize, info] = get_1d_filter_info(projspk_scaled, projrand_scaled)
% 
% 
%    subplot(3,2,1);
%    hist(projspk,25);
% 
%    subplot(3,2,3);
%    hist(projrand,25);
% 
%    subplot(3,2,5);
%    [nspk, xspk] = hist(projspk_scaled, 25);
%    [nr, xr] = hist(projrand_scaled, 25);
%    plot(xspk, nspk./sum(nspk), 'r-', xr, nr./sum(nr), 'k-');
%    legend('spk', 'no spk');
% 
%    subplot(3,2,2);
%    plot(binsize, info, 'ko-');
%    title('versus bin size');
% 
%    subplot(3,2,4);
%    plot(1 ./ binsize, info, 'ko-');
%    title('versus 1 / bin size');
% 
% 
%    % Now we have raw projection values. Next we want to normalize the
%    % values with respect to the standard deviation of p(proj|no spike)
%    mn = mean(projrand);
%    sd = std(projrand);
%    projspk_scaled = (projspk - mn) ./ sd;
%    projrand_scaled = (projrand - mn) ./ sd;
% 
% 
%    % Finally, we calculate probability distributions for different
%    % bin sizes
%    [binsize, info] = get_2d_filter_info(projspk_scaled, projrand_scaled)







function [binsize, info] = get_1d_filter_info(projspk_scaled, projrand_scaled)

   % Finally, we calculate probability distributions for different
   % bin sizes

   bins9 = linspace(-7,7,9);
   bins11 = linspace(-7,7,11);
   bins13 = linspace(-7,7,13);
   bins15 = linspace(-7,7,15);
   bins17 = linspace(-7,7,17);
   bins19 = linspace(-7,7,19);
   bins21 = linspace(-7,7,21);
   bins25 = linspace(-7,7,25);

   [nspk9] = hist(projspk_scaled, bins9);
   [nrand9] = hist(projrand_scaled, bins9);
   pxspk9 = nspk9 ./ sum(nspk9);
   px9 = nrand9 ./ sum(nrand9);
   index = pxspk9>0 & px9>0;
   info9 = sum( pxspk9(index) .* log2( pxspk9(index)./px9(index) ) );


   [nspk11] = hist(projspk_scaled, bins11);
   [nrand11] = hist(projrand_scaled, bins11);
   pxspk11 = nspk11 ./ sum(nspk11);
   px11 = nrand11 ./ sum(nrand11);
   index = pxspk11>0 & px11>0;
   info11 = sum( pxspk11(index) .* log2( pxspk11(index)./px11(index) ) );

   [nspk13] = hist(projspk_scaled, bins13);
   [nrand13] = hist(projrand_scaled, bins13);
   pxspk13 = nspk13 ./ sum(nspk13);
   px13 = nrand13 ./ sum(nrand13);
   index = pxspk13>0 & px13>0;
   info13 = sum( pxspk13(index) .* log2( pxspk13(index)./px13(index) ) );

   [nspk15] = hist(projspk_scaled, bins15);
   [nrand15] = hist(projrand_scaled, bins15);
   pxspk15 = nspk15 ./ sum(nspk15);
   px15 = nrand15 ./ sum(nrand15);
   index = pxspk15>0 & px15>0;
   info15 = sum( pxspk15(index) .* log2( pxspk15(index)./px15(index) ) );

   [nspk17] = hist(projspk_scaled, bins17);
   [nrand17] = hist(projrand_scaled, bins17);
   pxspk17 = nspk17 ./ sum(nspk17);
   px17 = nrand17 ./ sum(nrand17);
   index = pxspk17>0 & px17>0;
   info17 = sum( pxspk17(index) .* log2( pxspk17(index)./px17(index) ) );

   [nspk19] = hist(projspk_scaled, bins19);
   [nrand19] = hist(projrand_scaled, bins19);
   pxspk19 = nspk19 ./ sum(nspk19);
   px19 = nrand19 ./ sum(nrand19);
   index = pxspk19>0 & px19>0;
   info19 = sum( pxspk19(index) .* log2( pxspk19(index)./px19(index) ) );

   [nspk21] = hist(projspk_scaled, bins21);
   [nrand21] = hist(projrand_scaled, bins21);
   pxspk21 = nspk21 ./ sum(nspk21);
   px21 = nrand21 ./ sum(nrand21);
   index = pxspk21>0 & px21>0;
   info21 = sum( pxspk21(index) .* log2( pxspk21(index)./px21(index) ) );


   [nspk25] = hist(projspk_scaled, bins25);
   [nrand25] = hist(projrand_scaled, bins25);
   pxspk25 = nspk25 ./ sum(nspk25);
   px25 = nrand25 ./ sum(nrand25);
   index = pxspk25>0 & px25>0;
   info25 = sum( pxspk25(index) .* log2( pxspk25(index)./px25(index) ) );

   d9 = min(diff(bins9));
   d11 = min(diff(bins11));
   d13 = min(diff(bins13));
   d15 = min(diff(bins15));
   d17 = min(diff(bins17));
   d19 = min(diff(bins19));
   d21 = min(diff(bins21));
   d25 = min(diff(bins25));
   binsize = [d9 d11 d13 d15 d17 d19 d21 d25];
   info = [info9 info11 info13 info15 info17 info19 info21 info25];

return;




%function mHist = hist2d ([vY, vX], vYEdge, vXEdge)
%2 Dimensional Histogram
%Counts number of points in the bins defined by vYEdge, vXEdge.
%size(vX) == size(vY) == [n,1]
%size(mHist) == [length(vYEdge) -1, length(vXEdge) -1]
%
%EXAMPLE
%   mYX = rand(100,2);
%   vXEdge = linspace(0,1,10);
%   vYEdge = linspace(0,1,20);
%   mHist2d = hist2d(mYX,vYEdge,vXEdge);
%
%   nXBins = length(vXEdge);
%   nYBins = length(vYEdge);
%   vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
%   vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
%   pcolor(vXLabel, vYLabel,mHist2d); colorbar
function mHist = hist2d (mX, vYEdge, vXEdge)
nCol = size(mX, 2);
if nCol < 2
    error ('mX has less than two columns')
end

nRow = length (vYEdge)-1;
nCol = length (vXEdge)-1;

vRow = mX(:,1);
vCol = mX(:,2);

mHist = zeros(nRow,nCol);

for iRow = 1:nRow
    rRowLB = vYEdge(iRow);
    rRowUB = vYEdge(iRow+1);
    
    vColFound = vCol((vRow > rRowLB) & (vRow <= rRowUB));
    
    if (~isempty(vColFound))
        
        
        vFound = histc(vColFound, vXEdge);
        
        nFound = (length(vFound)-1);
        
        if (nFound ~= nCol)
            disp([nFound nCol])
            error ('hist2d error: Size Error')
        end
        
        [nRowFound, nColFound] = size(vFound);
        
        nRowFound = nRowFound - 1;
        nColFound = nColFound - 1;
        
        if nRowFound == nCol
            mHist(iRow, :)= vFound(1:nFound)';
        elseif nColFound == nCol
            mHist(iRow, :)= vFound(1:nFound);
        else
            error ('hist2d error: Size Error')
        end
    end
    
end

return;



