function whole_cell_mid_submission(sprpath)
% Cluster MID script for kyunghee's whole cell data
% 
%    whole_cell_mid_submission(sprpath)
% 
%        sprpath : directory where spr files are stored.
%        .spr, _param.mat, and -matrix.mat files must be
%        in this location.


library('midbox');
library('strfbox');


if nargin == 0
    sprpath = 'C:\Ripple_Noise';
end

d = dir('*_x0.mat');
load(d.name, 'x0');

d = dir('*.isk');
iskfile = d.name;


index_dash = findstr(iskfile, '-');
index_site = findstr(iskfile, 'site');
index_isk = findstr(iskfile, '.isk');

exp = str2num(iskfile(1:index_dash(1)-1));
site = str2num(iskfile(index_site+4:index_dash(2)-1));
cellnum = str2num(iskfile(index_dash(end)+1:index_isk-1));
stimtype = iskfile(index_dash(4)+1:index_dash(5)-1); 

sprfile = sprintf('%s-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt10_DFf5.spr', stimtype);
sprfile_path = fullfile(sprpath, sprfile);


% Load ripple, get dimensions
[stimulus] = SprFile2Matrix(sprfile_path);
[NF, nTrials] = size(stimulus);


%Creating submission scripts and parameter files
nTestReps = 4;
for testrep = 1:nTestReps

    %Creating parameter file
    midparampath = sprintf('rip_params_%d_%d_%d_%d.txt', exp, site, cellnum, testrep);

    if ( ~exist(midparampath,'file') )
        fprintf('Creating parameter file: %s . . . \n', midparampath);
    
        fid_midparam = fopen(midparampath,'w');
        fprintf(fid_midparam,'%s\n',sprfile);  %ripple file
        fprintf(fid_midparam,'%s\n',iskfile);  %isk file
        fprintf(fid_midparam,'%d\n',exp);   %experiment date
        fprintf(fid_midparam,'%d\n',site);     %penetration site number
        fprintf(fid_midparam,'%d\n',cellnum);  %spike index number
        fprintf(fid_midparam,'%d\n',NF);       %number of frequency bands
        fprintf(fid_midparam,'%d\n',x0);       % index into spr matrix
        fprintf(fid_midparam,'%d\n',nTrials);  %number of samples in ripple
        fprintf(fid_midparam,'%d\n',1000);     %number of iterations for MID1
        fprintf(fid_midparam,'%d\n',4000);     %total number of iterations
        fprintf(fid_midparam,'%d\n',testrep);  %testrep index
        fclose(fid_midparam);
        fprintf('done.\n');
        fprintf('Created: %s \n', midparampath);
    else
        fprintf('Already created: %s . . . \n', midparampath);
    end
end


%Creating submission script
shscriptpath = sprintf('script_%d_%d_%d_%d.sh', exp, site, cellnum, testrep);

if ( 1 ) %~exist(shscriptpath,'file') )

    fprintf('Creating submission script: %s . . . \n', shscriptpath);
    
    fid_shscript = fopen(shscriptpath,'w');
    fprintf(fid_shscript,'#!/bin/bash\n');
    fprintf(fid_shscript,'#$ -S /bin/bash\n');
    fprintf(fid_shscript,'#$ -o /netapp/home/craiga/WholeCell/Logs\n');
    fprintf(fid_shscript,'#$ -e /netapp/home/craiga/WholeCell/Errors\n');
    fprintf(fid_shscript,'#$ -cwd\n');
    fprintf(fid_shscript,'#$ -r n\n');
    fprintf(fid_shscript,'#$ -j y\n');
    fprintf(fid_shscript,'#$ -l mem_free=1G\n');
    fprintf(fid_shscript,'#$ -l arch=linux-x64\n');
    fprintf(fid_shscript,'#$ -l netapp=1G,scratch=1G\n');
    fprintf(fid_shscript,'#$ -l h_rt=336:00:00\n');
    
    fprintf(fid_shscript,'#$ -t 1-%d\n',nTestReps);
    fprintf(fid_shscript,'date\n');
    fprintf(fid_shscript,'hostname\n');
    
    fprintf(fid_shscript, './rippledir1_detailed rip_params_%d_%d_%d_$SGE_TASK_ID.txt\n', exp, site, cellnum);
    fprintf(fid_shscript, './rippledir2_detailed rip_params_%d_%d_%d_$SGE_TASK_ID.txt\n', exp, site, cellnum);
    
    fclose(fid_shscript);
    fprintf('Created: %s \n', shscriptpath);
    
else
    fprintf('Already created: %s . . . \n', shscriptpath);
end

fprintf('done.\n\n');


return


function [stimulus] = SprFile2Matrix(sprfile)
%[stimulus] = SprFile2Matrix(sprfile)
%
%
%       FILE NAME       : RT WSTRF DB
%       DESCRIPTION     : Real Time spectro-temporal receptive field
%			  Uses Lee/Schetzen Aproach via Specto-Temporal Envelope
%			  For dB Amplitude Sound distributions 
%
% SpecFile : Spectral Profile File
% sprtype : SPR File Type : 'float' or 'int16'
%         Default=='float'. This function assumes it is a float file.	
%
% stimulus : mxn matrix, which is the entire dmr envelope file.
%
% To get things to match up with Tatyana's code you use the following line
% to compute the STA:
%
% [taxis,faxis,STRF1,PP,Wo1,No1,SPLN] = get_sta_filter(file, 0.005, 0.095, spet, trigger, fs, SPL, MdB, tbins, fbins);
%
% tbins = 20 and fbins = 25 in every case
%
% caa 12/15/06

Sound = 'MR';
ModType = 'dB';
sprtype='float';


%Loading Parameter Data
index = findstr(sprfile,'.spr');
paramFile = [sprfile(1:index(1)-1) '_param.mat'];
f = ['load(''' paramFile ''')'];
eval(f);
clear App  MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 f2 
clear Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn

 
%Opening Spectral Profile File
fid = fopen(sprfile);
frewind(fid);
MdB = 40;
RMSP = -MdB/2;
stimulus = [];
while ( ~feof(fid) )
   [s1, count] = fread(fid,NT*NF,'float');
   if ( count )
      s1 = reshape(s1,NF,NT);
      s1 = MdB * s1 - RMSP;
      stimulus = [stimulus s1];
   end % (if)
end % (while)

[nf, ntrials] = size(stimulus);

fprintf('\n#Freqs = %.0f, #trials = %.0f\n\n', nf, ntrials);

%Closing all opened files
fclose('all');

return;



function [taxis, faxis, locator, numspikes, averate] = ...
   SpetSpr2Locator(specfile, t1, t2, spet, trigger, fss)
% [taxis, faxis, locator, numspikes, averate] = ...
%     SpetSpr2Locator(specfile, t1, t2, spet, trigger, fss)
%
% Obtain spike train locator for Maximally Informative Dimensions Analysis.
%
%  specfile       :  Spectral Profile File ( *.spr file )
%  t1, t2         :  Evaluation delay interval for STRF(T,F), in seconds
%                    T E [- T1 , T2 ], Note that T1 and T2 > 0
%  spet           :  Array of spike event times in sample number
%  trigger        :  Array of Trigger Times
%  fss            :  Sampling Rate for TRIGGER and SPET
%
%  RETURNED VALUES 
%
%  taxis      : Time Axis (is sec)
%  faxis      : Frequency Axis (in Hz)
%  locator    : spike train binned at the resolution of specfile
%  numspikes  : Number of spikes used in strf
%  averate    : average firing rate (spikes / sec)
%
%
% To obtain the output variables two files need to be placed in the same
% folder. One is the specfile, which is named as *.spr. The other file
% is a parameter file that is named as *_param.mat file. 
%
%
% Here specfile is the spr envelope file, 0.05 specifies the noncausal
% extent of the strf, in seconds, and 0.2 specifies the causal portion 
% of the strf, in seconds. For a faster estimate with fewer dimensions the
% values may be changed to 0.0 and 0.080.
%
% spet is the spike times in sample number, created via the
% following command: spet = spk(i).spiketimes / 1000 * fss, where spk is
% a data structure holding spike time information, with spike times in
% msec's.
%
% trigger is a vector of trigger times, in sample number.
% The triggers are included because the sampling rate of the sound and the 
% spike times are different, and also because the envelope file has been
% downsampled. The triggers help us line up the spike times with the
% appropriate bins of the envelope, from which we estimate get the 
% binned spike train.

numtrig = length(trigger);

% Loading Parameter Data
index = findstr(specfile,'.spr');
paramfile = [specfile(1:index(1)-1) '_param.mat'];
f = ['load ' paramfile];
eval(f);

% Clear possible useless variables in the param file
clear App MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 
clear f2 Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn


% Flipping Trig and Spet for channel 1 STRF
mintime = min([trigger spet]);
maxtime = max([trigger spet]);
spet = spet - mintime + 1;
trigger = trigger - mintime + 1;

% Converting Temporal Delays to Sample Numbers
n1 = round( t1 * Fs / DF);
n2 = round( t2 * Fs / DF);

% % Initializing Some Variables
% numspikes = 0;  % Number of Spikes

Ntrials = NT*numtrig; % newnt=320;

locator = zeros(Ntrials,1);

for trigcount = 2:length(trigger)-1 

    % Finding SPET in between triggers and resampling spet relative to the 
    %    Spectral Profile samples
    index = find( spet>=trigger(trigcount) & spet<trigger(trigcount+1) );
    spettrig = ceil( (spet(index)-trigger(trigcount)+1) * Fs / fss / DF );
   

	index2 = find(spet>trigger(numtrig-trigcount+1) & spet<=trigger(numtrig-trigcount+2));
	spettrig2 = ceil( (trigger(numtrig-trigcount+2)+1-spet(index2)) * Fs / fss /DF );

    for k = 1:length(spettrig)
        spike = spettrig(k);
        locator(spike+(trigcount-1)*NT) = locator(spike+(trigcount-1)*NT)+1;
    end
    
end % (while ~feof(fid) & trigcount<length(trigger)-1 )

% Closing all opened files
fclose('all');

[n, nt, nf] = ripple_stim_length(specfile);

index_min = min( [length(locator) n/nf] );
locator = locator( 1:index_min );

numspikes = sum(locator);
stim_duration = ( max(trigger) - min(trigger) ) / fss;
averate = numspikes / stim_duration;

% get time and frequency axes for the receptive field
faxis = faxis;
taxis = (-n1:n2-1) / (Fs/DF);

return;



function [n, NT, NF] = ripple_stim_length(specfile)
% n = ripple_stim_length(specfile)
%
% Determine length of ripple stimulus - How many time trials were there? 
%
%  INPUT VALUES
%
%  specfile       :  Spectral Profile File ( *.spr file )
%
%
%  RETURNED VALUES 
%
%  n : number of envelope samples
%  NT : Number of time samples in a ripple envelope segment
%  NF : Number of freq samples in a ripple segment


% Loading Parameter Data
index = findstr(specfile,'.spr');
paramfile = [specfile(1:index(1)-1) '_param.mat'];
f = ['load ' paramfile];
eval(f);

% Clear possible useless variables in the param file
clear App MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 
clear f2 Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn

fid = fopen(specfile); % open ripple envelop file
frewind(fid); % make sure the file is set to the beginning

n = 0; 

while ( ~feof(fid) )
      [env, count] = fread(fid, NT * NF, 'float');
      n = n + count;
end % ( while ~feof(fid) )

% Closing all opened files
fclose('all');

return;





