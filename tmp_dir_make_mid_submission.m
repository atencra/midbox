function tmp_dir_make_mid_submission(stimfolder, exp, site)
% tmp_dir_mid_submission Create computer cluster queue scripts for MID analysis
%
% tmp_dir_mid_submission(stimfolder)
% Takes spike trains in the current folder and prepares them to be used in 
% MID analysis. The preparation involves creating a file with the binned spike
% train, creating parameter files for the MID analysis, and creating a 
% script that will be executed by the supercomputer cluster and will run
% the MID code.
%
% Spike trains are held in 'spk' struct arrays. The struct arrays are stored
% in  '*-site*-*um-*db-rn*-fs*-spk-strfcmb.mat' files.
%
% stimfolder is holds stimulus files that will be used to calculate MIDs.
% The folder holds downsampled versions of ripple stimuli. The stimuli
% were downsampled to allow the MID to calculate filters with a reasonable 
% amount of parameters. 
%
% For rat plasticity experiments, stimuli were downsampled to by a factor
% of 10 along the time axis and a factor of 8 along the frequency axis. This
% 5 ms resolution in time, and 6 frequencies per octave spectrally.
% 
%       stimfolder = 'D:\TM_Plasticity_Project_Physiology\Stimuli_Downsampled';
%
% Spike trains in the 'spk' struct array are downsampled to have the same 
% temporal resolution as the stimulus. After binning, the spike train is
% saved to an integer '.isk' file.
%
% One the process is completed, there will be one .isk file, four '.txt'
% files that hold parameters for the MID analysis. There are four files
% because their are four training/test data sets that are used. Finally,
% there is one '.sh' which gives the commands to the computer cluster.
%
% 

library('stimbox');
library('midbox');

narginchk(0,3);

if nargin == 0
    stimfolder = 'D:\TM_Plasticity_Project_Physiology\Stimuli_Downsampled';
end

if nargin >= 1 && nargin <= 3
    if isempty(stimfolder)
        stimfolder = 'D:\TM_Plasticity_Project_Physiology\Stimuli_Downsampled';
    end
end


if ~exist(stimfolder, 'dir')
    error(sprintf('stimfolder %s does not exist.', stimfolder));
end


if ( nargin < 2 ) 
    d1 = dir('*-site*-*um-*db-rn*-fs*-spk-strfcmb.mat'); % for single units
    d2 = [];
elseif ( nargin == 2 )
    if isempty(exp)
        d1 = dir('*-site*-*um-*db-rn*-fs*-spk-strfcmb.mat'); % for single units
        d2 = [];
    else
        d1 = dir(sprintf('%s-site*-*um-*db-rn*-fs*-spk-strfcmb.mat',exp)); 
        d2 = [];
    end
else ( nargin == 3 )
    if ( isempty(exp) && isempty(site) )
        d1 = dir('*-site*-*um-*db-rn*-fs*-spk-strfcmb.mat'); % for single units
        d2 = [];
    elseif ( isempty(exp) && ~isempty(site) )
        d1 = dir(sprintf('*-site%.0f-*um-*db-rn*-fs*-spk-strfcmb.mat',site)); 
        d2 = [];
    elseif ( ~isempty(exp) && isempty(site) )
        d1 = dir(sprintf('%s-site*-*um-*db-rn*-fs*-spk-strfcmb.mat', exp)); 
        d2 = [];
    else
        d1 = dir(sprintf('%s-site%.0f-*um-*db-rn*-fs*-spk-strfcmb.mat', exp, site)); 
        d2 = [];
    end
end

d = [d1(:)' d2(:)'];



dft = 10;
dff = 8;
t1 = 0.02;
t2 = 0.2;

fig = figure;

% Process every spk train file and every stimulus for every df(i) value.
for i = 1:length(d)

    fprintf('Processing %.0f of %.0f\n', i, length(d));

    infile = d(i).name;

    if ( ~isempty(findstr(infile,'spk-strfcmb')) )
        str = load(infile, 'spk', 'trigger');
        spk = str.spk;
        trigger = str.trigger;
        index = findstr(infile,'.mat');
        basename = infile(1:index-1);
        x0file = sprintf('%s-x0.mat', basename);
        str = load(x0file, 'x0');
        x0 = str.x0;
        if ( length(spk) ~= length(x0) )
            error('Need to choose x0 values.');
        end
    elseif ( ~isempty(findstr(infile,'thresh-strf-gd')) )
        str = load(infile, 'thresh', 'trigger');
        spk = str.thresh;
        trigger = str.trigger;
    else        
        error(sprintf('%s is not a recognized spike train file.',infile));
    end

    exp = spk(1).exp;
    site = spk(1).site;
    fss = spk(1).fs;
    stim = spk(1).stim;        

    % Get stimulus file name for the experiment
    [sprfile, paramfile, dsfile] = ...
        stim_get_rat_experiment_downsampled_specfile(exp, stim, dft, dff);

    sprfile_full = fullfile(stimfolder, sprfile);
    dsfile_full = fullfile(stimfolder, dsfile);

    % Check to make sure that the locator and the stimulus have ...
    % the same number of trials
    load(dsfile_full,'stimulus');

    for j = 1:length(spk)

        x0index = x0(j);

        spiketimes = spk(j).spiketimes;
        spet = round(spiketimes / 1000 * fss);


        % Get binned spike train
        [taxis, faxis, locator, numspikes, averate, locator_ipsi] = ...
           get_locator_for_mid_analysis(sprfile_full, t1, t2, spet, trigger, fss);


        if ( length(locator) ~= length(stimulus(1,:)) )
            error('spike train and stimulus must have same number of trials.');
        end


        % Write spike train to .isk file
        index = findstr(infile, '.mat');
        basename = infile(1:index-1);
        iskfile = sprintf('%s-%.0f.isk', basename, j);

        if ( ~exist(iskfile,'file') )
            fprintf('\n');
            fprintf('iskfile = %s\n', iskfile);
            fid = fopen(iskfile, 'w');
            count = fprintf(fid, '%u \n', locator);
            fclose(fid);
            fprintf('%s saved.\n', iskfile);
        else
            fprintf('Already saved: %s\n', iskfile);
        end



        % Organize parameters for submission script
        spkindex = j;
        NF = length(stimulus(:,1));
        nTrials = length(locator);
        nIterMid1 = 1000;
        nIterMid2 = 4000;
        nTestReps = 4;

        sta = get_sta_from_locator(locator, stimulus, 20);
        clf(fig); 
        imagesc(sta);
        hold on;
        plot([1 20], [x0index x0index], 'k-');
        plot([1 20], [x0index+25-1 x0index+25-1], 'k-');
        set(gcf,'position', [1233 555 560 420]);
        %imagesc(sta(x0:x0+25-1,:));
        pause(2);

        % Write submission script
        tmp_mid_qb3_submission_script(sprfile, iskfile, exp, site, stim, ...
            spkindex, NF, x0index, nTrials, nIterMid1, nIterMid2, nTestReps)

    end % (for j)

end % (for i)


return;



function tmp_mid_qb3_submission_script(sprfile, iskfile, exp, site, stim,...
    spkindex, NF, x0, nTrials, nIterMid1, nIterMid2, nTestReps)

narginchk(9,12);

if nargin == 9
    nIterMid1 = 1000;
    nIterMid2 = 4000;
    nTestReps = 4;
end

if nargin == 10
    nIterMid2 = 4000;
    nTestReps = 4;
end

if nargin == 11
    nTestReps = 4;
end


% At this point, we have the spr file name, we've written the 
% isk file to disk

if ( strcmp(stim, 'rn1') )
    part1 = '10';
elseif strcmp(stim,'rn11') 
    part1 = '11';
elseif strcmp(stim,'rn12')
    part1 = '12'; 
elseif ( strcmp(stim, 'rn4') )
    part1 = '40';
elseif ( strcmp(stim, 'rn41') )
    part1 = '41';
elseif ( strcmp(stim, 'rn42') )
    part1 = '42';
elseif ( strcmp(stim, 'rn8') )
    part1 = '80';
elseif ( strcmp(stim, 'rn81') )
    part1 = '81';
elseif ( strcmp(stim, 'rn82') )
    part1 = '82';
elseif ( strcmp(stim, 'rn16') )
    part1 = '160';
elseif ( strcmp(stim, 'rn161') )
    part1 = '161';
elseif ( strcmp(stim, 'rn162') )
    part1 = '162';
else
    error('Wrong file. RN must be 1, 4, 8, or 16');
end


% Create cell number so stimulus type is taken into account
if ( spkindex < 10 )
    part2 = sprintf('00%.0f', spkindex);
else
    part2 = sprintf('0%.0f', spkindex);
end
cellnum = str2double(sprintf('%s%s',part1,part2));


% Create experiment number so it has the form 20020619 or 20031029
index1 = findstr(exp,'-');
year = exp(1:index1(1)-1);
month = str2num(exp(index1(1)+1:index1(2)-1));
day = str2num(exp(index1(end)+1:end));

part1 = year;
if ( month < 10 )
    part2 = sprintf('0%.0f',month);
else
    part2 = sprintf('%.0f',month);
end

if ( day < 10 )
    part3 = sprintf('0%.0f',day);
else
    part3 = sprintf('%.0f',day);
end

expNum = str2double(sprintf('%s%s%s',part1,part2,part3));


%Creating submission scripts and parameter files
for testrep = 1:nTestReps

    %Creating parameter file
    midparampath = sprintf('rip_params_%d_%d_%d_%d.txt', expNum, site, cellnum, testrep);

    if ( ~exist(midparampath,'file') )
        fprintf('Creating parameter file: %s . . . ', midparampath);
    
        fid_midparam = fopen(midparampath,'w');
        fprintf(fid_midparam,'%s\n',sprfile);  %ripple file
        fprintf(fid_midparam,'%s\n',iskfile);  %isk file
        fprintf(fid_midparam,'%d\n',expNum);      %experiment date
        fprintf(fid_midparam,'%d\n',site);     %penetration site number
        fprintf(fid_midparam,'%d\n',cellnum);  %spike index number
        fprintf(fid_midparam,'%d\n',NF);       %number of frequency bands
        fprintf(fid_midparam,'%d\n',x0);       %dummy variable (used to be width of image)
        fprintf(fid_midparam,'%d\n',nTrials);  %number of samples in ripple
        fprintf(fid_midparam,'%d\n',1000);     %number of iterations for MID1
        fprintf(fid_midparam,'%d\n',4000);     %total number of iterations
        fprintf(fid_midparam,'%d\n',testrep);  %testrep index
        fclose(fid_midparam);
        fprintf('done.\n');
    else
        fprintf('Already created: %s . . . ', midparampath);
    end
end


%Creating submission script
shscriptpath = sprintf('script_%d_%d_%d_%d.sh', expNum, site, cellnum, testrep);

if ( ~exist(shscriptpath,'file') )
    fprintf('Creating submission script: %s . . . ', shscriptpath);
    
    fid_shscript = fopen(shscriptpath,'w');
    fprintf(fid_shscript,'#!/bin/bash\n');
    fprintf(fid_shscript,'#$ -S /bin/csh\n');
    fprintf(fid_shscript,'#$ -o /netapp/home/craiga/TMP/Logs\n');
    fprintf(fid_shscript,'#$ -e /netapp/home/craiga/TMP/Errors\n');
    fprintf(fid_shscript,'#$ -cwd\n');
    fprintf(fid_shscript,'#$ -r n\n');
    fprintf(fid_shscript,'#$ -j y\n');
    fprintf(fid_shscript,'#$ -l mem_free=1G\n');
    fprintf(fid_shscript,'#$ -l arch=linux-x64\n');
    fprintf(fid_shscript,'#$ -l netapp=1G,scratch=1G\n');
    fprintf(fid_shscript,'#$ -l h_rt=336:00:00\n');
    
    % edits to batch the 4 repetitions into one submission script. 
    fprintf(fid_shscript,'#$ -t 1-%d\n',nTestReps);
    fprintf(fid_shscript,'date\n');
    fprintf(fid_shscript,'hostname\n');
    
    fprintf(fid_shscript,...
    './rippledir1_detailed rip_params_%d_%d_%d_$SGE_TASK_ID.txt\n', expNum, site, cellnum);
    
    fprintf(fid_shscript,...
    './rippledir2_detailed rip_params_%d_%d_%d_$SGE_TASK_ID.txt\n', expNum, site, cellnum);

    fclose(fid_shscript);
else
    fprintf('Already created: %s . . . ', shscriptpath);
end


fprintf('done.\n');
        


return;

        






