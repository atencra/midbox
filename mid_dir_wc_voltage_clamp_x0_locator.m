function mid_dir_wc_current_clamp_x0_locator(stimtype, x0)

sprpath = 'C:\Ripple_Noise';


library('midbox');
library('strfbox');


% Change the following to downsampled version for mid analysis
% on C: drive on work computer
%sprfile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr';
if nargin == 0
    p = pwd;
    [pathstr, name, ext] = fileparts(p);
    index = findstr(name, '_');
    stimtype = name(index(1)+1:end);
end


sprfile = sprintf('%s-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt10_DFf5.spr', ...
    stimtype);

sprfile_path = fullfile(sprpath, sprfile);


fprintf('SPR file: %s\n', sprfile_path);





dspike = dir('*_th*.mat');

dtotal = [dspike(:)'];
files = {dtotal.name};

for i = 1:length(files)

    infile = files{i};
    load(infile, 'spk', 'trigger');

    if ( nargin < 2 )
        [x0, locator] = mid_get_x0_locator_for_mid_submission(sprfile_path, ...
                                                              spk, ...
                                                              trigger);
    else
        [~, locator] = mid_get_x0_locator_for_mid_submission(sprfile_path, ...
                                                             spk, ...
                                                             trigger, ...
                                                             x0);
    end


    % Make the file name for locator
    exp = spk(1).exp;
    site = spk(1).site;
    depth = spk(1).depth;
    atten = spk(1).atten;
    fsad = spk(1).fs;

    index = findstr(sprfile, '-');
    stim = sprfile(1:index(1)-1);


    if ( ~isempty(findstr(infile, 'epsps')) )
        index = 1;
        epsp = 1;
        spikes = 0;
    elseif ( ~isempty(findstr(infile, 'spikes')) )
        index = 2; 
        epsp = 0;
        spikes = 1;
    else
        index = 3; 
        epsp = 0;
        spikes = 0;

    end

    iskfile = sprintf('%s-site%.0f-%.0fum-%.0fdb-%s-fs%.0f-%u.isk', ...
        exp, site, depth, atten, stim, fsad, index);


    fprintf('iskfile = %s\n', iskfile);


    % Write the locator to the file
    fid = fopen( iskfile, 'w' );
    count = fprintf( fid, '%u \n', locator );
    fclose( fid );

    fprintf( 'Ntrials = %.0f\n\n', length(locator) );


    x0file = strrep(iskfile, '.isk', '-x0.mat');
    save(x0file, 'x0', 'epsp', 'spikes');

end



return;












