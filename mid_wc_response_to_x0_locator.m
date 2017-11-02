function whole_cell_response_to_locator(sprpath)

if nargin == 0
    sprpath = 'E:\Ripple_Noise';
end


library('midbox');
library('strfbox');


d = dir('*_processed.mat');
pfile = d.name;
parts = strsplit(pfile, '_');
spkfile = sprintf('%s.mat', parts{1});

load(spkfile, 'spk', 'trigger');


% Change the following to downsampled version for mid analysis
% on C: drive on work computer
%sprfile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min.spr';
p = pwd;
[pathstr, name, ext] = fileparts(p);
index = findstr(name, '_');
stimtype = name(index(1)+1:index(2)-1);
sprfile = sprintf('%s-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt10_DFf5.spr', ...
    stimtype);

%sprfile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt10_DFf5.spr';
    
sprfile_path = fullfile(sprpath, sprfile);




strf = calculate_strf(spk, trigger, 0, sprfile_path); 
plot_strf_single(strf,trigger, [5000 40000]);



fsad = spk(1).fs;
spet = round(spk(1).spiketimes / 1000 * fsad);
[taxis, faxis, locator, numspikes, averate, ~] = ...
    get_locator_for_mid_analysis(sprfile_path, 0, 0.1, spet, trigger, fsad);


matrix_file = strrep(sprfile_path, '.spr', '-matrix.mat');
load(matrix_file, 'stimulus');

sta = get_sta_from_locator(locator, stimulus);

figure;
imagesc(sta);
size(stimulus)





% Make the file name for locator
exp = spk(1).exp;
site = spk(1).site;
depth = spk(1).depth;
atten = spk(1).atten;

index = findstr(sprfile, '-');
stim = sprfile(1:index(1)-1);


spkfile = sprintf('%s-site%.0f-%.0fum-%.0fdb-%s-fs%.0f-%u.isk', ...
    exp, site, depth, atten, stim, fsad, 1);

fprintf('spkfile = %s\n', spkfile);

% Write the locator to the file
fid = fopen( spkfile, 'w' );
count = fprintf( fid, '%u \n', locator );
fclose( fid );

fprintf( 'Ntrials = %.0f\n\n', length(locator) );


return;


















