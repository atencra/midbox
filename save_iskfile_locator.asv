function locStr = save_iskfile_locator
% get_spkloc_from_isk_spk  Add locator to spk struct array 
% 
%    spkloc = get_spkloc_from_isk_spk(spk) reads in .isk files for
%    the corresponding elements in spk and returns a struct away
%    with two new fields: iskfile and locator. 
% 
%    iskfile is a binary file that holds the spike train of the 
%    neuron
% 
%    locator is a vector of 0s and 1s, corresponding to the spiketimes
% 
%    The data is coordinated was used by Tatyana Sharpee to estimate
%    MIDs in Atencio et al (2008) and Atencio et al (2009).
% 
%    spkloc now has the locator vector. To estimate an STA the
%    stimulus matrix may be read in and correlated with the
%    spike train
% 
%    Craig Atencio
%    4/10/12

% error(nargchk(1,1,nargin));


% Put channels in ascending order
% pos = [spk.position];
% [y, index] = sort(pos);
% spk = spk(index);
% 
% % Specs for recording
% exp = spk(1).exp;
% site = spk(1).site;

% isk files for recording
% d = dir( sprintf('%s-site%.0f-*.isk', exp, site) );
d = dir( sprintf('*.isk') );

% if ( length(d) ~= length(spk) )
%    error('Mismatch between # isk files and # spk elements.');
% end
% 
% spkloc = spk;

for i = 1:length(d)
   
   iskfile = d(i).name;
   index_dash = findstr(iskfile, '-');
   index_dot = findstr(iskfile, '.');

   unit = str2double ( iskfile(index_dash(end)+1:index_dot-1) );

   locator = locator_from_isk_file(iskfile);

   locStr(i).iskfile = iskfile;
   locStr(i).locator = locator;

end % (for i)



return;






