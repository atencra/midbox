function [filestruct] = get_files_for_mid_analysis(stimfolder)
% get_files_for_mid_analysis - put MID data files into a struct array
%    so data can be analysed
%
% filestruct = get_files_for_mid_analysis(location)
%
% The folder C:\MATLAB65\work\tatyana\Filters contains all the data from
% Tatyana's maximally informative dimensions analysis. This function goes
% through that folder, extracts all the files for one location, and then
% creates a struct array of the file information, where each element of the
% struct array corresponds to one single unit in the location.
%
% location : a scalar. May be 500, 515, 516, 517, 519, 532, 534, 537, 602,
% 603, 604, 608, 609, 611, 614, or 616.
%
%
% caa 4/11/06


if ~exist('stimfolder','var')
    %stimfolder = 'I:\Ripple_Noise\downsampled_for_MID';
end

files = gfn('rpsta_*', 0);
re = regexp(files, '(?<=(rpsta_))\S+(?=(_\dx\d{2}x\d{2}))', 'match','once');
uniquefiles = unique(re);
details = regexp(uniquefiles, '_', 'split');
exp = cellfun(@(x) x{1}, details, 'UniformOutput', 0);
site = cellfun(@(x) x{2}, details, 'UniformOutput', 0);
unitnum = cellfun(@(x) x{3}, details, 'UniformOutput', 0);
    

temp = files{1};
splitexp = regexp(temp, 'x', 'split');
fbins = str2double(splitexp{2});
tbins = str2double(regexp(splitexp{3}, '^\d{1,3}', 'match','once'));

filestruct(length(unitnum)).exp = [];

for i = 1:length(unitnum)
    
    
    % get x0 from input .txt files
    txtfile = gfn(['*' exp{i} '*' site{i} '_' unitnum{i} '_1.txt'], 0);
    f = fopen(txtfile, 'r');
    C = textscan(f, '%s', 'Delimiter', '\n');
    fclose('all');
    x0 = str2double(C{1}{7});
    stim = regexp(C{1}{1}, '^\w+', 'match', 'once');      
    %sprfile = fullfile(stimfolder, C{1}{1});
    sprfile = C{1}{1};
    iskfile = C{1}{2};
    locator = isk_file_to_locator(iskfile);


    % Define a struct to hold the results
    % filestruct = struct(...
    % 'location',      [], ...
    % 'site',          [], ...
    % 'unit',          [], ...
    % 'tbins',         [], ...
    % 'fbins',         [], ...
    % 'rpsta',     [], ...
    % 'rpx1pxpxt_sta', [], ...
    % 'rpdbest1_v1',    [], ...
    % 'rpdbest1_v2',    [], ...
    % 'rpdbest2_v1',    [], ...
    % 'rpdbest2_v2',    [], ...
    % 'rpdtest1_v1',    [], ...
    % 'rpdtest1_v2',    [], ...
    % 'rpdtest2_v1',    [], ...
    % 'rpdtest2_v2',    [], ...
    % 'rpdx1x2px_pxt_1',   [], ...
    % 'rpdx1x2px_pxt_2',   []);



   filestruct(i).exp = exp{i};
   filestruct(i).site = site{i};
   filestruct(i).unit = unitnum{i};
   filestruct(i).stim = stim;
   filestruct(i).sprfile = sprfile;
   filestruct(i).iskfile = iskfile;
   filestruct(i).locator = locator;
   filestruct(i).tbins = tbins;
   filestruct(i).fbins = fbins;
   filestruct(i).nh = 25;
   filestruct(i).nv = 1;
   filestruct(i).nlags = tbins;
   filestruct(i).x0 = x0;
   
%    txtfile = gfn([fullfile(currfold, '\Inputs\') '*' exp '*_0' num2str(site) '_' num2str(i) '_1.txt'], 0);
%    f = fopen(fullfile('Inputs\',txtfile{1}), 'r');
%    C = textscan(f, '%s', 'Delimiter', '\n');
%    label = regexp(C{1}{2},'(?<=(spk-))\w*_\d{1,3}(?=(.isk$))','match','once');
%    filestruct(i).label = regexprep(label, '_', '-');
%    fclose('all');


   % The RPSTA Files
   %=================================================================

   % Get the STA files
   dfile = dir(['rpsta_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpsta{1} = dfile(1).name;
      filestruct(i).rpsta{2} = dfile(2).name;
      filestruct(i).rpsta{3} = dfile(3).name;
      filestruct(i).rpsta{4} = dfile(4).name;
   end


   % Get the STA input/output function files
   dfile = dir(['rpx1pxpxt_sta_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpx1pxpxt_sta{1} = dfile(1).name;
      filestruct(i).rpx1pxpxt_sta{2} = dfile(2).name;
      filestruct(i).rpx1pxpxt_sta{3} = dfile(3).name;
      filestruct(i).rpx1pxpxt_sta{4} = dfile(4).name;
   end


   % The RPDTEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v1_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdtest1_v1{1} = dfile(1).name;
      filestruct(i).rpdtest1_v1{2} = dfile(2).name;
      filestruct(i).rpdtest1_v1{3} = dfile(3).name;
      filestruct(i).rpdtest1_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v2_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdtest1_v2{1} = dfile(1).name;
      filestruct(i).rpdtest1_v2{2} = dfile(2).name;
      filestruct(i).rpdtest1_v2{3} = dfile(3).name;
      filestruct(i).rpdtest1_v2{4} = dfile(4).name;
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v1_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdtest2_v1{1} = dfile(1).name;
      filestruct(i).rpdtest2_v1{2} = dfile(2).name;
      filestruct(i).rpdtest2_v1{3} = dfile(3).name;
      filestruct(i).rpdtest2_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v2_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdtest2_v2{1} = dfile(1).name;
      filestruct(i).rpdtest2_v2{2} = dfile(2).name;
      filestruct(i).rpdtest2_v2{3} = dfile(3).name;
      filestruct(i).rpdtest2_v2{4} = dfile(4).name;
   end



   % The RPDBEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v1_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdbest1_v1{1} = dfile(1).name;
      filestruct(i).rpdbest1_v1{2} = dfile(2).name;
      filestruct(i).rpdbest1_v1{3} = dfile(3).name;
      filestruct(i).rpdbest1_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v2_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdbest1_v2{1} = dfile(1).name;
      filestruct(i).rpdbest1_v2{2} = dfile(2).name;
      filestruct(i).rpdbest1_v2{3} = dfile(3).name;
      filestruct(i).rpdbest1_v2{4} = dfile(4).name;
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v1_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdbest2_v1{1} = dfile(1).name;
      filestruct(i).rpdbest2_v1{2} = dfile(2).name;
      filestruct(i).rpdbest2_v1{3} = dfile(3).name;
      filestruct(i).rpdbest2_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v2_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdbest2_v2{1} = dfile(1).name;
      filestruct(i).rpdbest2_v2{2} = dfile(2).name;
      filestruct(i).rpdbest2_v2{3} = dfile(3).name;
      filestruct(i).rpdbest2_v2{4} = dfile(4).name;
   end



   % Get the MID input/output function files - First and Second Most
   % Informative Dimensions
   dfile = dir(['rpdx1x2px_pxt_1_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdx1x2px_pxt_1{1} = dfile(1).name;
      filestruct(i).rpdx1x2px_pxt_1{2} = dfile(2).name;
      filestruct(i).rpdx1x2px_pxt_1{3} = dfile(3).name;
      filestruct(i).rpdx1x2px_pxt_1{4} = dfile(4).name;
   end


   dfile = dir(['rpdx1x2px_pxt_2_' exp{i} '_' site{i} '_' unitnum{i} '_*.dat']);
   if ( ~isempty(dfile) && length(dfile) == 4 )
      filestruct(i).rpdx1x2px_pxt_2{1} = dfile(1).name;
      filestruct(i).rpdx1x2px_pxt_2{2} = dfile(2).name;
      filestruct(i).rpdx1x2px_pxt_2{3} = dfile(3).name;
      filestruct(i).rpdx1x2px_pxt_2{4} = dfile(4).name;
   end

end % (for ii)



% row=3;
% start_column=1;
% fname_first=sprintf('rpdtest%u_v1_%u_%u_1x%ux%u_%u',optlevel, location, cell, Nv, Nh, cy)
% Nparts=4;
% [v1,coeff1,projection1,mtx1] = plot_a_vector_audSTRF(fname_first,Nh,Nv,nlags,Nwiny,Nwinx,row,start_column,between_column,nlag_start,nlag_end,Nparts);
% axis tight;
% ylabel('MID_1');
% 
% fname_first=sprintf('rpdtest%u_v2_%u_%u_1x%ux%u_%u',optlevel,location,cell,Nv,Nh,cy);
% row=4;
% [v2,coeff2,projection1,mtx2]=plot_a_vector_audSTRF(fname_first,Nh,Nv,nlags,Nwiny,Nwinx,row,start_column,...
%     between_column,nlag_start,nlag_end,Nparts);
% axis tight;
% ylabel('MID_2');
% 
% fname=sprintf('rpdx1x2px_pxt_%u_%u_%u_1x%ux%u_%u',optlevel,location,cell,Nv,Nh,cy);
% start_column=2;
% row=3;
% Nbins_medium=15;
% plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,nlags+1,coeff1,coeff2,Nbins_medium,Nparts,Nbins_short);
% 
% row=4;
% plot_an_1db_from2dior_improved(fname,Nwiny,Nwinx,row,nlags+1,coeff1,coeff2,Nbins_medium,Nparts,Nbins_short);
% xlabel('projection in units of standard deviation')
% h=gca;
% 
% row=2;
% [indx,indy]=plot_an_2dior_improvedjust(fname,Nwiny/2,Nwinx,row,nlags+2,coeff1,coeff2,Nparts,Nbins_medium);
% title('P(spike|2D proj)')
% 
% colormap('hot')
% 
% set(hfig,'PaperOrientation','landscape');
% set(hfig,'PaperPosition',[0.25 0.614603 10.5 7.27079])





