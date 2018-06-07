function filestruct = get_files_for_mid_analysis(location)
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

if ( nargin~=1 )
   error('You need one input arg.');
end

current_directory = pwd;

if ( isempty(findstr(current_directory, 'C:\MATLAB65\work\tatyana\Filters')) )
   error('You need to run this function in C:\MATLAB65\work\tatyana\Filters');
end

dfile = dir(['rpsta_' num2str(location) '*.dat']);

if ( isempty(dfile) )
   error('That location does not exist.');
end

unitnum = [];
for i = 1:length(dfile)
   file = dfile(i).name;
   index = findstr('_', file);
   unitnum = [unitnum str2num(file(index(2)+1:index(3)-1))];
end
unitnum = sort(unique(unitnum));

temp = dfile(1).name;
index = findstr(file,'x');
fbins = str2num( file(index(1)+1:index(2)-1) );
tbins = str2num( file(index(2)+1:index(2)+2) );


% Define a struct to hold the results
filestruct = struct(...
'location',      [], ...
'unit',          [], ...
'tbins',         [], ...
'fbins',         [], ...
'rpsta',     [], ...
'rpx1pxpxt_sta', [], ...
'rpdbest1_v1',    [], ...
'rpdbest1_v2',    [], ...
'rpdbest2_v1',    [], ...
'rpdbest2_v2',    [], ...
'rpdtest1_v1',    [], ...
'rpdtest1_v2',    [], ...
'rpdtest2_v1',    [], ...
'rpdtest2_v2',    [], ...
'rpdx1x2px_pxt_1',   [], ...
'rpdx1x2px_pxt_2',   []);


for i = 1:length(unitnum)

   filestruct(i).location = location;
   filestruct(i).unit = unitnum(i);
   filestruct(i).tbins = tbins;
   filestruct(i).fbins = fbins;


   % The RPSTA Files
   %=================================================================

   % Get the STA files
   dfile = dir(['rpsta_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpsta{1} = dfile(1).name;
      filestruct(i).rpsta{2} = dfile(2).name;
      filestruct(i).rpsta{3} = dfile(3).name;
      filestruct(i).rpsta{4} = dfile(4).name;
   end


   % Get the STA input/output function files
   dfile = dir(['rpx1pxpxt_sta_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpx1pxpxt_sta{1} = dfile(1).name;
      filestruct(i).rpx1pxpxt_sta{2} = dfile(2).name;
      filestruct(i).rpx1pxpxt_sta{3} = dfile(3).name;
      filestruct(i).rpx1pxpxt_sta{4} = dfile(4).name;
   end


   % The RPDTEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdtest1_v1{1} = dfile(1).name;
      filestruct(i).rpdtest1_v1{2} = dfile(2).name;
      filestruct(i).rpdtest1_v1{3} = dfile(3).name;
      filestruct(i).rpdtest1_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdtest1_v2{1} = dfile(1).name;
      filestruct(i).rpdtest1_v2{2} = dfile(2).name;
      filestruct(i).rpdtest1_v2{3} = dfile(3).name;
      filestruct(i).rpdtest1_v2{4} = dfile(4).name;
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdtest2_v1{1} = dfile(1).name;
      filestruct(i).rpdtest2_v1{2} = dfile(2).name;
      filestruct(i).rpdtest2_v1{3} = dfile(3).name;
      filestruct(i).rpdtest2_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdtest2_v2{1} = dfile(1).name;
      filestruct(i).rpdtest2_v2{2} = dfile(2).name;
      filestruct(i).rpdtest2_v2{3} = dfile(3).name;
      filestruct(i).rpdtest2_v2{4} = dfile(4).name;
   end



   % The RPDBEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdbest1_v1{1} = dfile(1).name;
      filestruct(i).rpdbest1_v1{2} = dfile(2).name;
      filestruct(i).rpdbest1_v1{3} = dfile(3).name;
      filestruct(i).rpdbest1_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdbest1_v2{1} = dfile(1).name;
      filestruct(i).rpdbest1_v2{2} = dfile(2).name;
      filestruct(i).rpdbest1_v2{3} = dfile(3).name;
      filestruct(i).rpdbest1_v2{4} = dfile(4).name;
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdbest2_v1{1} = dfile(1).name;
      filestruct(i).rpdbest2_v1{2} = dfile(2).name;
      filestruct(i).rpdbest2_v1{3} = dfile(3).name;
      filestruct(i).rpdbest2_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdbest2_v2{1} = dfile(1).name;
      filestruct(i).rpdbest2_v2{2} = dfile(2).name;
      filestruct(i).rpdbest2_v2{3} = dfile(3).name;
      filestruct(i).rpdbest2_v2{4} = dfile(4).name;
   end



   % Get the MID input/output function files - First and Second Most
   % Informative Dimensions
   dfile = dir(['rpdx1x2px_pxt_1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
      filestruct(i).rpdx1x2px_pxt_1{1} = dfile(1).name;
      filestruct(i).rpdx1x2px_pxt_1{2} = dfile(2).name;
      filestruct(i).rpdx1x2px_pxt_1{3} = dfile(3).name;
      filestruct(i).rpdx1x2px_pxt_1{4} = dfile(4).name;
   end


   dfile = dir(['rpdx1x2px_pxt_2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) & length(dfile) == 4 )
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





