function filestruct = get_rippledir_detailed_mid_icc_files(location, x0, nh, nlags, strf)
% get_files_for_mid_icc_analysis - put MID data files into a struct array
%    so data can be analysed
%
% filestruct = get_rippledir_detailed_mid_icc_files(location, x0, nh, nlags, strf)
%
% The folder C:\MATLAB65\work\tatyana\Filters contains all the data from
% Tatyana's maximally informative dimensions analysis. This function goes
% through that folder, extracts all the files for one location, and then
% creates a struct array of the file information, where each element of the
% struct array corresponds to one single unit in the location.
%
% location, for auditory cortex : a scalar. May be 500, 515, 516, 517, 519, 
% 532, 534, 537, 602, 603, 604, 608, 609, 611, 614, or 616.
%
% x0 : starting frequency bin for mid analysis.
%
% For ICC:
%
%     Experiment/Site         Location       x0
%     -----------------------------------------
%     2003-11-24 site7        707            20
%     2003-11-24 site8        708            20
%     2003-11-24 site9        709            20
%     2003-11-24 site14       714            35
%     2003-11-19 site1        701            1
%
% nh = 25 and nlags = 20, in all cases.
%
% strf : strf data from which the original middata were originally extracted
%
% caa 1/28/10

filestruct = [];

if ( nargin ~= 5 )
   error('You need 5 input args.');
end

if ( nargin == 2 )
   nh = 25;
   nlags = 20;
end

if ( nargin == 3 )
   nlags = 20;
end



nv = 1;


% dfile = dir(['rpsta_' num2str(location) '*.dat']);

% If rpdtest2_v2 exists, then we have run the code for two mids
dfile = dir(['rpdtest2_v2_' num2str(location) '*.dat']);

% rpdtest2_v2_707_1_1x25x20_1_1.dat

if ( isempty(dfile) )
   error('That location does not exist.');
end

unitnum = [];
for i = 1:length(dfile)
   file = dfile(i).name;
   index = findstr('_', file);
   unitnum = [unitnum str2num(file(index(3)+1:index(4)-1))];
end
unitnum = sort(unique(unitnum));



temp = dfile(1).name;
index = findstr(file,'x');
numfbins = str2num( file(index(1)+1:index(2)-1) );
numtbins = str2num( file(index(2)+1:index(2)+2) );


% Define a struct to hold the results
filestruct = struct(...
'exp',             [], ...
'site',            [], ...
'chan',            [], ...
'model',           [], ...
'depth',           [], ...
'position',        [], ...
'stim',            [], ...
'atten',           [], ...
'spl',             [], ...
'sm',              [], ...
'tm',              [], ...
'mdb',             [], ...
'location',        [], ...
'unit',            [], ...
'x0',              [], ...
'nh',              [], ...
'nv',              [], ...
'nlags',           [], ...
'numtbins',        [], ...
'numfbins',        [], ...
'rpsta',           [], ...
'rpx1pxpxt_sta',   [], ...
'rpdbest1_v1',     [], ...
'rpdbest1_v2',     [], ...
'rpdbest2_v1',     [], ...
'rpdbest2_v2',     [], ...
'rpdtest1_v1',     [], ...
'rpdtest1_v2',     [], ...
'rpdtest2_v1',     [], ...
'rpdtest2_v2',     [], ...
'rpdx1x2px_pxt_1', [], ...
'rpdx1x2px_pxt_2', []);


for i = 1:length(unitnum)

   unit = unitnum(i);

   filestruct(i).exp = strf(unit).exp;
   filestruct(i).site = strf(unit).site;
   filestruct(i).chan = strf(unit).chan;
   filestruct(i).model = strf(unit).model;
   filestruct(i).depth = strf(unit).depth;
   filestruct(i).position = strf(unit).position;
   filestruct(i).stim = strf(unit).stim;
   filestruct(i).atten = strf(unit).atten;
   filestruct(i).spl = strf(unit).spl;
   filestruct(i).sm = strf(unit).sm;
   filestruct(i).tm = strf(unit).tm;
   filestruct(i).mdb = strf(unit).mdb;

   filestruct(i).location = location;
   filestruct(i).unit = unit;

   filestruct(i).x0 = x0;

   filestruct(i).nh = nh;
   filestruct(i).nv = nv;
   filestruct(i).nlags = nlags;

   filestruct(i).numtbins = numtbins;
   filestruct(i).numfbins = numfbins;


   % The RPSTA Files
   %=================================================================

   % Get the STA files
   dfile = dir(['rpsta_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpsta{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpsta{1} = dfile(1).name;
%       filestruct(i).rpsta{2} = dfile(2).name;
%       filestruct(i).rpsta{3} = dfile(3).name;
%       filestruct(i).rpsta{4} = dfile(4).name;
   end


   % Get the STA input/output function files
   dfile = dir(['rpx1pxpxt_sta_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpx1pxpxt_sta{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpx1pxpxt_sta{1} = dfile(1).name;
%       filestruct(i).rpx1pxpxt_sta{2} = dfile(2).name;
%       filestruct(i).rpx1pxpxt_sta{3} = dfile(3).name;
%       filestruct(i).rpx1pxpxt_sta{4} = dfile(4).name;
   end


   % The RPDTEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest1_v1{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdtest1_v1{1} = dfile(1).name;
%       filestruct(i).rpdtest1_v1{2} = dfile(2).name;
%       filestruct(i).rpdtest1_v1{3} = dfile(3).name;
%       filestruct(i).rpdtest1_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest1_v2{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdtest1_v2{1} = dfile(1).name;
%       filestruct(i).rpdtest1_v2{2} = dfile(2).name;
%       filestruct(i).rpdtest1_v2{3} = dfile(3).name;
%       filestruct(i).rpdtest1_v2{4} = dfile(4).name;
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest2_v1{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdtest2_v1{1} = dfile(1).name;
%       filestruct(i).rpdtest2_v1{2} = dfile(2).name;
%       filestruct(i).rpdtest2_v1{3} = dfile(3).name;
%       filestruct(i).rpdtest2_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest2_v2{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdtest2_v2{1} = dfile(1).name;
%       filestruct(i).rpdtest2_v2{2} = dfile(2).name;
%       filestruct(i).rpdtest2_v2{3} = dfile(3).name;
%       filestruct(i).rpdtest2_v2{4} = dfile(4).name;
   end



   % The RPDBEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest1_v1{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdbest1_v1{1} = dfile(1).name;
%       filestruct(i).rpdbest1_v1{2} = dfile(2).name;
%       filestruct(i).rpdbest1_v1{3} = dfile(3).name;
%       filestruct(i).rpdbest1_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest1_v2{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdbest1_v2{1} = dfile(1).name;
%       filestruct(i).rpdbest1_v2{2} = dfile(2).name;
%       filestruct(i).rpdbest1_v2{3} = dfile(3).name;
%       filestruct(i).rpdbest1_v2{4} = dfile(4).name;
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest2_v1{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdbest2_v1{1} = dfile(1).name;
%       filestruct(i).rpdbest2_v1{2} = dfile(2).name;
%       filestruct(i).rpdbest2_v1{3} = dfile(3).name;
%       filestruct(i).rpdbest2_v1{4} = dfile(4).name;
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest2_v2{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdbest2_v2{1} = dfile(1).name;
%       filestruct(i).rpdbest2_v2{2} = dfile(2).name;
%       filestruct(i).rpdbest2_v2{3} = dfile(3).name;
%       filestruct(i).rpdbest2_v2{4} = dfile(4).name;
   end


   % Get the MID input/output function files - First and Second Most
   % Informative Dimensions
   dfile = dir(['rpdx1x2px_pxt_1_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdx1x2px_pxt_1{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdx1x2px_pxt_1{1} = dfile(1).name;
%       filestruct(i).rpdx1x2px_pxt_1{2} = dfile(2).name;
%       filestruct(i).rpdx1x2px_pxt_1{3} = dfile(3).name;
%       filestruct(i).rpdx1x2px_pxt_1{4} = dfile(4).name;
   end


   dfile = dir(['rpd_x1x2px_pxt_2_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdx1x2px_pxt_2{j} = dfile(j).name;
      end % (for j)
%       filestruct(i).rpdx1x2px_pxt_2{1} = dfile(1).name;
%       filestruct(i).rpdx1x2px_pxt_2{2} = dfile(2).name;
%       filestruct(i).rpdx1x2px_pxt_2{3} = dfile(3).name;
%       filestruct(i).rpdx1x2px_pxt_2{4} = dfile(4).name;
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





