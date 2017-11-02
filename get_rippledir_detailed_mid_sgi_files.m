function filestruct = get_rippledir_detailed_mid_sgi_files(location, strf)
% filestruct = get_rippledir_detailed_mid_sgi_files(location, x0, nh,
% nlags, strf)
% get_files_for_mid_icc_analysis - put MID data files into a struct array
%    so data can be analysed
%
% filestruct = get_rippledir_detailed_mid_sgi_files(location, x0, nh, nlags, strf)
%
% Get mid files from the sgi cluster. Since the cluster reads in a
% parameter file, we will read through this file to get the parameters
% that were used, and then we will save the data to a struct array.
%
% This function requires that a ripple_params_*_.txt file is present in
% the directory with the .dat files.
%
% This differs from earlier versions of the code, since older versions
% didn't use an input parameter file.
%
% location, for auditory cortex : a scalar. It is the penetration number
% for the experiment. Usually something like 32, 34, or 38.
%
% x0 : starting frequency bin for mid analysis.
%
% nh = 25 and nlags = 20, in all cases.
%
% strf : strf data from which the original mid data were originally
% extracted. The function will still run if strf is not supplied. If strf
% is supplied, then the parameters within in will be added to the output
% struct array. These parameters are TM, SM MdB, experiment, site, chan,
% model, etc. for the neuron.
%
% caa 1/28/10

error(nargchk(1,2,2));

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
'expdate',         [], ...
'location',        [], ...
'unit',            [], ...
'dimx',            [], ...
'x0',              [], ...
'ntrials',         [], ...
'maxcount1',       [], ...
'maxcount2',       [], ...
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

nv = 1;
nh = 25;
nlags = 20;
numtbins = nlags;
numfbins = nh;

txtfile = dir(['ripple_params_*_*_' num2str(location) '_*_1.txt']);

for i = 1:length(txtfile)

   fid = fopen(txtfile(i).name, 'r');

   stimfile = fgets(fid);
   iskfile = fgets(fid);
   expdate = str2double(fgets(fid));
   site = str2double(fgets(fid));
   unit = str2double(fgets(fid));
   dimx = str2double(fgets(fid));
   x0 = str2double(fgets(fid));
   ntrials = str2double(fgets(fid));
   maxcount1 = str2double(fgets(fid));
   maxcount2 = str2double(fgets(fid));

   fclose all;

%    unit = unitnum(i);

   if ( nargin == 2 )
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
   end

   filestruct(i).expdate = expdate;
   filestruct(i).location = location;
   filestruct(i).unit = unit;
   filestruct(i).dimx = dimx;
   filestruct(i).x0 = x0;
   filestruct(i).ntrials = ntrials;
   filestruct(i).maxcount1 = maxcount1;
   filestruct(i).maxcount2 = maxcount2;

   filestruct(i).nh = nh;
   filestruct(i).nv = nv;
   filestruct(i).nlags = nlags;
   filestruct(i).numtbins = numtbins;
   filestruct(i).numfbins = numfbins;


   % The RPSTA Files
   %=================================================================

   % Get the STA files
   dfile = dir(['rpsta_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
%    dfile = dir(['rpsta_' num2str(location) '_' num2str(unitnum(i)) '_*.dat']);
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
   dfile = dir(['rpx1pxpxt_sta_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpx1pxpxt_sta{j} = dfile(j).name;
      end % (for j)
   end


   % The RPDTEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v1_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest1_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v2_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest1_v2{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v1_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest2_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v2_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdtest2_v2{j} = dfile(j).name;
      end % (for j)
   end



   % The RPDBEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v1_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest1_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v2_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest1_v2{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v1_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest2_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v2_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdbest2_v2{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID input/output function files - First FIO
   dfile = dir(['rpdx1x2px_pxt_1_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdx1x2px_pxt_1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID input/output function files - First and Second Most FIO
   dfile = dir(['rpd_x1x2px_pxt_2_' num2str(expdate) '_' num2str(location) '_' num2str(unit) '_*.dat']);
   if ( ~isempty(dfile) )
%    if ( ~isempty(dfile) & length(dfile) == 4 )
      for j = 1:length(dfile)
         filestruct(i).rpdx1x2px_pxt_2{j} = dfile(j).name;
      end % (for j)
   end

end % (for ii)


return;







