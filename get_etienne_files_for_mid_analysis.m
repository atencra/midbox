function filestruct = get_etienne_files_for_mid_analysis(expdate)
% get_files_for_mid_analysis - put MID data files into a struct array
%    so data can be analysed
%
% filestruct = get_files_for_mid_analysis(expdate)
%
% expdate : Either 20121031 or 20121220.
%
%
% caa 4/11/06

if ( nargin~=1 )
   error('You need one input arg.');
end

current_directory = pwd;


dfile = dir(['rpsta_' num2str(expdate) '*.dat']);

length(dfile)
return;

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
'exp',            [], ...
'unit',           [], ...
'tbins',          [], ...
'fbins',          [], ...
'rpsta',          [], ...
'rpx1pxpxt_sta',  [], ...
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

%{
for i = 1:length(unitnum)

   filestruct(i).exp = expdate;
   filestruct(i).unit = unitnum(i);
   filestruct(i).tbins = tbins;
   filestruct(i).fbins = fbins;


   % The RPSTA Files
   %=================================================================

   % Get the STA files
   dfile = dir(['rpsta_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) ) 
      for j = 1:length(dfile)
         filestruct(i).rpsta{j} = dfile(j).name;
      end % (for j)
   end


   % Get the STA input/output function files
   dfile = dir(['rpx1pxpxt_sta_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpx1pxpxt_sta{j} = dfile(j).name;
      end % (for j)
   end


   % The RPDTEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v1_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdtest1_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest1_v2_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdtest1_v2{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v1_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdtest2_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdtest2_v2_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdtest2_v2{j} = dfile(j).name;
      end % (for j)
   end



   % The RPDBEST Files
   %=================================================================

   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v1_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdbest1_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest1_v2_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdbest1_v2{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID1 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v1_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdbest2_v1{j} = dfile(j).name;
      end % (for j)
   end


   % Get the MID2 files - First Most Informative Dimension
   dfile = dir(['rpdbest2_v2_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdbest2_v2{j} = dfile(j).name;
      end % (for j)
   end



   % Get the MID input/output function files - First and Second Most
   % Informative Dimensions
   dfile = dir(['rpdx1x2px_pxt_1_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdx1x2px_pxt_1{j} = dfile(j).name;
      end % (for j)
   end


   dfile = dir(['rpdx1x2px_pxt_2_' num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);
   if ( ~isempty(dfile) )
      for j = 1:length(dfile)
         filestruct(i).rpdx1x2px_pxt_2{j} = dfile(j).name;
      end % (for j)
   end

end % (for ii)


return;

%}

filetype = {'rpsta_', ...
'rpx1pxpxt_sta_', ...
'rpdtest1_v1_', ...
'rpdtest1_v2_', ...
'rpdtest2_v1_', ...
'rpdtest2_v2_', ...
'rpdbest1_v1_', ...
'rpdbest1_v2_', ...
'rpdbest2_v1_', ...
'rpdbest2_v2_', ...
'rpdx1x2px_pxt_1_', ...
'rpdx1x2px_pxt_2_'}


savetype = {'filestruct(i).rpsta{j} = dfile(j).name;', ...
'filestruct(i).rpx1pxpxt_sta{j} = dfile(j).name;', ...
'filestruct(i).rpdtest1_v1{j} = dfile(j).name;', ...
'filestruct(i).rpdtest1_v2{j} = dfile(j).name;', ...
'filestruct(i).rpdtest2_v1{j} = dfile(j).name;', ...
'filestruct(i).rpdtest2_v2{j} = dfile(j).name;', ...
'filestruct(i).rpdbest1_v1{j} = dfile(j).name;', ...
'filestruct(i).rpdbest1_v2{j} = dfile(j).name;', ...
'filestruct(i).rpdbest2_v1{j} = dfile(j).name;', ...
'filestruct(i).rpdbest2_v2{j} = dfile(j).name;', ...
'filestruct(i).rpdx1x2px_pxt_1{j} = dfile(j).name;', ...
'filestruct(i).rpdx1x2px_pxt_2{j} = dfile(j).name;'}


for i = 1:length(unitnum)

   filestruct(i).exp = expdate;
   filestruct(i).unit = unitnum(i);
   filestruct(i).tbins = tbins;
   filestruct(i).fbins = fbins;


   for ii = 1:length(filetype)

      % Get the files
      dfile = dir([filetype{ii} num2str(expdate) '_' num2str(unitnum(i)) '_*.dat']);

      if ( ~isempty(dfile) ) % save the files to the struct
         for j = 1:length(dfile)
            eval(savetype{ii})
         end % (for j)
      end

   end % (for ii)

end % (for i)


return;

















