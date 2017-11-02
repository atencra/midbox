function [filtstr] = get_sgi_filters_from_dat_files(filestruct, prefix, paramfile)
% get_filters_from_dat_files - all filters from mid analysis in one data
%
% [filtstr] = get_filters_from_dat_files(filestruct, prefix, paramfile)
% -----------------------------------------------------------------------------
%
% filestruct : struct array holding files from the mid analysis. Has the 
%              following fields:
%
%     location
%     unit
%     x0
%     nh
%     nv
%     nlags
%     tbins
%     fbins
%     rpsta
%     rpx1pxpxt_sta
%     rpdbest2_v1
%     rpdbest2_v2
%     rpdtest2_v1
%     rpdtest2_v2
%     rpdx1x2px_pxt_2
%
%     location, unit, x0, nh, nb, nlags, tbins, and fbins are scalars.
%     The other fields are 1x4 cell arrays, with each element a string that
%     specifies the *.dat files for the corresponding data type.
%
%     In some cases *sta files do not exist. In this case the current function
%     calculates the STAs.
%
%     filestruct is obtained by:
%
%     filestruct = get_rippledir_detailed_mid_files(location, x0, nh, nlags)
%
% prefix : path to the folder containing the *.dat files specified in
% filestruct. If not included, then the default is 'dat_files\'
%
% paramfile : dmr parameter file. Needed because it holds the frequency and
% time vectors for the filters. If not included then the default is
%
% D:\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat
%
% caa 1/22/10



if ( nargin == 0 )
   error('You need between 1 and 4 input args.');
end

if ( nargin == 1 )
   prefix = 'dat_files\';
   paramfile = 'D:\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat';
end

if ( nargin == 2 )
   paramfile = 'D:\stimuli\20031124\dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8_param.mat';
end




nh = filestruct(1).nh;
nv = filestruct(1).nv;
nlags = filestruct(1).nlags;
numtbins = filestruct(1).tbins; %20;
numfbins = filestruct(1).fbins; %25;
x0 = filestruct(1).x0; %20;
index_freq = (x0):(numfbins-1+x0);

s = load(paramfile,'taxis', 'faxis');
taxis = s.taxis;
time = taxis(1:numtbins); % time axis for filters
faxis = s.faxis;
freq = faxis(index_freq); % frequency axis for filters



% if ( length(locator) ~= size(stimulus,2) )
%    error('Spike train and envelope file have different number of trials.');
% end

filtstr = filestruct;

for i = 1:length(filestruct)

% 	if ( isfield(filestruct, 'exp') )
% 	   filtstr(i).exp = filestruct(i).exp;
% 		filtstr(i).site = filestruct(i).site;
% 		filtstr(i).chan = filestruct(i).chan;
% 		filtstr(i).model = filestruct(i).model;
% 		filtstr(i).depth = filestruct(i).depth;
% 		filtstr(i).position = filestruct(i).position;
% 		filtstr(i).stim = filestruct(i).stim;
% 		filtstr(i).atten = filestruct(i).atten;
% 		filtstr(i).spl = filestruct(i).spl;
% 		filtstr(i).sm = filestruct(i).sm;
% 		filtstr(i).tm = filestruct(i).tm;
% 		filtstr(i).mdb = filestruct(i).mdb;
% 	else
% 	   filtstr(i).exp = [];
% 		filtstr(i).site = [];
% 		filtstr(i).chan = [];
% 		filtstr(i).model = [];
% 		filtstr(i).depth = [];
% 		filtstr(i).position = [];
% 		filtstr(i).stim = [];
% 		filtstr(i).atten = [];
% 		filtstr(i).spl = [];
% 		filtstr(i).sm = [];
% 		filtstr(i).tm = [];
% 		filtstr(i).mdb = [];
% 	end

%    filtstr(i).location = filestruct(i).location;
%    filtstr(i).unit = filestruct(i).unit;
%    filtstr(i).x0 = filestruct(i).x0;
%    filtstr(i).nh = filestruct(i).nh;
%    filtstr(i).nv = filestruct(i).nv;
%    filtstr(i).nlags = filestruct(i).nlags;
%    filtstr(i).numtbins = numtbins;
%    filtstr(i).numfbins = numfbins;

   filtstr(i).time = time;
   filtstr(i).freq = freq;


   %--------------------------------------------------------------------
   %   STA Filter from MID Code
   %--------------------------------------------------------------------

   if ( ~isempty(filestruct(i).rpsta) && length(filestruct(i).rpsta)==4 )

      if ( isempty(prefix) )
         file_sta = filestruct(i).rpsta;
      else
         for j = 1:length( filestruct(i).rpdtest2_v1 )
            file_sta{j} = sprintf('%s%s', prefix, filestruct(i).rpsta{j});
         end
      end

      [v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_filter(file_sta, nh, nv, nlags);
   
      filtstr(i).v_sta = v_sta;
      filtstr(i).mtx_sta = mtx_sta;

   elseif ( ~isempty(filestruct(i).rpsta) && length(filestruct(i).rpsta)~=4 )
      filtstr(i).v_sta = [];
      filtstr(i).mtx_sta = [];
   else
      error('Code not ready to compute STAs just yet.');
   end





   %--------------------------------------------------------------------
   %   MID1 Filter from MID Code
   %--------------------------------------------------------------------

   if ( ~isempty(filestruct(i).rpdtest2_v1) && length(filestruct(i).rpdtest2_v1)==4 )

      if ( isempty(prefix) )
         file_v1 = filestruct(i).rpdtest2_v1;
      else
         for j = 1:length( filestruct(i).rpdtest2_v1 )
            file_v1{j} = sprintf('%s%s', prefix, filestruct(i).rpdtest2_v1{j});
         end
      end

      [v1, coeff_v1, projection_v1, mtx_v1] = get_auditory_filter(file_v1, nh, nv, nlags);
   
      filtstr(i).v1 = v1;
      filtstr(i).mtx_v1 = mtx_v1;

   elseif ( ~isempty(filestruct(i).rpdtest2_v1) && length(filestruct(i).rpdtest2_v1)~=4 )
      filtstr(i).v1 = [];
      filtstr(i).mtx_v1 = [];
   else
      fprintf('You need to run the mid code to compute the second filter.\n');
      %error('You need to run the mid code to compute the second filter.');
   end



   %--------------------------------------------------------------------
   %   MID2 Filter from MID Code
   %--------------------------------------------------------------------

   if ( ~isempty(filestruct(i).rpdtest2_v2) && length(filestruct(i).rpdtest2_v2)==4 )

      if ( isempty(prefix) )
         file_v2 = filestruct(i).rpdtest2_v2;
      else
         for j = 1:length( filestruct(i).rpdtest2_v2 )
            file_v2{j} = sprintf('%s%s', prefix, filestruct(i).rpdtest2_v2{j});
         end
      end

      [v2, coeff_v2, projection_v2, mtx_v2] = get_auditory_filter(file_v2, nh, nv, nlags);
   
      filtstr(i).v2 = v2;
      filtstr(i).mtx_v2 = mtx_v2;

   elseif ( ~isempty(filestruct(i).rpdtest2_v2) && length(filestruct(i).rpdtest2_v2)~=4 )
      filtstr(i).v2 = [];
      filtstr(i).mtx_v2 = [];
   else
      fprintf('You need to run the mid code to compute the second filter.\n');
      %error('You need to run the mid code to compute the second filter.');
   end

end

return;







