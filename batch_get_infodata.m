function extrap_info = batch_get_infodata(sitenum)
% batch_get_infodata - Calculate and save extrapolated
% mutual information values for STA and MID filters.
%
% batch_get_infodata(sitenum)
% ========================================================
%
% This function must be run inside the directory
%
%       C:\MATLAB65\work\mid\
%
%
% sitenum : optional. If specified it tells which site(s)
%   to process. This is done because in some cases
%   different tuning curve parameters were used for
%   different sites. Should have the form [1 4 8 12].
%
%
% caa 9/16/07


if ( nargin == 0 )
   dsite = dir('site*');
else
   if ( isempty(sitenum) )
      dsite = dir('site*');
   else
      for i = 1:length(sitenum)
         dsite(i).name = ['site' num2str(sitenum(i))];
      end
   end
end

extrap_info = [];

   for i = 1:length(dsite)

      % get the file that contains the spike times and trigger
      site = dsite(i).name;
      dfile = dir([site '\'  'projection_information_data_*.mat']);

      if ( ~isempty(dfile) )

         infile = dfile.name

         outfile = [site '_extrapolation_information_data']

         load([site '\' infile]); % load in the spike times and the trigger

         if ( exist('proj')==1 & exist('information')==1 ) % make sure appropriate variables are in the file

            infodata = plot_proj_info(information);

            extrap_info = [extrap_info infodata];

            save([site '\' outfile], 'infodata');

         end % (if)

         clear proj info information infodata  % get rid of workspace variables

   end % (if)

end % (for i)


% save('extrapolation_information_database', 'extrap_info');