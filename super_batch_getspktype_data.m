function database = super_batch_get_strf_puretone_params_data(experiments)
% super_batch_get_strf_puretone_params_data - compute auto, cross correlations, spectrums 
%    between electrode channels in one big process
%
% database = super_batch_get_strf_puretone_params_data(experiments)
%
% experiments : a vector, optional. If specified it tells which experiment(s)
%   to process. 
%
% If the experiments input argument is not specified then the 
% following experiments will be processed:
%
%     20020619
%     20020730
%     20020826
%     20030305
%     20030408
%     20030506
%     20030527
%     20031028
%     20031112
%     20040114
%
% caa 1/11/06

if ( nargin==0 )
   experiments = [20020619; 20020730; 20020826; 20030305; 20030408; ...
                  20030506; 20030527; 20031028; 20031112; 20040114];
end

stim{1} = 'dmr1';
stim{2} = 'dmr2';

database = [];

for ii = 1:length(experiments)

   dsite = dir([num2str(experiments(ii)) '\site*']);

   for i = 1:length(dsite)

      expsite = [num2str(experiments(ii)) '\' dsite(i).name];

      % Process DMR1 spktype data
      % -----------------------------------------------
      dfile = dir([expsite '\*-dmr1-fs*-spktype.mat']);

      if ( ~isempty(dfile) )

         infile = dfile.name;
         load( [expsite '\' infile] );

         if ( exist('spktype')==1 )
            database = [database spktype];
         end % (if)

         clear spktype % get rid of workspace variables

      end % (if-elseif)


      % Process DMR2 spktype data
      % -----------------------------------------------
      dfile = dir([expsite '\*-dmr2-fs*-spktype.mat']);

      if ( ~isempty(dfile) )

         infile = dfile.name;
         load( [expsite '\' infile] );

         if ( exist('spktype')==1 )
            database = [database spktype];
         end % (if)

         clear spktype  % get rid of workspace variables

      end % (if-elseif)

   end % (for i)

end % (for ii)





