function filtfiostr = mid_dat_filestr_to_filter_fio_struct(filestr)
% mid_dir_sm_file_struct_to_filters Filters/nonlinearities from file names
%
% filtfiostr = mid_dir_sm_file_struct_to_filters (filestr)
% --------------------------------------------------------------------------------
%
% Goes throught the struct array filestr and loads the filters and nonlinearities
% from the files stored in filestr.
%
% The data are returned in a struct array of the same length as filestr.
%
%


narginchk(1,1);

filepath = '.';

filtfiostr = [];

for i = 1:length(filestr)

    fprintf('Processing %.0f of %.0f\n', i, length(filestr));

        
    data = filestr(i);

    example_file = data.rpsta_files{1};

    index = findstr(example_file, '_');

    filter_specs = example_file(index(end-2)+1:index(end-1)-1);

    index = findstr(filter_specs, 'x');

    nf_filter = str2double(filter_specs(index(1)+1:index(2)-1));

    nt_filter = str2double(filter_specs(index(2)+1:end));


    % Assign intial parameters from MID analysis
    %--------------------------------------------------------------------
    data.nf_filter = nf_filter;
    data.nt_filter = nt_filter;

    
    %   STA Filter from MID Code
    %--------------------------------------------------------------------
    sta_files = fullfile(filepath, data.rpsta_files);

    [sta_filter, coeff_sta, projection_sta, filter_matrix_sta] = ...
        get_auditory_filter(sta_files, nf_filter, nt_filter);
    data.filter_mean_sta = sta_filter;
    data.coeff_sta = coeff_sta;
    data.projection_sta = projection_sta;
    data.filter_matrix_sta = filter_matrix_sta;

    

    %--------------------------------------------------------------------
    %   STA Nonlinearity from MID Code
    %--------------------------------------------------------------------
    sta_files = fullfile(filepath, data.rpx1pxpxt_sta_files);
    fio_sta = mid_get_dat_sta_fio(sta_files, coeff_sta);
    data.fio_sta = fio_sta;

    data.fio_sta
    pause
   


    % Get MID filters
    %--------------------------------------------------------------------
    filter_file_types = {...
        'rpdbest1_v1_files', ...
        'rpdbest1_v2_files', ...
        'rpdtest1_v1_files', ...
        'rpdtest1_v2_files', ...
        'rpdbest2_v1_files', ...
        'rpdbest2_v2_files', ...
        'rpdtest2_v1_files', ...
        'rpdtest2_v2_files'};

    filter_data_names = {...
        'best1_v1', ...
        'best1_v2', ...
        'test1_v1', ...
        'test1_v2', ...
        'best2_v1', ...
        'best2_v2', ...
        'test2_v1', ...
        'test2_v2'};
        

    for j = 1:length(filter_file_types)
        s0 = sprintf('length(data.%s)', filter_file_types{j});
        nfiles = eval(s0);

        if nfiles == 4

            total_files = {};
            for k = 1:4
                s1 = sprintf('data.%s{%.0f}', filter_file_types{j},k);
                s2 = sprintf('fullfile(filepath,%s);', s1);
                total_files{k} = eval(s2);
            end

            [filter_mean, coeff, projection, filter_matrix] = ...
                get_auditory_filter(total_files, nf_filter, nt_filter);

            string_data = sprintf('data.filter_mean_%s', filter_data_names{j});
            eval(sprintf('%s = filter_mean;', string_data));

            string_data = sprintf('data.coeff_%s', filter_data_names{j});
            eval(sprintf('%s = coeff;', string_data));

            string_data = sprintf('data.projection_%s', filter_data_names{j});
            eval(sprintf('%s = projection;', string_data));

            string_data = sprintf('data.filter_matrix_%s', filter_data_names{j});
            eval(sprintf('%s = filter_matrix;', string_data));

        else

            string_data = sprintf('data.filter_mean_%s', filter_data_names{j});
            eval(sprintf('%s = [];', string_data));

            string_data = sprintf('data.coeff_%s', filter_data_names{j});
            eval(sprintf('%s = [];', string_data));

            string_data = sprintf('data.projection_%s', filter_data_names{j});
            eval(sprintf('%s = [];', string_data));

            string_data = sprintf('data.filter_matrix_%s', filter_data_names{j});
            eval(sprintf('%s = [];', string_data));

        end

    end


    % Get MID nonlinearities
    %--------------------------------------------------------------------
    [fio_mid1] = get_dat_mid1_fio(data.rpdx1x2px_pxt_2_files, ...
        data.coeff_test2_v1, data.coeff_test2_v2);

    [fio_mid2] = get_dat_mid2_fio(data.rpdx1x2px_pxt_2_files, ...
        data.coeff_test2_v1, data.coeff_test2_v2);

    [fio_mid12] = get_dat_mid12_fio(data.rpdx1x2px_pxt_2_files, ...
        data.coeff_test2_v1, data.coeff_test2_v2);

    data.fio_mid1 = fio_mid1;
    data.fio_mid2 = fio_mid2;
    data.fio_mid12 = fio_mid12;

    filtfiostr = [filtfiostr data];
    clear('data');

end % (for i)


return;







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


for i = 1:length(filestruct)

	if ( isfield(filestruct, 'exp') )
	   filtstr(i).exp = filestruct(i).exp;
		filtstr(i).site = filestruct(i).site;
		filtstr(i).chan = filestruct(i).chan;
		filtstr(i).model = filestruct(i).model;
		filtstr(i).depth = filestruct(i).depth;
		filtstr(i).position = filestruct(i).position;
		filtstr(i).stim = filestruct(i).stim;
		filtstr(i).atten = filestruct(i).atten;
		filtstr(i).spl = filestruct(i).spl;
		filtstr(i).sm = filestruct(i).sm;
		filtstr(i).tm = filestruct(i).tm;
		filtstr(i).mdb = filestruct(i).mdb;
	else
	   filtstr(i).exp = [];
		filtstr(i).site = [];
		filtstr(i).chan = [];
		filtstr(i).model = [];
		filtstr(i).depth = [];
		filtstr(i).position = [];
		filtstr(i).stim = [];
		filtstr(i).atten = [];
		filtstr(i).spl = [];
		filtstr(i).sm = [];
		filtstr(i).tm = [];
		filtstr(i).mdb = [];
	end

   filtstr(i).location = filestruct(i).location;
   filtstr(i).unit = filestruct(i).unit;
   filtstr(i).x0 = filestruct(i).x0;
   filtstr(i).nh = filestruct(i).nh;
   filtstr(i).nv = filestruct(i).nv;
   filtstr(i).nlags = filestruct(i).nlags;
   filtstr(i).numtbins = numtbins;
   filtstr(i).numfbins = numfbins;
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







