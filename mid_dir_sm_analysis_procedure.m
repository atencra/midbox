function mid_dir_sm_analysis_procedure(varargin)

library('tmpmidbox');

options = struct('batch', 0);

options = input_options(options, varargin);


if ( options.batch )

    d = dir('.');
    isub = [d(:).isdir]; %# returns logical vector
    nameFolds = {d(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    outerfolder = pwd;

    for i = 1:length(nameFolds)

        cd(nameFolds{i});

        mid_dir_sm_file_names_to_file_struct('savepath', '.');

        mid_dir_sm_file_struct_to_filters;

        mid_dir_plot_filter_fio;

        mid_dir_plot_filter_fio2d;

        %mid_dir_pdf_plot_filter_fio;

        cd(outerfolder);

    end

else

    mid_dir_file_names_to_file_struct('savepath', '.');

    mid_dir_file_struct_to_filters;

    mid_dir_plot_filter_fio;

    mid_dir_plot_filter_fio2d;

    %mid_dir_pdf_plot_filter_fio;

end




return;
