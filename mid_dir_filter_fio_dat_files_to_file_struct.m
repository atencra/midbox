function filestr = mid_dir_filter_fio_dat_files_to_file_struct(varargin)
%  mid_dir_filter_fio_dat_files_to_file_struct - Collect MID dat files 
%
%  filestr = mid_dir_filter_fio_dat_files_to_file_struct
%  ===============================================================
%  Searches through the current directory and finds the file names
%  for each single unit that was analyzed, and saves the file names
%  to a struct array.
%
%  Example filestr struct:
%
%                            site: '602'
%                            unit: '1'
%                     rpsta_files: {1x4 cell}
%             rpx1pxpxt_sta_files: {1x4 cell}
%               rpdbest1_v1_files: {1x4 cell}
%               rpdbest1_v2_files: {1x4 cell}
%               rpdtest1_v1_files: {1x4 cell}
%               rpdtest1_v2_files: {1x4 cell}
%               rpdbest2_v1_files: {1x4 cell}
%               rpdbest2_v2_files: {1x4 cell}
%               rpdtest2_v1_files: {1x4 cell}
%               rpdtest2_v2_files: {1x4 cell}
%           rpdx1x2px_pxt_1_files: {1x4 cell}
%          rpd_x1x2px_pxt_2_files: {1x4 cell}
%
%

options = struct('savepath', []);

options = input_options(options, varargin);



d = dir('rpdtest2_v1_*.dat');

if ( isempty(d) )
    error('No test2 v1 files in current directory.');
end

site_unit = {};
for i = 1:length(d)
    filename = d(i).name;
    index = strfind(filename, '_');
    site_unit{i} = filename( index(end-4)+1 : index(end-2)-1 );
end

site_unit_list = unique(site_unit);
clear('site_unit');


filestr = [];

for i = 1:length(site_unit_list)

    site_unit = site_unit_list{i};

    parts = strsplit(site_unit,'_');

    site = parts{1};
    unit = parts{2};

    df.site = site;
    df.unit = unit;


    % STA filter file names
    %============================================================
    rpsta_files = dir(sprintf('rpsta_*%s_*.dat', site_unit));

    if ( length(rpsta_files) ~= 4 )
        warning(sprintf('length(rpsta_files) = %.0f', ...
            length(rpsta_files)));
    end
    rpsta_files = {rpsta_files.name};
    df.rpsta_files = rpsta_files;



    % STA nonlinearity file names
    %============================================================
    rpx1pxpxt_sta_files = dir(sprintf('rpx1pxpxt_sta_*%s_*.dat', site_unit));

    if ( length(rpx1pxpxt_sta_files) ~= 4 )
        warning(sprintf('length(rpx1pxpxt_sta_files) = %.0f', ...
            length(rpx1pxpxt_sta_files)));
    end
    rpx1pxpxt_sta_files = {rpx1pxpxt_sta_files.name};
    df.rpx1pxpxt_sta_files = rpx1pxpxt_sta_files;



    % Best1 MID files
    %============================================================
    rpdbest1_v1_files = dir(sprintf('rpdbest1_v1_*%s_*.dat', site_unit));

    if ( length(rpdbest1_v1_files) ~= 4 )
        warning(sprintf('length(rpdbest1_v1_files) = %.0f', ...
            length(rpdbest1_v1_files)));
    end
    rpdbest1_v1_files = {rpdbest1_v1_files.name};
    df.rpdbest1_v1_files = rpdbest1_v1_files;


    rpdbest1_v2_files = dir(sprintf('rpdbest1_v2_*%s_*.dat', site_unit));

    if ( length(rpdbest1_v2_files) ~= 4 )
        warning(sprintf('length(rpdbest1_v2_files) = %.0f', ...
            length(rpdbest1_v2_files)));
    end
    rpdbest1_v2_files = {rpdbest1_v2_files.name};
    df.rpdbest1_v2_files = rpdbest1_v2_files;



    % Test1 MID files
    %============================================================
    rpdtest1_v1_files = dir(sprintf('rpdtest1_v1_*%s_*.dat', site_unit));

    if ( length(rpdtest1_v1_files) ~= 4 )
        warning(sprintf('length(rpdtest1_v1_files) = %.0f', ...
            length(rpdtest1_v1_files)));
    end
    rpdtest1_v1_files = {rpdtest1_v1_files.name};
    df.rpdtest1_v1_files = rpdtest1_v1_files;



    rpdtest1_v2_files = dir(sprintf('rpdtest1_v2_*%s_*.dat', site_unit));

    if ( length(rpdtest1_v2_files) ~= 4 )
        warning(sprintf('length(rpdtest1_v2_files) = %.0f', ...
            length(rpdtest1_v2_files)));
    end
    rpdtest1_v2_files = {rpdtest1_v2_files.name};
    df.rpdtest1_v2_files = rpdtest1_v2_files;



    % Best2 MID files
    %============================================================
    rpdbest2_v1_files = dir(sprintf('rpdbest2_v1_*%s_*.dat', site_unit));

    if ( length(rpdbest2_v1_files) ~= 4 )
        warning(sprintf('length(rpdbest2_v1_files) = %.0f', ...
            length(rpdbest2_v1_files)));
    end
    rpdbest2_v1_files = {rpdbest2_v1_files.name};
    df.rpdbest2_v1_files = rpdbest2_v1_files;



    rpdbest2_v2_files = dir(sprintf('rpdbest2_v2_*%s_*.dat', site_unit));

    if ( length(rpdbest2_v2_files) ~= 4 )
        warning(sprintf('length(rpdbest2_v2_files) = %.0f', ...
            length(rpdbest2_v2_files)));
    end
    rpdbest2_v2_files = {rpdbest2_v2_files.name};
    df.rpdbest2_v2_files = rpdbest2_v2_files;


    % Test2 MID files
    %============================================================
    rpdtest2_v1_files = dir(sprintf('rpdtest2_v1_*%s_*.dat', site_unit));

    if ( length(rpdtest2_v1_files) ~= 4 )
        warning(sprintf('length(rpdtest2_v1_files) = %.0f', ...
            length(rpdtest2_v1_files)));
    end
    rpdtest2_v1_files = {rpdtest2_v1_files.name};
    df.rpdtest2_v1_files = rpdtest2_v1_files;


    rpdtest2_v2_files = dir(sprintf('rpdtest2_v2_*%s_*.dat', site_unit));

    if ( length(rpdtest2_v2_files) ~= 4 )
        warning(sprintf('length(rpdtest2_v2_files) = %.0f', ...
            length(rpdtest2_v2_files)));
    end
    rpdtest2_v2_files = {rpdtest2_v2_files.name};
    df.rpdtest2_v2_files = rpdtest2_v2_files;



    % MID optimization #1 nonlinearity files
    %============================================================
    rpdx1x2px_pxt_1_files = dir(sprintf('rpdx1x2px_pxt_1_*%s_*.dat', site_unit));

    if ( length(rpdx1x2px_pxt_1_files) ~= 4 )
        warning(sprintf('length(rpdx1x2px_pxt_1_files) = %.0f', ...
            length(rpdx1x2px_pxt_1_files)));
    end
    rpdx1x2px_pxt_1_files  = {rpdx1x2px_pxt_1_files.name};
    df.rpdx1x2px_pxt_1_files = rpdx1x2px_pxt_1_files;



    % MID optimization #2 nonlinearity files
    %============================================================
    rpdx1x2px_pxt_2_files = dir(sprintf('rpd_x1x2px_pxt_2_*%s_*.dat', site_unit));

    if ( length(rpdx1x2px_pxt_2_files) ~= 4 )
        warning(sprintf('length(rpd_x1x2px_pxt_2_files) = %.0f', ...
            length(rpdx1x2px_pxt_2_files)));
    end
    rpdx1x2px_pxt_2_files = {rpdx1x2px_pxt_2_files.name};
    df.rpdx1x2px_pxt_2_files = rpdx1x2px_pxt_2_files;


    filestr = [filestr df];

    if ~isempty(options.savepath)
        outfile = sprintf('%s-filestr.mat', site_unit);
        outfile = fullfile(options.savepath, outfile);
        save(outfile, 'filestr');
        filestr = [];
    end

    clear('df');

end


return;




