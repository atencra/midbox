function mid_dir_plot_train_test_filters_nonlinearities(batch)
%mid_dir_plot_train_test_filters_nonlinearities
%
%    mid_dir_plot_train_test_filters_nonlinearities(batch)
% 
%    mid_dir_pdf_plot_train_test_filters_nonlinearities(batch) goes through 
%    the folder datapath and searches for *-filter-fio-proj-info.mat files.
%    These files hold filter, nonlinearity, and information results for
%    the MID analysis.
%
%    It then plots the filters and nonlineaities within the folder,
%    converts the figures to pdf files, and saves the pdfs in the
%    folder figpath. The individual pdfs are appended, and then they
%    are deleted. Thus, there is one multiple figure pdf for each
%    neuron.
%
%    batch : if 1, then the function assumes it is outside the folders
%       holding the data. It then goes in to each folder, makes plots, and 
%       comes back out. Default = 0. 
% 
%    This function requires that the Ghostscript executable be installed.
%    It is assumed to have the name and path:
%       gsfile = 'gswin64c.exe';
%       gsdir = 'C:\Program Files\gs\gs9.20\bin';
%    


library('export_fig');


narginchk(0,1);

if nargin == 0
    batch = 0;
end

if isempty(batch)
    batch = 0;
end


figpath = '.';
datapath = '.';
gsdir = 'C:\Program Files\gs\gs9.20\bin';
gsfile = 'gswin64c.exe';
gspath = fullfile(gsdir, gsfile);



if ( batch )
    folders = get_directory_names;
else
    folders = {'.'};
end

startdir = pwd;

for ii = 1:length(folders)

    cd(folders{ii});

    dfilters = dir(fullfile(datapath, '*-filter-fio-proj-info.mat'));
    if isempty(dfilters)
        fprintf('\nNo *-filter-fio-proj-info.mat files in %s.\n\n', folders{ii});
    end

    matfiles = {dfilters.name};


    for i = 1:length(matfiles)
       
        fprintf('\nPlotting %s\n', matfiles{i});
        load(matfiles{i}, 'fio', 'projinfo');

        eps_file_filter = strrep(matfiles{i}, '.mat', '-filter.eps');
        pdf_file_filter = strrep(matfiles{i}, '.mat', '-filter.pdf');

        eps_file_fio = strrep(matfiles{i}, '.mat', '-fio.eps');
        pdf_file_fio = strrep(matfiles{i}, '.mat', '-fio.pdf');

        eps_file_info0 = strrep(matfiles{i}, '.mat', '-info0.eps');
        pdf_file_info0 = strrep(matfiles{i}, '.mat', '-info0.pdf');

        eps_file_info1 = strrep(matfiles{i}, '.mat', '-info1.eps');
        pdf_file_info1 = strrep(matfiles{i}, '.mat', '-info1.pdf');

        eps_file_info2 = strrep(matfiles{i}, '.mat', '-info2.eps');
        pdf_file_info2 = strrep(matfiles{i}, '.mat', '-info2.pdf');

        eps_file_info12 = strrep(matfiles{i}, '.mat', '-info12.eps');
        pdf_file_info12 = strrep(matfiles{i}, '.mat', '-info12.pdf');

        iskfile = strrep(fio.iskfile, '_', '-');

        crop = 1;
        append = 0;
        gray = 0;
        quality = 1000;


        mid_plot_train_fio_1d_2d(fio);
        fig2eps(eps_file_fio);
        pause(0.5);
        %eps2pdf1(eps_file, gspath);
        eps2pdf(eps_file_fio, pdf_file_fio, crop, append, gray, quality);
        pause(0.5);
        close all;
        pause(1)


        mid_plot_train_fio_filters(fio);
        fig2eps(eps_file_filter);
        pause(0.5);
        %eps2pdf1(eps_file, gspath);
        eps2pdf(eps_file_filter, pdf_file_filter, crop, append, gray, quality);
        pause(0.5);
        close all;
        pause(1)



        fraction = [90 92.5 95 97.5 100];

        mid_plot_projinfo_train_test_information(fraction, projinfo.ifrac0_mtx_train, ...
           projinfo.ifrac0_mn_train, projinfo.ifrac0_std_train, ...
           projinfo.ifrac0_mtx_test, projinfo.ifrac0_mn_test, projinfo.ifrac0_std_test, ...
           'STA Info Analysis');
        fig2eps(eps_file_info0);
        pause(0.5);
        eps2pdf(eps_file_info0, pdf_file_info0, crop, append, gray, quality);
        pause(0.5);
        close all;
        pause(1)



        mid_plot_projinfo_train_test_information(fraction, projinfo.ifrac1_mtx_train, ...
           projinfo.ifrac1_mn_train, projinfo.ifrac1_std_train, ...
           projinfo.ifrac1_mtx_test, projinfo.ifrac1_mn_test, projinfo.ifrac1_std_test, ...
           'MID1 Info Analysis');
        fig2eps(eps_file_info1);
        pause(0.5);
        eps2pdf(eps_file_info1, pdf_file_info1, crop, append, gray, quality);
        pause(0.5);
        close all;
        pause(1)



        mid_plot_projinfo_train_test_information(fraction, projinfo.ifrac2_mtx_train, ...
           projinfo.ifrac2_mn_train, projinfo.ifrac2_std_train, ...
           projinfo.ifrac2_mtx_test, projinfo.ifrac2_mn_test, projinfo.ifrac2_std_test, ...
           'MID2 Info Analysis');
        fig2eps(eps_file_info2);
        pause(0.5);
        eps2pdf(eps_file_info2, pdf_file_info2, crop, append, gray, quality);
        pause(0.5);
        close all;
        pause(1)



        mid_plot_projinfo_train_test_information(fraction, projinfo.ifrac12_mtx_train, ...
           projinfo.ifrac12_mn_train, projinfo.ifrac12_std_train, ...
           projinfo.ifrac12_mtx_test, projinfo.ifrac12_mn_test, projinfo.ifrac12_std_test, ...
           'MID12 Info Analysis');
        fig2eps(eps_file_info12);
        pause(0.5);
        eps2pdf(eps_file_info12, pdf_file_info12, crop, append, gray, quality);
        pause(0.5);
        close all;
        pause(1)



        output_pdf = strrep(matfiles{i}, '.mat', '.pdf');
        delete(output_pdf);

        append_pdfs(output_pdf, ...
                    pdf_file_filter, ...
                    pdf_file_fio, ...
                    pdf_file_info0, ...
                    pdf_file_info1, ...
                    pdf_file_info2, ...
                    pdf_file_info12 ); 

        fprintf('\nCreated %s\n\n', output_pdf);

        delete(eps_file_filter);
        delete(pdf_file_filter);

        delete(eps_file_fio);
        delete(pdf_file_fio);

        delete(eps_file_info0);
        delete(pdf_file_info0);

        delete(eps_file_info1);
        delete(pdf_file_info1);

        delete(eps_file_info2);
        delete(pdf_file_info2);

        delete(eps_file_info12);
        delete(pdf_file_info12);

    end % (for i)

    cd(startdir)

end % (for ii)


return;






function mid_plot_train_fio_filters(fio)

numfbins = fio.nf_filter;
numtbins = fio.nt_filter;

sta{1} = reshape( fio.filter_matrix_sta(:,1), numfbins, numtbins );
sta{2} = reshape( fio.filter_matrix_sta(:,2), numfbins, numtbins );
sta{3} = reshape( fio.filter_matrix_sta(:,3), numfbins, numtbins );
sta{4} = reshape( fio.filter_matrix_sta(:,4), numfbins, numtbins );

mid1{1} = reshape( fio.filter_matrix_test2_v1(:,1), numfbins, numtbins );
mid1{2} = reshape( fio.filter_matrix_test2_v1(:,2), numfbins, numtbins );
mid1{3} = reshape( fio.filter_matrix_test2_v1(:,3), numfbins, numtbins );
mid1{4} = reshape( fio.filter_matrix_test2_v1(:,4), numfbins, numtbins );

mid2{1} = reshape( fio.filter_matrix_test2_v2(:,1), numfbins, numtbins );
mid2{2} = reshape( fio.filter_matrix_test2_v2(:,2), numfbins, numtbins );
mid2{3} = reshape( fio.filter_matrix_test2_v2(:,3), numfbins, numtbins );
mid2{4} = reshape( fio.filter_matrix_test2_v2(:,4), numfbins, numtbins );


% Plot the filters
% ----------------------------------------------------------

figure;

for i = 1:4
    subplot(3,4,i);
    imagesc( sta{i} );
    minmin = min(min(sta{i}));
    maxmax = max(max(sta{i}));
    boundary = max([abs(minmin) abs(maxmax)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    cmap = brewmaps('rdbu',21);
    colormap(cmap);
    if i == 1
        ylabel('STA');
    end
    title(sprintf('Train Set %.0f',i));
end % (for i)




for i = 1:4
    subplot(3,4,i+4);
    imagesc( mid1{i} );
    minmin = min(min(mid1{i}));
    maxmax = max(max(mid1{i}));
    boundary = max([abs(minmin) abs(maxmax)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    cmap = brewmaps('rdbu',21);
    colormap(cmap);
    if i == 1
        ylabel('MID1');
    end
end % (for i)


for i = 1:4
    subplot(3,4,i+8);
    imagesc( mid2{i} );
    minmin = min(min(mid2{i}));
    maxmax = max(max(mid2{i}));
    boundary = max([abs(minmin) abs(maxmax)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    cmap = brewmaps('rdbu',21);
    colormap(cmap);
    if i == 1
        ylabel('MID2');
    end
end % (for i)

return;








function mid_plot_train_fio_1d_2d(fio)

iskfile = fio.iskfile;
nspk = sum(fio.locator);
ms = 3;

figure;


% Plot the STA nonlinearities
% ----------------------------------------------------------

xbins = fio.x0bins;
pspk = fio.pspk;
pspkx = fio.pspkx0;
        
maxmax = max([max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

for i = 1:4

    subplot(4,4,i);
    hold on;
    plot(xbins{i}, pspkx{i}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
    plot([min(xbins{i}) max(xbins{i})], [pspk{i} pspk{i}], 'k--');
    xlim(1.1*[min(xbins{i}) max(xbins{i})]);
    ylim([0 maxmax]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

    if i == 1
        ylabel('STA');
    end

end % (for i)


%subplot(4,4,2);
%hold on;
%plot(xbins{2}, pspkx{2}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{2}) max(xbins{2})], [pspk{2} pspk{2}], 'k--');
%xlim(1.1*[min(xbins{2}) max(xbins{2})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%
%
%subplot(4,4,3);
%hold on;
%plot(xbins{3}, pspkx{3}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{3}) max(xbins{3})], [pspk{3} pspk{3}], 'k--');
%xlim(1.1*[min(xbins{3}) max(xbins{3})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%
%
%subplot(4,4,4);
%hold on;
%plot(xbins{4}, pspkx{4}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{4}) max(xbins{4})], [pspk{4} pspk{4}], 'k--');
%xlim(1.1*[min(xbins{4}) max(xbins{4})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);




% Plot the MID1 nonlinearities
% ----------------------------------------------------------

xbins = fio.x1bins;
pspkx = fio.pspkx1;
        
maxmax = max([max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);

for i = 1:4
    subplot(4,4,i+4);
    hold on;
    plot(xbins{i}, pspkx{i}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
    plot([min(xbins{i}) max(xbins{i})], [pspk{1} pspk{1}], 'k--');
    xlim(1.1*[min(xbins{i}) max(xbins{i})]);
    ylim([0 maxmax]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    if i == 1
        ylabel('MID1');
    end
end % (for i)

%subplot(4,4,6);
%hold on;
%plot(xbins{2}, pspkx{2}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{2}) max(xbins{2})], [pspk{2} pspk{2}], 'k--');
%xlim(1.1*[min(xbins{2}) max(xbins{2})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%
%subplot(4,4,7);
%hold on;
%plot(xbins{3}, pspkx{3}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{3}) max(xbins{3})], [pspk{3} pspk{3}], 'k--');
%xlim(1.1*[min(xbins{3}) max(xbins{3})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%
%subplot(4,4,8);
%hold on;
%plot(xbins{4}, pspkx{4}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{4}) max(xbins{4})], [pspk{4} pspk{4}], 'k--');
%xlim(1.1*[min(xbins{4}) max(xbins{4})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

%if (  ~isempty(titlestr) )
%   suptitle(titlestr);
%end
%
%set(gcf, 'position', [293 267 1095 571]);




% Plot the nonlinearities
% ----------------------------------------------------------

xbins = fio.x2bins;
pspkx = fio.pspkx2;
        
maxmax = max([max(pspkx{1}) max(pspkx{2}) max(pspkx{3}) max(pspkx{4}) ]);
for i = 1:4
    subplot(4,4,i+8);
    hold on;
    plot(xbins{i}, pspkx{i}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
    plot([min(xbins{i}) max(xbins{i})], [pspk{1} pspk{1}], 'k--');
    xlim(1.1*[min(xbins{i}) max(xbins{i})]);
    ylim([0 maxmax]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    if i == 1
        ylabel('MID2');
    end
end % (for i)


%subplot(4,4,10);
%hold on;
%plot(xbins{2}, pspkx{2}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{2}) max(xbins{2})], [pspk{2} pspk{2}], 'k--');
%xlim(1.1*[min(xbins{2}) max(xbins{2})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%
%subplot(4,4,11);
%hold on;
%plot(xbins{3}, pspkx{3}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{3}) max(xbins{3})], [pspk{3} pspk{3}], 'k--');
%xlim(1.1*[min(xbins{3}) max(xbins{3})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
%
%subplot(4,4,12);
%hold on;
%plot(xbins{4}, pspkx{4}, 'ko-', 'markerfacecolor', 'k', 'markersize', ms);
%plot([min(xbins{4}) max(xbins{4})], [pspk{4} pspk{4}], 'k--');
%xlim(1.1*[min(xbins{4}) max(xbins{4})]);
%ylim([0 maxmax]);
%set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);



% 2D nonlinearity
%------------------------------


x1binedges = fio.x1binedges;
x2binedges = fio.x2binedges;
pspkx1x2 = fio.pspkx1x2;

for i = 1:4
    subplot(4,4,i+12);
    mxmx = max(pspkx1x2{i}(:));
    imagesc(edge2center(x1binedges{i}), edge2center(x2binedges{i}),  log10(pspkx1x2{i}));
    set(gca,'clim', [-3 log10(mxmx)]);
    set(gca,'tickdir', 'out');
    axis('xy');
    xlabel('MID1 Proj');
    if i == 1
        ylabel('MID2 Proj');
    end
end % (for i)


%subplot(4,4,14);
%mxmx = max(pspkx1x2{2}(:));
%imagesc(edge2center(x1binedges{2}), edge2center(x2binedges{2}),  log10(pspkx1x2{2}));
%set(gca,'clim', [-3 log10(mxmx)]);
%set(gca,'tickdir', 'out');
%axis('xy');
%xlabel('MID1 Proj');
%
%
%subplot(4,4,15);
%mxmx = max(pspkx1x2{3}(:));
%imagesc(edge2center(x1binedges{3}), edge2center(x2binedges{3}),  log10(pspkx1x2{3}));
%set(gca,'clim', [-3 log10(mxmx)]);
%set(gca,'tickdir', 'out');
%axis('xy');
%xlabel('MID1 Proj');
%
%
%
%subplot(4,4,16);
%mxmx = max(pspkx1x2{4}(:));
%imagesc(edge2center(x1binedges{4}), edge2center(x2binedges{4}),  log10(pspkx1x2{3}));
%set(gca,'clim', [-3 log10(mxmx)]);
%set(gca,'tickdir', 'out');
%axis('xy');
%xlabel('MID1 Proj');



%
%hc = colorbar('EastOutside');
%
%set(hc,'ytick', [log10(0.001) log10(mxmx)], ...
%       'yticklabel', [0.001 roundsd(mxmx,3)], ...
%       'Position', [0.5 0.11 0.02 0.15], 'tickdir', 'out');
%
%

suptitle(strrep(iskfile, '_', '-'));

print_mfilename(mfilename);


return;


function mid_plot_projinfo_train_test_information(fraction, ifrac_train_mtx, ifrac_train_mn, ...
    ifrac_train_std, ifrac_test_mtx, ifrac_test_mn, ifrac_test_std, titlestr)
% mid_plot_projinfo_train_test_information - Show MID information results
%
% Information values were calculated for different fractions of the data.
% We fit a line to the information values versus the inverse of the data
% fraction. Information values were calculated for 8 different data set. 
% 4 training sets and 4 test sets.
%
% Input arguments:
%
% fraction : the fractions of the data that were used. Usually something
% like [80 85 90 92.5 95 97.5 100]
%
% ifrac_mtx_train : matrix of information values for each training set and
% all the data fractions. A 4xlength(fraction) matrix. Each row represents
% a training set, and each column is a data fraction.
%
% ifrac_mn_train : a vector of mean information values over the 4 training
% sets for each data fraction.
%
% ifrac_mtx_test : matrix of information values for each test set and
% all the data fractions. A 4xlength(fraction) matrix. Each row represents
% a test set, and each column is a data fraction.
%
% ifrac_mn_test : a vector of mean information values over the 4 test
% sets for each data fraction.
%
% Output arguments: none.
%


if ( nargin < 7 )
   error('You need 7 or 8 input args.');
end

if ( nargin ~= 8 )
   titlestr = [];
end


xmin = min(1./fraction);
xmax = max(1./fraction);
xrange = xmax - xmin;
xlim_vec = [xmin-0.1*xrange xmax+0.1*xrange];

x = 1./fraction;
xfit = linspace(0, max(x), 1000);

y = ifrac_train_mtx(1,:);
p1train = polyfit(x,y,1);
y1train = polyval(p1train, xfit);

y = ifrac_train_mtx(2,:);
p2train = polyfit(x,y,1);
y2train = polyval(p2train, xfit);

y = ifrac_train_mtx(3,:);
p3train = polyfit(x,y,1);
y3train = polyval(p3train, xfit);

y = ifrac_train_mtx(4,:);
p4train = polyfit(x,y,1);
y4train = polyval(p4train, xfit);


y = ifrac_test_mtx(1,:);
p1test = polyfit(x,y,1);
y1test = polyval(p1test, xfit);

y = ifrac_test_mtx(2,:);
p2test = polyfit(x,y,1);
y2test = polyval(p2test, xfit);

y = ifrac_test_mtx(3,:);
p3test = polyfit(x,y,1);
y3test = polyval(p3test, xfit);

y = ifrac_test_mtx(4,:);
p4test = polyfit(x,y,1);
y4test = polyval(p4test, xfit);




figure; 

% Plot the train data set information
subplot(2,5,1);
hold on;
plot( 1./fraction, ifrac_train_mtx(1,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y1train, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_train_mtx(1,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Train Set 1');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);
ylabel('Info [bits/spike]');



subplot(2,5,2);
hold on;
plot( 1./fraction, ifrac_train_mtx(2,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y2train, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_train_mtx(2,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Train Set 2');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);

subplot(2,5,3);
hold on;
plot( 1./fraction, ifrac_train_mtx(3,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y3train, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_train_mtx(3,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Train Set 3');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);

subplot(2,5,4);
hold on;
plot( 1./fraction, ifrac_train_mtx(4,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y4train, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_train_mtx(4,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Train Set 4');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);



% Plot the test data set information
subplot(2,5,6);
hold on;
plot( 1./fraction, ifrac_test_mtx(1,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y1test, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_test_mtx(1,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Test Set 1');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);
ylabel('Info [bits/spike]');


subplot(2,5,7);
hold on;
plot( 1./fraction, ifrac_test_mtx(2,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y2test, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_test_mtx(2,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Test Set 2');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);


subplot(2,5,8);
hold on;
plot( 1./fraction, ifrac_test_mtx(3,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y3test, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_test_mtx(3,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Test Set 3');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);
xlabel('Data Percentage');


subplot(2,5,9);
hold on;
plot( 1./fraction, ifrac_test_mtx(4,:), 'ko', 'markerfacecolor', 'k');
plot(xfit, y4test, 'k-');
xlim(xlim_vec);
ylim([0 1.1*max( ifrac_test_mtx(4,:) )]);
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Test Set 4');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);


y = ifrac_train_mn;
ptrain = polyfit(x,y,1);
ytrain = polyval(ptrain, xfit);

y = ifrac_test_mn;
ptest = polyfit(x,y,1);
ytest = polyval(ptest, xfit);


subplot(2,5,5);
hold on;
%errorbar( 1./fraction, ifrac_train_mn, ifrac_train_std/2, 'ko', 'markerfacecolor', 'k');
plot( 1./fraction, ifrac_train_mn, 'ko', 'markerfacecolor', 'k');
plot(xfit, ytrain, 'k-');
xlim(xlim_vec);
ylim( [0 1.1*max(ifrac_train_mn)+ max(ifrac_train_std)] );
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Train Set Mean');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);


subplot(2,5,10);
hold on;
%errorbar( 1./fraction, ifrac_test_mn, ifrac_test_std/2, 'ko', 'markerfacecolor', 'k');
plot( 1./fraction, ifrac_test_mn, 'ko', 'markerfacecolor', 'k');
plot(xfit, ytest, 'k-');
xlim(xlim_vec);
ylim( [0 1.1*max(ifrac_test_mn)+ max(ifrac_test_std)] );
set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
title('Test Set Mean');
xticklabel = roundsd(1./get(gca,'xtick'), 2);
set(gca,'xticklabel', xticklabel);

if (  ~isempty(titlestr) )
   suptitle(titlestr);
end

set(gcf, 'position', [133   295   760   420]);

return;


