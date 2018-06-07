function mid_cat_plot_fio_mid_database(fio, middb)


graphics_toolkit('gnuplot');
%graphics_toolkit('qt');

library('export_fig');


taxis = fio(1).time * 1000;
faxis = round(fio(1).freq / 1000);
nt = fio(1).tbins;
nf = fio(1).fbins;
cmap = brewmaps('rdbu',21);


for i = 1:length(fio)

    try

    exp = fio(i).exp;
    site = fio(i).site;
    unit = fio(i).unit;
    
    basefile = (sprintf('%s_site%s_unit%s',exp,site,unit));
    eps_file = sprintf('%s.eps', basefile);
    pdf_file = sprintf('%s.pdf', basefile); 

    if exist(fullfile('.', pdf_file),'file')
        fprintf('\n%s already exists! Skipping...\n', pdf_file)
        continue
    end       


    %h = figure('units', 'normalized','outerposition',[0 0 1 1]);
    h = figure('units', 'pixels');
 

    experiment = '2004-1-14';
    site = str2double(fio(i).site);

    if ( site == 16 )
        site = 616;
    elseif ( site == 2 )
        site = 602;
    end

    unit = str2double(fio(i).unit);

    ind = 0;

    for j = 1:length(middb)
        if ( middb(j).location == site )
            if ( middb(j).unit == unit )
                ind = j;
            end
        end
    end % (for j)

    if ~ind
        error('Did not match data structures.');
    end

    olddata = middb(ind);


    
    % plot New STA and nonlinearity
    subplot(3,3,1)
    hold on
    sta = fio(i).v_sta;
    imagesc(sta);
    minmin = min(min(sta));
    maxmax = max(max(sta));
    boundary = max([abs(minmin) abs(maxmax)]);
    axis xy
    xlim([-size(sta,2) size(sta,2)]);
    ylim([-5 size(sta,1)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', [], 'xticklabel', '');
    colormap(cmap);
    title('New STA');
    box on



    % plot Original STA and nonlinearity
    subplot(3,3,2)
    hold on
    sta = reshape(olddata.rpsta.filter, olddata.fbins, olddata.tbins);
    imagesc(sta);
    minmin = min(min(sta));
    maxmax = max(max(sta));
    boundary = max([abs(minmin) abs(maxmax)]);
    axis xy
    xlim([1 size(sta,2)]);
    ylim([1 size(sta,1)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', [], 'xticklabel', '');
    colormap(cmap);
    title('Orig STA');
    box on
   
  
    % Plot nonlinearities for both  
    subplot(3,3,3)
    nl = fio(i).sta_pspkx;
    xbins = fio(i).sta_xbins;
    hold on;
    plot(xbins, 2*nl, 'ko-', 'markerfacecolor', 'k', ...
        'markersize', 2, 'markerfacecolor', 'k');
    plot(olddata.rpx1pxpxt_sta.x, olddata.rpx1pxpxt_sta.ior_mean, 'ro-', ...
        'markersize', 2, 'markerfacecolor', 'r');
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    title('Nonlinearity');
    legend('New', 'Orig', 'location', 'NorthWest');


    % New MID1
    subplot(3,3,4)
    hold on
    mid1 = fio(i).v1;
    imagesc(mid1);
    minmin = min(min(mid1));
    maxmax = max(max(mid1));
    boundary = max([abs(minmin) abs(maxmax)]);
    axis xy
    xlim([-size(mid1,2) size(mid1,2)]);
    ylim([-5 size(mid1,1)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', [], 'xticklabel', '');
    title('New MID 1')
    box on


    % Original MID1
    subplot(3,3,5)
    hold on
    mid1 = reshape(olddata.rpdtest2_v1.filter, olddata.fbins, olddata.tbins);
    imagesc(mid1);
    minmin = min(min(mid1));
    maxmax = max(max(mid1));
    boundary = max([abs(minmin) abs(maxmax)]);
    axis xy
    xlim([1 size(mid1,2)]);
    ylim([1 size(mid1,1)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', [], 'xticklabel', '');
    title('Orig MID 1')
    box on
    

    subplot(3,3,6)
    new_fio = fio(i).mid1_pspkx;
    new_x = fio(i).mid1_xbins;

    orig_x = olddata.rpdx1x2px_pxt_2.x1;
    orig_fio = olddata.rpdx1x2px_pxt_2.ior1_mean;


    hold on;
    plot(new_x, 2*new_fio, 'ko-', 'markerfacecolor', 'k', ...
        'markersize', 2, 'markerfacecolor', 'k');
    plot(orig_x, orig_fio, 'ro-', ...
        'markersize', 2, 'markerfacecolor', 'r');
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);



    % New MID2
    subplot(3,3,7)
    hold on
    mid2 = fio(i).v2;
    imagesc(mid2);
    minmin = min(min(mid2));
    maxmax = max(max(mid2));
    boundary = max([abs(minmin) abs(maxmax)]);
    axis xy
    xlim([-size(mid2,2) size(mid2,2)]);
    ylim([-5 size(mid2,1)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', [], 'xticklabel', '');
    title('New MID 2')
    box on


    % Original MID2
    subplot(3,3,8)
    hold on
    mid2 = reshape(olddata.rpdtest2_v2.filter, olddata.fbins, olddata.tbins);
    imagesc(mid2);
    minmin = min(min(mid2));
    maxmax = max(max(mid2));
    boundary = max([abs(minmin) abs(maxmax)]);
    axis xy
    xlim([1 size(sta,2)]);
    ylim([1 size(sta,1)]);
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', [], 'xticklabel', '');
    title('Orig MID 2')
    box on
    

    subplot(3,3,9)
    new_fio = fio(i).mid2_pspkx;
    new_x = fio(i).mid2_xbins;

    orig_x = olddata.rpdx1x2px_pxt_2.x2;
    orig_fio = olddata.rpdx1x2px_pxt_2.ior2_mean;

    hold on;
    plot(new_x, 2*new_fio, 'ko-', 'markerfacecolor', 'k', ...
        'markersize', 2, 'markerfacecolor', 'k');
    plot(orig_x, orig_fio, 'ro-', ...
        'markersize', 2, 'markerfacecolor', 'r');
    set(gca,'tickdir', 'out', 'ticklength', [0.02 0.02]);

    
    suptitle(sprintf('%s-site%.0f-unit%.0f', exp, site, unit));
    
    print_mfilename(mfilename);

    set(h, 'position', [100 100 1300 1400]);

    fig2eps(eps_file);
    pause(1);
    crop = 0;
    append = 0;
    gray = 0;
    quality = 1000;
    eps2pdf(eps_file, pdf_file, crop, append, gray, quality);
    pause(1);
    fprintf('Figure saved in: %s\n\n', pdf_file);
    close all;

    catch
        fprintf('Found an error with %s\n\n', basefile);
    end
end




