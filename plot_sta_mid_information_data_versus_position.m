function varargout = plot_sta_mid_information_data_versus_position(info)
% plot_sta_mid_information_data_versus_positions -compare information values
% for sta and mids, and also determine degree of nonlinearity and synergy.
%
% layerstr = plot_sta_mid_information_data_versus_position(info)
%
% info - struct array holding information data
%
% layerstr - struct holding statistical tests for the plotted information
%     quantities at different cortical depths.
%
% If layerstr is not listed as an output argument then only the profiles
% are displayed.
%
% You can make pretty pre-publication eps plots via the following commands:
% options = struct('color','rgb','height',6,'width',2.75,'fontmode','fixed','fontsize', 8);
% exportfig(gcf,'*_profile.eps', options);
%
% caa 6/1/06

set(0,'defaultAxesFontName', 'Palatino')
set(0,'defaultAxesFontSize', 10)

position = [];
info_sta = [];
info1 = [];
info2 = [];
info_both = [];

for i = 1:length(info)

   position = [position info(i).position];
   info_sta = [info_sta mean(info(i).sta.information)];
   info1 = [info1 mean(info(i).mid1.information)];
   info2 = [info2 mean(info(i).mid2.information)];
   info_both = [info_both mean(info(i).mid12.information)];

end % (for i)


% Error check to get rid of bad INF values:

index = ~isinf(info_sta) & ~isinf(info1) & ~isinf(info2) & ~isinf(info_both);

info_sta = info_sta(index);
info1 = info1(index);
info2 = info2(index);
info_both = info_both(index);
position = position(index);

% info1_info12 = 100 * info1 ./ info_both; % linearity index
info1_info1_info2 = 100 * info1 ./ (info1 + info2); % linearity index
% info12_info1_info2 = 100 * info_both ./ (info1 + info2); % synergy index
% info2_info1 = 100 * info2 ./ info1; % second filter contribution
% infosta_info1 = 100 * info_sta ./ info1; % sta compared to first mid
% infosta_info12 = 100 * info_sta ./ info_both; % sta compared to 
% info2_info12 = 100 * info2 ./ info_both;
info2_info1_info2 = 100 * info2 ./ (info1 + info2);
% info12_minus_info1_info2 = 100 * ( info_both - (info1 + info2) ) / info_both;

% [info1_info1_info2(:) info2_info1_info2(:) info1_info1_info2(:)+info2_info1_info2(:)]
% pause

layer2 = 200;
layer3 = 400;
layer4 = 800;
layer5 = 1100;
layer6 = 1500;
whitematter = 2000;

titlestr{1} = '100 * I(v_1) / I(v_1,v_2)';
titlestr{2} = '100 * I(v_1) / ( I(v_1) + I(v_2) )';
titlestr{3} = '100 * I(v_1,v_2) / ( I(v_1) + I(v_2) )';
titlestr{4} = '100 * I(v_2) / I(v_1)';
titlestr{5} = '100 * I(sta) / I(v_1)';
titlestr{6} = '100 * I(sta) / I(v_1,v_2)';
titlestr{7} = '100 * I(v_2) / I(v_1,v_2)';
titlestr{8} = '100 * I(v_2) / ( I(v_1) + I(v_2) )';
titlestr{9} = '100 * ( I(v_1,v_2) - I(v_1) - I(v_2) ) / I(v_1,v_2)';

layerbins = 0:200:2400;
centerbins = 100:200:2300;

close all;

layerstr(9) = struct('type', [], 'mntest', [], 'data', []);

for i = 1:9

   if ( i == 1 )
      data = 100 * info1 ./ info_both;
   elseif ( i == 2 )
      data = 100 * info1 ./ (info1 + info2);
   elseif ( i == 3 )
      data = 100 * info_both ./ (info1 + info2);
   elseif ( i == 4 )
      data = 100 * info2 ./ info1;
   elseif ( i == 5 )
      data = 100 * info_sta ./ info1;
   elseif ( i == 6 )
      data = 100 * info_sta ./ info_both;
   elseif ( i == 7 )
      data = 100 * info2 ./ info_both;
   elseif ( i == 8 )
      data = 100 * info2 ./ (info1 + info2);
   elseif ( i==9 )
      data = 100 * ( info_both - info1 - info2 ) ./ info_both;
   end

   xmin = -5;
   if ( i==1 )
      xmin = 35;
      xmax = 105;
   elseif ( i==2 )
      xmin = 40;
      xmax = 105;
   elseif ( i==3 )
      xmin = 50;
      xmax = 225;
   elseif ( i==4 )
      xmin = 0;
      xmax = 50;
   elseif ( i==5 )
      xmin = 35;
      xmax = 105;
   elseif ( i==6 )
      xmin = 15;
      xmax = 105;
   elseif ( i==7 )
      xmin = -5;
      xmax = 55;
   elseif ( i==8 )
      xmin = -5;
      xmax = 55;
   elseif ( i==9 )
      xmin = -30;
      xmax = 55;
   end


   if ( nargout==0 )

      % Plot the data

      if ( i==1 )
         figure; 
      end;

      subplot(3,3,i);
      hold on;
      plot(data, position, 'ko', 'markerfacecolor','k', 'markersize', 3);
      posdata = min(min(unique(position))):max(max(unique(position)));
      fitmodel = fit(position(:), data(:), 'smoothingspline', 'SmoothingParam', 0.1e-5);
      [fitdata] = feval(fitmodel, posdata);
      plot(fitdata, posdata, 'r-', 'linewidth', 1);

      xlim([xmin xmax]);
      ylim([0 2000]);

      box on;
      plot([xmin xmax],[layer2 layer2],'k:');
      plot([xmin xmax],[layer3 layer3],'k:');
      plot([xmin xmax],[layer4 layer4],'k:');
      plot([xmin xmax],[layer5 layer5],'k:');
      plot([xmin xmax],[layer6 layer6],'k:');

      if ( ismember(i,[1 2 5]) )
         xtick = 40:20:100;
      elseif ( i==3 )
         xtick = 50:50:200;
      elseif ( ismember(i, [4 7 8]) )
         xtick = 0:10:50;
      elseif ( i==6 )
         xtick = 20:20:100;
      elseif ( i==9 )
         xtick = -20:20:40;
      end

      set(gca,'xtick', xtick, 'xticklabel', xtick);
      set(gca,'ytick',0:200:2000);

%       set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);

      xlabel([titlestr{i}]);

      if ( i==1 | i==4 | i==7 )
         ylabel('Position (um)')
      end

      set(gca,'ydir','rev');
      set(gca,'tickdir', 'out');

   else % nargout==1

      % Plot the data

      figure;

      subplot(3,1,1);
      hold on;
      plot(data, position, 'ko', 'markerfacecolor','k', 'markersize', 2);
      posdata = min(min(unique(position))):max(max(unique(position)));
      fitmodel = fit(position(:), data(:), 'smoothingspline', 'SmoothingParam', 0.1e-5);
      [fitdata] = feval(fitmodel, posdata);
      plot(fitdata, posdata, 'r-', 'linewidth', 2);

      xlim([xmin xmax]);
      ylim([0 2400]);

      box on;
      plot([xmin xmax],[layer2 layer2],'k:');
      plot([xmin xmax],[layer3 layer3],'k:');
      plot([xmin xmax],[layer4 layer4],'k:');
      plot([xmin xmax],[layer5 layer5],'k:');
      plot([xmin xmax],[layer6 layer6],'k:');
      plot([xmin xmax],[whitematter whitematter],'k-');

      set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
      xlabel([titlestr{i}]);
      ylabel('Position (um)')
      set(gca,'ydir','rev');
      set(gca,'tickdir', 'out');


      % Now we plot the mean parameters values as a function of position. These
      % results will be cleaner than the previous subplot.

      mn = [];
      sd = [];
      stermn = [];

      layerbins = 0:200:2000;
      centerbins = 100:200:1900;
      datastr = struct('position', [], 'data', [], 'mn', [], 'sd', [], 'stermn', []);

      for j = 1:(length(layerbins)-1)
         index = find(position > layerbins(j) & position <= layerbins(j+1) );
         temp = data(index);
         mn = [mn mean(temp)];
         sd = [sd std(temp)];
         stermn = [stermn std(temp)/sqrt(length(temp))];

         datastr(j).position = centerbins(j);
         datastr(j).data = temp;
         datastr(j).mn = mn;
         datastr(j).sd = sd;
         datastr(j).stermn = stermn;      
      end % (for j)


      subplot(3,1,2);
      hold on;
      h = errorbar(centerbins, mn, stermn, 'ko-');
      set(h, 'markersize', 3, 'markerfacecolor', 'k', 'linewidth', 1);
      set(gca,'xtick', centerbins, 'xticklabel', centerbins);
      set(gca,'tickdir', 'out');

%       edata = get(h(1), 'xdata');
%       edata(4:9:end) = edata(1:9:end) - 0;
%       edata(7:9:end) = edata(1:9:end) - 0;
%       edata(5:9:end) = edata(1:9:end) + 0;
%       edata(8:9:end) = edata(1:9:end) + 0;
%       set(h(1), 'xdata', edata);

      xlim([0 2000]);


      if ( i==1 )
         ylim([50 70]);
         ytick = 50:5:70;
      elseif ( i==2 )
         ylim([50 90]);
         ytick = 50:10:90;
      elseif ( i==3 )
         ylim([100 160]);
         ytick = 100:10:160;
      elseif ( i==4 )
         ylim([15 45]);
         ytick = 15:10:45;
      elseif ( i==5 )
         ylim([50 80]);
         ytick = 50:10:80;
      elseif ( i==6 )
         ylim([25 55]);
         ytick = 25:10:55;
      elseif ( i==7 )
         ylim([10 30]);
         ytick = 10:5:30;
      elseif ( i==8 )
         ylim([15 32.5]);
         ytick = 15:5:30;
      elseif ( i==9 )
         ylim([7.5 37.5]);
         ytick = 10:5:35;
      end
      set(gca,'ytick', ytick, 'yticklabel', ytick);


      subplot(3,1,3);
      bar(centerbins, mn);
      colormap([0.75 0.75 0.75]);
      xlim([-50 2050]);
      if ( i==1 )
         ylim([50 67.5]);
         ytick = 50:5:65;
      elseif ( i==2 )
         ylim([50 90]);
         ytick = 50:10:90;
      elseif ( i==3 )
         ylim([100 160]);
         ytick = 100:10:160;
      elseif ( i==4 )
         ylim([15 40]);
         ytick = 15:5:40;
      elseif ( i==5 )
         ylim([50 80]);
         ytick = 50:10:80;
      elseif ( i==6 )
         ylim([25 50]);
         ytick = 25:5:50;
      elseif ( i==7 )
         ylim([10 25]);
         ytick = 10:5:25;
      elseif ( i==8 )
         ylim([15 30]);
         ytick = 15:5:30;
      elseif ( i==9 )
         ylim([10 35]);
         ytick = 10:5:35;
      end
      set(gca,'ytick', ytick, 'yticklabel', ytick);
      box off;
      set(gca,'tickdir', 'out');
      xlabel('Position (um)')
      ylabel([titlestr{i}]);

      mntestmat = [];
      cmb = nchoosek(1:length(datastr),2);
      mntestmat = zeros(size(cmb,1), 6);
      for k = 1:size(cmb,1)
         pos1 = datastr(cmb(k,1)).position;
         pos2 = datastr(cmb(k,2)).position;
         data1 = datastr( cmb(k,1) ).data;
         data2 = datastr( cmb(k,2) ).data;
         mn1 = mean(data1);
         mn2 = mean(data2);
         [ha, pvalue, ci, stats] = ttest2(data1, data2, 0.05, 0);
         mntestmat(k,:) = [pos1./1000 pos2./1000 mn1 mn2 ha pvalue];
      end % (for k)

      layerstr(i).type = titlestr{i};
      layerstr(i).mntest = mntestmat;
      layerstr(i).data = datastr;

   end % (if)

end % (for i)


if ( nargout==1 )
   varargout{1} = layerstr;
end % (if)

return;


% 
% 
% close all;
% 
% data(1).title = 'STA ASI';
% data(1).pop = sta_asi_pop;
% data(1).mn = sta_asi_mn;
% data(1).se = sta_asi_se;
% data(1).ytick = [0.65:0.05:0.9];
% data(1).yax = [0.65 0.9];
% 
% data(2).title = 'MID_1 ASI';
% data(2).pop = v1_asi_pop;
% data(2).mn = v1_asi_mn;
% data(2).se = v1_asi_se;
% data(2).ytick = [0.4:0.1:0.80];
% data(2).yax = [0.35 0.75];
% 
% data(3).title = 'MID_2 ASI';
% data(3).pop = v2_asi_pop;
% data(3).mn = v2_asi_mn;
% data(3).se = v2_asi_se;
% data(3).ytick = [-0.2:0.1:0.2];
% data(3).yax = [-0.2 0.2];
% 
% data(4).title = '2D Nonlinearity Inseparability Index';
% data(4).pop = sepindex_pop;
% data(4).mn = sepindex_mn;
% data(4).se = sepindex_se;
% data(4).ytick = [0.2:0.05:0.4];
% data(4).yax = [0.2 0.4];
% 
% 
% plot_depth_nonlinearity_params_circle_niceyticks(count, depth, data);
% 
% nonlinearity_params_ttest(depth, data);
% 
% 
% function plot_depth_nonlinearity_params_circle_niceyticks(count, depth, data)
% 
% figure;
% 
% for i = 1:4
% 
%    pop = data(i).pop;
%    mn = data(i).mn;
%    se = data(i).se;
%    ytick = data(i).ytick;
%    yax = data(i).yax;
%    titlestr = data(i).title;
% 
%    subplot(4,1,i);
%    h = errorbar( depth, mn, se, 'ko-' );
%    set(h, 'linewidth', 0.5, 'markersize', 3, 'markerfacecolor', [0.7 0.7 0.7], ...
%    'markeredgecolor', [0.7 0.7 0.7]);
% 
% get(h(2))
% 
%    edata = get(h(1), 'xdata');
%    edata(4:9:end) = edata(1:9:end) - 0;
%    edata(7:9:end) = edata(1:9:end) - 0;
%    edata(5:9:end) = edata(1:9:end) + 0;
%    edata(8:9:end) = edata(1:9:end) + 0;
%    set(h(1), 'xdata', edata);
% 
% 
% %    hc = get(h,'children')
% %    edata = get(hc(2), 'xdata');
% %    edata(4:9:end) = edata(1:9:end) - 0;
% %    edata(7:9:end) = edata(1:9:end) - 0;
% %    edata(5:9:end) = edata(1:9:end) + 0;
% %    edata(8:9:end) = edata(1:9:end) + 0;
% %    set(hc(2), 'xdata', edata);
% 
%    xrange = max(depth)-min(depth);
%    xmax = max(depth);
%    xmin = min(depth);
%    xlim([xmin - 0.055*xrange xmax+0.055*xmax]);
%    set(gca,'xtick', depth, 'xticklabel', depth/1000);
% 
%    temp_plus = mn + se;
%    temp_minus = mn - se;
%    yrange = max(temp_plus) - min(temp_minus);
%    ymax = max(temp_plus) + 0.06*yrange;
%    ymin = min(temp_minus) - 0.06*yrange;
% %    ylim([ymin ymax]);
%    ylim(yax);
% %    ytick = linspace(ymin,ymax,5);
%    set(gca, 'ytick', ytick, 'yticklabel', ytick);
% 
%    set(gca, 'tickdir', 'out');
%    set(gca, 'ticklength', [0.02 0.02]);
% 
%    set(gca, 'box', 'off');
%    xlabel('Position (mm)');
%    ylabel(titlestr);
%    if ( i == 1 )
%       title(sprintf('%.0f  ', count));
%    end
% 
% end
% 
% return;
% 
% 
% 
% 
% 
% function nonlinearity_params_ttest(depth, data)
% 
% layerstr(4) = struct('type', [], 'ttest', [], 'data', []); 
% 
% 
% for i = 4:4
% 
%    pop = data(i).pop;
%    titlestr = data(i).title;
% 
%    datastr(length(depth)) = struct('position', [], 'data', []);
% 
%    for j = 1:length(depth)
%       datastr(j).position = depth(j);
%       datastr(j).data = pop{j};
%    end % (for j)
% 
%    ttestmat = [];
%    rstestmat = [];
%    kstestmat = [];
%    cmb = nchoosek(1:length(datastr),2);
% 
%    ttestmat = zeros(size(cmb,1), 3);
% 
%    for k = 1:size(cmb,1)
%       pos1 = datastr(cmb(k,1)).position;
%       pos2 = datastr(cmb(k,2)).position;
%       data1 = datastr( cmb(k,1) ).data;
%       data2 = datastr( cmb(k,2) ).data;
%       [htt, ptt] = ttest2(data1, data2, 0.05, 0);
%       ttestmat(k,:) = [pos1/1000 pos2/1000 ptt];
%    end % (for k)
% 
%    layerstr(i).type = titlestr;
%    layerstr(i).ttest = ttestmat;
%    layerstr(i).rstest = rstestmat;
%    layerstr(i).kstest = kstestmat;
%    layerstr(i).data = datastr;
% 
% end
% 
% % [ttestmat(:,3) 28.*ttestmat(:,3) kstestmat(:,3) 28.*kstestmat(:,3)]
% 
% % ttestmat
% % [pid, pin] = fdr(ttestmat(:,3), 0.05);
% % 
% % pid
% % pin
% 
% 
% [ttestmat(:,3) length(depth).* ttestmat(:,3)  <= 0.05 ttestmat(:,1:2)]
% 
% return;
% 
% 
% 
% varargout{1} = options;
% if ( nargout==2 )
%    varargout{2} = layerstr;
% end
% 
















