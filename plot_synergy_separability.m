function options = plot_synergy_separability(info, fio, mid)
%
% options = plot_synergy_separability(info, fio, mid)
%
% options is a struct which holds plotting preferences. It can
% be used with the following command
%
% exportfig(gcf,'***.eps', options); 
%
% to make nice EPS figures which can be adjusted in
% Adobe Illustrator.
%
% caa 5/1/06

set(0,'defaultAxesFontName', 'Palatino')
set(0,'defaultAxesFontSize', 10);
options = struct('color','rgb','height',6,'width',2.75,'fontmode','fixed','fontsize', 8);


if ( length(info)~=length(mid) & length(info)~=length(fio) & ...
      length(fio) ~= length(mid) )
   error('info, fio, mid must be the same length.');
end


position = [];
info_sta = [];
info1 = [];
info2 = [];
info_both = [];
sta_asi = [];
v1_asi = [];
v2_asi = [];
sepindex = [];

for i = 1:length(info)

   position = [position info(i).position];
   info_sta = [info_sta mean(info(i).sta.information)];
   info1 = [info1 mean(info(i).mid1.information)];
   info2 = [info2 mean(info(i).mid2.information)];
   info_both = [info_both mean(info(i).mid12.information)];

   x = fio(i).sta.x;
   fx = fio(i).sta.fx;
   indexright = find(x>0);
   indexleft = find(x<0);
   right = sum( fx(indexright) - min(fx) );
   left = sum( fx(indexleft) - min(fx) );
   sta_asi = [sta_asi (right-left)/(right+left)];

   x = fio(i).v1.x;
   fx = fio(i).v1.fx;
   indexright = find(x>0);
   indexleft = find(x<0);
   right = sum( fx(indexright) - min(fx) );
   left = sum( fx(indexleft) - min(fx) );
   v1_asi = [v1_asi (right-left)/(right+left)];

   x = fio(i).v2.x;
   fx = fio(i).v2.fx;
   indexright = find(x>0);
   indexleft = find(x<0);
   right = sum( fx(indexright) - min(fx) );
   left = sum( fx(indexleft) - min(fx) );
   v2_asi = [v2_asi (right-left)/(right+left)];
 
   fx1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   fx2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   fx1x2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;
   [u,s,v] = svd(fx1x2);
   singvals = sum(s);
   eigvals = singvals .^ 2;
   sepindex = [sepindex 1-eigvals(1)/(sum(eigvals)+eps)];

end % (for)


% Error check to get rid of bad INF values:

index = ~isinf(info_sta) & ~isinf(info1) & ~isinf(info2) & ~isinf(info_both);

info_sta = info_sta(index);
info1 = info1(index);
info2 = info2(index);
info_both = info_both(index);
position = position(index);
sta_asi = sta_asi(index);
v1_asi = v1_asi(index);
v2_asi = v2_asi(index);
sepindex = sepindex(index);


info12_info1_info2 = 100 * info_both ./ (info1 + info2); % synergy index

synergy = 100 .* info_both ./ (info1 + info2); % synergy index

% info1_info12 = 100 * info1 ./ info_both; % linearity index
% info1_info1_info2 = 100 * info1 ./ (info1 + info2); % linearity index
% info12_info1_info2 = 100 * info_both ./ (info1 + info2); % synergy index
% info2_info1 = 100 * info2 ./ info1; % second filter contribution
% infosta_info1 = 100 * info_sta ./ info1; % sta compared to first mid
% infosta_info12 = 100 * info_sta ./ info_both; % sta compared to 
% info2_info12 = 100 * info2 ./ info_both;
% info2_info1_info2 = 100 * info2 ./ (info1 + info2);
% info12_minus_info1_info2 = 100 * ( info_both - (info1 + info2) ) / info_both;

min(synergy)
max(synergy)




%    plot(data, position, 'ko', 'markerfacecolor','k', 'markersize', 2);
%    posdata = min(min(unique(position))):max(max(unique(position)));
%    fitmodel = fit(position(:), data(:), 'smoothingspline', 'SmoothingParam', 0.1e-5);
%    [fitdata] = feval(fitmodel, posdata);
%    plot(fitdata, posdata, 'r-', 'linewidth', 2);


close all;

figure;

[beta, s] = polyfit(sepindex, synergy, 1);
x = linspace(min(sepindex), max(sepindex), 100);
y = polyval(beta, x);
[r, pvalue] = corrcoef(sepindex, synergy);

xtick = 0:0.1:0.7;
ytick = 75:25:200;

subplot(2,1,1);
hold on;
plot(sepindex, synergy, 'ko'); %, 'markerfacecolor','k', 'markersize', 3);
plot(x, y, 'k-');
xlim([0 0.7]);
ylim([75 205]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'tickdir', 'out');
xlabel('Separability Index');
ylabel('Synergy');
title(sprintf('r=%.3f, p=%.5f', r(1,2), pvalue(1,2)));

subplot(2,1,2);
hold on;
plot(sepindex, synergy, 'ks', 'markerfacecolor','k', 'markersize', 3);
plot(x, y, 'k-');
xlim([0 0.7]);
ylim([75 205]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', ytick, 'yticklabel', ytick);
set(gca,'tickdir', 'out');
xlabel('Separability Index');
ylabel('Synergy');
title(sprintf('r=%.3f, p=%.5f', r(1,2), pvalue(1,2)))

return;










layer2 = 200;
layer3 = 400;
layer4 = 800;
layer5 = 1100;
layer6 = 1500;
whitematter = 2000;

layerbins = 0:200:2400;
centerbins = 100:200:2300;

datatype{1} = 'sta_fx_asi';
datatype{2} = 'v1_fx_asi';
datatype{3} = 'v2_fx_asi';
datatype{4} = 'sepindex';

titlestr{1} = 'STA Asymmetry Index';
titlestr{2} = 'V_1 Asymmetry Index';
titlestr{3} = 'V_2 Asymmetry Index';
titlestr{4} = 'Separability Index for P(spike|x_1,x_2)';

close all;

for i = 1:4

   data = eval([ datatype{i} ]);

   figure;

   if ( i==4 )
      xmin = -0.1;
      xmax = 1.1;
   else
      xmin = -0.25;
      xmax = 1.1;
   end

   % Plot the data and the fit
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
   title([titlestr{i}]);
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
   he = errorbar(centerbins, mn, stermn, 'ko-');
   set(he, 'markersize', 3, 'markerfacecolor', 'k', 'linewidth', 1);
   set(gca,'xtick', centerbins, 'xticklabel', centerbins);
   set(gca,'tickdir', 'out');
   xlim([0 2000]);

%    ymax = max(mn+stermn)
%    ymin = min(mn-stermn)
%    yrange = ymax - ymin
%    ylim([ymin-0.05*yrange ymax+0.05*yrange]);

   if ( i==1 | i==2 )
      ylim([0.35 1]);
      ytick = 0.4:0.1:1;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
   elseif ( i==3 )
      ylim([-0.5 0.25]);
      ytick = -0.5:0.125:0.25;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
   elseif ( i==4 )
      ylim([0.2 0.4]);
      ytick = 0.2:0.05:0.4;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
   end

   subplot(3,1,3);
   bar(centerbins, mn);
   colormap([0.75 0.75 0.75]);
%    ymax = max(mn);
%    axis([-50 2050 ymin 1.05*ymax]);

   xlim([-50 2050]);
   if ( i==1 | i==2 )
      ylim([0.3 1]);
      ytick = 0.4:0.1:1;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
   elseif ( i==3 )
      ylim([-0.25 0.25]);
      ytick = -0.25:0.0625:0.25;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
   elseif ( i==4 )
      ylim([0.2 0.35]);
      ytick = 0.2:0.05:0.35;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
   end
   box off;
   set(gca,'tickdir', 'out');
   xlabel('Position (um)')
   ylabel([titlestr{i}]);


   mntestmat = [];
   cmb = nchoosek(1:length(datastr),2);
   for k = 1:size(cmb,1)
      [ha, pvalue, ci, stats] = ttest2(datastr( cmb(k,1) ).data, datastr( cmb(k,2) ).data, 0.05, 0);
      mntestmat = [mntestmat; datastr(cmb(k,1)).position datastr(cmb(k,2)).position ha pvalue];
   end % (for k)

   layerstr(i).type = titlestr{i};
   layerstr(i).mntest = mntestmat;
   layerstr(i).data = datastr;

end % (for j)



