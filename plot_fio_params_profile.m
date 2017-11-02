function varargout = plot_fio_params_profile(fio_params, mid)
%
% [options, layerstr] = plot_fio_params_profile(fio_params, mid)
%
% options is a struct which holds plotting preferences. It can
% be used with the following command
%
% layerstr - struct holding statistical tests for the plotted information
%     quantities at different cortical depths.
%
% If layerstr is not listed as an output argument then only the profiles
% are displayed.
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

if ( nargout > 2 )
   error('Wrong number of output arguments.');
end

if ( length(fio_params) ~= length(mid) )
   error('fio_params and mid must be the same length.');
end



for i = 1:length(fio_params)

   position(i) = mid(i).position;

   % STA f(x) params

   sta_fx_max(i) = max( fio_params(i).sta.fx );
   sta_fx_mod(i) = 100 * ( fio_params(i).sta.fx(end-1) - fio_params(i).sta.fx(2) ) / fio_params(i).sta.fx(end-1);
   sta_dfx_mn(i) = ( fio_params(i).sta.dfx(end) + fio_params(i).sta.dfx(end-1) ) / 2;
   sta_dfx_end(i) = fio_params(i).sta.dfx(end);

   x = fio_params(i).sta.x;
   fx = fio_params(i).sta.fx;
   indexright = find(x>0);
   indexleft = find(x<0);
   right = sum( fx(indexright) - min(fx) );
   left = sum( fx(indexleft) - min(fx) );
   sta_fx_asi(i) = (right - left) / (right + left);

   sta_ef(i) = fio_params(i).sta.ef;
   sta_ef2(i) = fio_params(i).sta.ef2;
   sta_eflogf(i) = fio_params(i).sta.eflogf;
   sta_efdf(i) = fio_params(i).sta.efdf;

   % MID1 f(x) params
   v1_fx_max(i) = max( fio_params(i).v1.fx );
   v1_fx_mod(i) = 100 * ( fio_params(i).v1.fx(end-1) - fio_params(i).v1.fx(2) ) / fio_params(i).v1.fx(end-1);
   v1_dfx_mn(i) = ( fio_params(i).v1.dfx(end) + fio_params(i).v1.dfx(end-1) ) / 2;
   v1_dfx_end(i) = fio_params(i).v1.dfx(end);

   x = fio_params(i).v1.x;
   fx = fio_params(i).v1.fx;
   indexright = find(x>0);
   indexleft = find(x<0);
   right = sum( fx(indexright) - min(fx) );
   left = sum( fx(indexleft) - min(fx) );
   v1_fx_asi(i) = (right - left) / (right + left);

   v1_ef(i) = fio_params(i).v1.ef;
   v1_ef2(i) = fio_params(i).v1.ef2;
   v1_eflogf(i) = fio_params(i).v1.eflogf;
   v1_efdf(i) = fio_params(i).v1.efdf;

   % MID2 f(x) params
   v2_fx_mod(i) = 100 * ( fio_params(i).v2.fx(end-1) - fio_params(i).v2.fx(2) ) / fio_params(i).v2.fx(end-1);
   v2_ef(i) = fio_params(i).v2.ef;
   v2_ef2(i) = fio_params(i).v2.ef2;
   v2_eflogf(i) = fio_params(i).v2.eflogf;
   v2_efdf(i) = fio_params(i).v2.efdf;

   x = fio_params(i).v2.x;
   fx = fio_params(i).v2.fx;
   indexright = find(x>0);
   indexleft = find(x<0);
   right = sum( fx(indexright) - min(fx) );
   left = sum( fx(indexleft) - min(fx) );
   v2_fx_asi(i) = (right - left) / (right + left);
 
   r_sta_v1(i) = fio_params(i).r_sta_v1;
   r_sta_v2(i) = fio_params(i).r_sta_v2;
   r_v1_v2(i) = fio_params(i).r_v1_v2;

   sta = mid(i).rpsta.filter;
   sta(abs(sta)<0.5) = 0;
   v1 = mid(i).rpdtest2_v1.filter;
   v1(abs(v1)<0.5) = 0;
   [c,p] = corrcoef(sta,v1);
   si(i) = c(1,2);

   fx1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   fx2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   fx1x2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;
   [u,s,v] = svd(fx1x2);
   singvals = sum(s);
   eigvals = singvals .^ 2;
   sepindex(i) = 1 - eigvals(1) / (sum(eigvals)+eps);

end % (for)


layer2 = 200;
layer3 = 400;
layer4 = 800;
layer5 = 1100;
layer6 = 1500;
whitematter = 2000;

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



   if ( nargout==0 )

      if ( i==1 )
         figure;
      end

      if ( i==4 )
         xmin = 0;
         xmax = 0.6;
      elseif ( i==3 )
         xmin = -1;
         xmax = 1;
      else
         xmin = -0.25;
         xmax = 1.1;
      end

      % Plot the data and the fit
      subplot(2,2,i);
      hold on;
      plot([xmin xmax],[layer2 layer2],'-', 'color', [0.85 0.85 0.85]);
      plot([xmin xmax],[layer3 layer3],'-', 'color', [0.85 0.85 0.85]);
      plot([xmin xmax],[layer4 layer4],'-', 'color', [0.85 0.85 0.85]);
      plot([xmin xmax],[layer5 layer5],'-', 'color', [0.85 0.85 0.85]);
      plot([xmin xmax],[layer6 layer6],'-', 'color', [0.85 0.85 0.85]);
      plot(data, position, 'ks', 'markerfacecolor','k', 'markersize', 3);
      posdata = min(min(unique(position))):max(max(unique(position)));
      fitmodel = fit(position(:), data(:), 'smoothingspline', 'SmoothingParam', 0.1e-5);
      [fitdata] = feval(fitmodel, posdata);
      plot(fitdata, posdata, 'r-', 'linewidth', 1);
      xlim([xmin xmax]);
      ylim([0 2000]);
      box on;

      if ( ismember(i,[1 2]) )
         xtick = -0.25:0.25:1.0;
      elseif ( i==3 )
         xtick = -1:0.25:1;
      elseif ( ismember(i, [4]) )
         xtick = 0:0.1:0.6;
      end

      set(gca,'xtick', xtick, 'xticklabel', xtick);
      set(gca,'ytick',0:200:2000);

      set(gca,'ytick',0:200:2000);
      xlabel([titlestr{i}]);
      if ( i==1 | i==3 )
         ylabel('Position (um)')
      end
      set(gca,'ydir','rev');
      set(gca,'tickdir', 'out');

   else % ( nargout>0 )

      if ( i==4 )
         xmin = -0.1;
         xmax = 1.1;
      else
         xmin = -0.25;
         xmax = 1.1;
      end

      figure;

      % Plot the data and the fit
      subplot(3,1,1);
      hold on;
      plot(data, position, 'ko', 'markerfacecolor','k', 'markersize', 3);
      posdata = min(min(unique(position))):max(max(unique(position)));
      fitmodel = fit(position(:), data(:), 'smoothingspline', 'SmoothingParam', 0.1e-5);
      [fitdata] = feval(fitmodel, posdata);
      plot(fitdata, posdata, 'r-', 'linewidth', 1);
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
%       ymax = max(mn);
%       axis([-50 2050 ymin 1.05*ymax]);

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
      xlabel('Position (um)');
      ylabel([titlestr{i}]);

      mntestmat = [];
      cmb = nchoosek(1:length(datastr),2);
      for k = 1:size(cmb,1)
         pos1 = datastr(cmb(k,1)).position;
         pos2 = datastr(cmb(k,2)).position;
         data1 = datastr( cmb(k,1) ).data;
         data2 = datastr( cmb(k,2) ).data;
         mn1 = mean(data1);
         mn2 = mean(data2);
         [ha, pvalue, ci, stats] = ttest2(data1, data2, 0.05, 0);
         mntestmat = [mntestmat; pos1./1000 pos2./1000 length(data1) length(data2) mn1 mn2 ha pvalue];
      end % (for k)

      layerstr(i).type = titlestr{i};
      layerstr(i).mntest = mntestmat;
      layerstr(i).data = datastr;

      varargout{1} = options;
      if ( nargout==2 )
         varargout{2} = layerstr;
      end

   end

end % (for j)



