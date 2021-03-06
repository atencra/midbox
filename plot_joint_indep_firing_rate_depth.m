function plot_joint_indep_firing_rate_depth(midpos)
%
% Figure 4 can be save using the following command:
%
% exportfig(gcf,'sta_v1_v2_asymmetry.eps', 'color', 'gray', ...
% 'height', 11, 'width', 8.5, 'fontmode', 'fixed', 'fontsize', 8);

prob1_plus2 = [];
prob12 = [];
position = [];
position_fr = [];

for i = 1:length(midpos)

   pos = midpos(i).position;
   
   fx1 = midpos(i).rpdx1x2px_pxt_2.ior1_mean;
   fx2 = midpos(i).rpdx1x2px_pxt_2.ior2_mean;
   fx1x2 = midpos(i).rpdx1x2px_pxt_2.ior12_mean;
   [u,s,v] = svd(fx1x2);
   singvals = sum(s);
   eigvals = singvals .^ 2;
   sepindex(i) = 1 - eigvals(1) / (sum(eigvals)+eps);

   position = [position pos];

   max1 = max(max(fx1));
   max2 = max(max(fx2));
   max12 = max(max(fx1x2));

   if ( ~isinf(max1) & ~isinf(max2) & ~isinf(max12) & ...
         max1<1 & max2<1 & max12<1 ) 

%       % Some debugging trick here. I don't remember how it originated,
%       % though there must have been some problem.
%       if ( max12 > 0.2499 & max12 < 0.25001 )
%          max12 = 0.25 + 0.5*rand(1);
%       end

      position_fr = [position_fr pos];
      prob1_plus2 = [prob1_plus2 max1+max2];
      prob12 = [prob12 max12];

   end

end % (for)


% [length(position) length(sepindex)]

% Change probabilities to firing rate. Each bin was 5 ms = 200 Hz
fr1plus2 = 200 .* prob1_plus2;
fr12 = 200 .* prob12;

fr_ratio = 100 .* fr12 ./ fr1plus2;


% prob12(prob12>0.2 & prob12<0.4)
% find(prob12>0.2499 & prob12<0.25001)
% pause


close all;

figure;

bins = 0:0.05:1;
xtick = 0:0.1:1;

subplot(1,2,1);
count = histc(sepindex, bins);
newcount = 100*count ./ sum(count);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.05 1.05 0 1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:10:50], 'yticklabel', [0:10:50]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('Percent of neurons (%)');
title('P(sp|x1,x2) InSep Index', 'fontsize', 10);


subplot(1,2,2)
hold on
xmax = max([max(fr12) max(fr1plus2)]);
plot(fr1plus2, fr12, 'ko')
fplot('x', 200.*[0 1.25], 'k-')
% fplot('x', [0 max([ max(prob12) max(prob1_plus2)])], 'k-')
% fplot('2*x', [0 max(info_both)], 'k-')
xlim(200.*[0.01 1.25]);
ylim(200.*[0.01 1.25]);
% xlim([0 1.1*xmax]);
% ylim([0 1.1*xmax]);
axis square
set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
xlabel('Sum of Peak Firing Rates')
ylabel('Joint Peak Firing Rate')
title('Comparison between Joint and Ind FR')

print_mfilename(mfilename);




boundary = 200:200:2000;
depth = 300:200:1900;

index = cell(1,length(boundary)-1);
fr_pop = cell(1,length(boundary)-1);
fr_mn = zeros(1,length(boundary)-1);
fr_se = zeros(1,length(boundary)-1);

count = zeros(1,length(boundary)-1);

for i = 1:length(boundary)-1

   index{i} = find( position_fr>boundary(i) & position_fr<=boundary(i+1) );
   count(i) = length(index{i});

   fr_pop{i} = fr_ratio(index{i});
   fr_mn(i) = mean( fr_ratio(index{i}) );
   fr_se(i) = std( fr_ratio(index{i}) ) / sqrt( length( fr_ratio(index{i}) ) );

end % (for i)


index = cell(1,length(boundary)-1);
sepindex_pop = cell(1,length(boundary)-1);
sepindex_mn = zeros(1,length(boundary)-1);
sepindex_se = zeros(1,length(boundary)-1);

for i = 1:length(boundary)-1

   index{i} = find( position>boundary(i) & position<=boundary(i+1) );

   sepindex_pop{i} = sepindex(index{i});
   sepindex_mn(i) = mean( sepindex(index{i}) );
   sepindex_se(i) = std( sepindex(index{i}) ) / sqrt( length( sepindex(index{i}) ) );

end % (for i)



data(1).title = '100 Joint FR / Ind FR';
data(1).pop = fr_pop;
data(1).mn = fr_mn;
data(1).se = fr_se;
data(1).ytick = [100:100:425];
data(1).yax = [100 425];

data(2).title = '2D Nonlinearity Inseparability Index';
data(2).pop = sepindex_pop;
data(2).mn = sepindex_mn;
data(2).se = sepindex_se;
data(2).ytick = [0.2:0.05:0.4];
data(2).yax = [0.2 0.4];


plot_fr_sepindex_depth_params(count, depth, data);

% nonlinearity_params_ttest(depth, data);





function plot_fr_sepindex_depth_params(count, depth, data)

figure;

for i = 1:length(data)

   pop = data(i).pop;
   mn = data(i).mn;
   se = data(i).se;
   ytick = data(i).ytick;
   yax = data(i).yax;
   titlestr = data(i).title;

   subplot(length(data),1,i);
   h = errorbar( depth, mn, se, 'ko-' );
   set(h, 'linewidth', 0.5, 'markersize', 3, 'markerfacecolor', [0.7 0.7 0.7], ...
   'markeredgecolor', [0.7 0.7 0.7]);

   hc = get(h,'children');
   edata = get(hc(2), 'xdata');
   edata(4:9:end) = edata(1:9:end) - 0;
   edata(7:9:end) = edata(1:9:end) - 0;
   edata(5:9:end) = edata(1:9:end) + 0;
   edata(8:9:end) = edata(1:9:end) + 0;
   set(hc(2), 'xdata', edata);

   xrange = max(depth)-min(depth);
   xmax = max(depth);
   xmin = min(depth);
   xlim([xmin - 0.055*xrange xmax+0.055*xmax]);
   set(gca,'xtick', depth, 'xticklabel', depth/1000);

   temp_plus = mn + se;
   temp_minus = mn - se;
   yrange = max(temp_plus) - min(temp_minus);
   ymax = max(temp_plus) + 0.06*yrange;
   ymin = min(temp_minus) - 0.06*yrange;
%    ylim([ymin ymax]);
   ylim(yax);
%    ytick = linspace(ymin,ymax,5);
   set(gca, 'ytick', ytick, 'yticklabel', ytick);

   set(gca, 'tickdir', 'out');
   set(gca, 'ticklength', [0.02 0.02]);

   set(gca, 'box', 'off');
   xlabel('Position (mm)');
   ylabel(titlestr);
   if ( i == 1 )
      title(sprintf('%.0f  ', count));
   end

end

print_mfilename(mfilename)

return;






function nonlinearity_params_ttest(depth, data)

layerstr(4) = struct('type', [], 'ttest', [], 'data', []); 


for i = 4:4

   pop = data(i).pop;
   titlestr = data(i).title;

   datastr(length(depth)) = struct('position', [], 'data', []);

   for j = 1:length(depth)
      datastr(j).position = depth(j);
      datastr(j).data = pop{j};
   end % (for j)

   ttestmat = [];
   rstestmat = [];
   kstestmat = [];
   cmb = nchoosek(1:length(datastr),2);

   ttestmat = zeros(size(cmb,1), 3);

   for k = 1:size(cmb,1)
      pos1 = datastr(cmb(k,1)).position;
      pos2 = datastr(cmb(k,2)).position;
      data1 = datastr( cmb(k,1) ).data;
      data2 = datastr( cmb(k,2) ).data;
      [htt, ptt] = ttest2(data1, data2, 0.05, 0);
      ttestmat(k,:) = [pos1/1000 pos2/1000 ptt];
   end % (for k)

   layerstr(i).type = titlestr;
   layerstr(i).ttest = ttestmat;
   layerstr(i).rstest = rstestmat;
   layerstr(i).kstest = kstestmat;
   layerstr(i).data = datastr;

end

% [ttestmat(:,3) 28.*ttestmat(:,3) kstestmat(:,3) 28.*kstestmat(:,3)]

ttestmat
[pid, pin] = fdr(ttestmat(:,3), 0.05);

pid
pin


[ttestmat(:,3) length(depth).* ttestmat(:,3)  <= 0.05 ttestmat(:,1:2)]

return;


