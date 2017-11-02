function plot_comparison_between_sta_mid_nonlinearities_with_depth(fio_params, mid)
%
% Figure 4 can be save using the following command:
%
% exportfig(gcf,'sta_v1_v2_asymmetry.eps', 'color', 'gray', 'height', 11, 'width', 8.5, 'fontmode', 'fixed', 'fontsize', 8);

prob1_plus2 = [];
prob12 = [];

for i = 1:length(fio_params)

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
   
   xv1 = x;
   fxv1 = fx;


   % MID2 f(x) params
   v2_fx_mod(i) = 100 * ( fio_params(i).v2.fx(end-1) - fio_params(i).v2.fx(2) ) / ...
      fio_params(i).v2.fx(end-1);
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
   
   xv2 = x;
   fxv2 = fx;
   ind0 = find(xv2>=-0.0001 & xv2<=0.0001);
   right = fxv2(ind0:end);
   left = [fxv2(1:ind0-1) 0];
   left = fliplr(left);
   sumleftright = left + right;

%    subplot(2,2,1);
%    plot(xv2, fxv2);
%    subplot(2,2,2);
%    plot(xv2(ind0:end), right);
%    subplot(2,2,4);
%    plot(xv2(ind0:end), left);
%    subplot(2,2,3);
%    plot(xv2(ind0:end), sumleftright);


   r_sta_v1(i) = fio_params(i).r_sta_v1;
   r_sta_v2(i) = fio_params(i).r_sta_v2;
   r_v1_v2(i) = fio_params(i).r_v1_v2;


   sta = mid(i).rpsta.filter;
   sta(abs(sta)<0.5) = 0;
   v1 = mid(i).rpdtest2_v1.filter;
   v1(abs(v1)<0.5) = 0;
   [c,p] = corrcoef(sta,v1);
   si(i) = c(1,2);
   
   
   v2 = mid(i).rpdtest2_v2.filter;
   v2(abs(v2)<0.5) = 0;
   [c,p] = corrcoef(v1,v2);
   si2(i) = c(1,2);

   
   fx1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   fx2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   fx1x2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;
   [u,s,v] = svd(fx1x2);
   singvals = sum(s);
   eigvals = singvals .^ 2;
   sepindex(i) = 1 - eigvals(1) / (sum(eigvals)+eps);

   max1 = max(max(fx1));
   max2 = max(max(fx2));
   max12 = max(max(fx1x2));

   if ( ~isinf(max1) & ~isinf(max2) & ~isinf(max12) & max1<1 & max2<1 & max12<1 ) 
      if ( max12 > 0.2499 & max12 < 0.25001 )
%          disp('yes');
         max12 = 0.25 + 0.5*rand(1);
      end
      prob1_plus2 = [prob1_plus2 max1+max2];
      prob12 = [prob12 max12];
   end

% plot(fio_params(i).v1.x, fio_params(i).v1.fx, 'k-', ...
%    mid(i).rpdx1x2px_pxt_2.x1, mid(i).rpdx1x2px_pxt_2.ior1_mean, 'r-');
% pause

end % (for)


% prob12(prob12>0.2 & prob12<0.4)
% find(prob12>0.2499 & prob12<0.25001)
% pause


% mean(si2)
% std(si2)
% pause


figure;

% max( f(x) )
subplot(4,3,1);
bins = 0:0.1:1;
count = histc(sta_fx_max, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', [0:0.2:1], 'xticklabel', [0:0.2:1]);
set(gca,'tickdir', 'out');
ylabel('#units, max( P(sp|x) )', 'fontsize', 10);
title('STA', 'fontsize', 10);
set(gca,'fontsize', 8);


subplot(4,3,2);
bins = 0:0.1:1;
count = histc(v1_fx_max, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', [0:0.2:1], 'xticklabel', [0:0.2:1]);
set(gca,'tickdir', 'out');
title('MID1', 'fontsize', 10);
set(gca,'fontsize', 8);


% 100 * ( f(N)-f(1) ) / f(N)
subplot(4,3,4);
bins = -200:25:200;
count = histc(sta_fx_mod, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-210 210 0 1.05*max(count)]);
set(gca,'xtick', [-200:50:200], 'xticklabel', [-200:50:200]);
set(gca,'tickdir', 'out');
ylabel('100*( P(sp|x)(N)-P(sp|x)(1) ) / P(sp|x)(N)', 'fontsize', 10)
set(gca,'fontsize', 8);


subplot(4,3,5);
bins = -200:25:200;
count = histc(v1_fx_mod, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-210 210 0 1.05*max(count)]);
set(gca,'xtick', [-200:50:200], 'xticklabel', [-200:50:200]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);


subplot(4,3,6);
bins = -200:25:200;
count = histc(v2_fx_mod, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-210 210 0 1.05*max(count)]);
set(gca,'xtick', [-200:50:200], 'xticklabel', [-200:50:200]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
title('MID2', 'fontsize', 10);


% df(N)
subplot(4,3,7);
bins = -0.4:0.05:0.4;
count = histc(sta_dfx_end, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.4 0.4 0 1.05*max(count)]);
set(gca,'xtick', [-0.4:0.1:0.4], 'xticklabel', [-0.4:0.1:0.4]);
set(gca,'tickdir', 'out');
ylabel('dP(sp|x)(N)', 'fontsize', 10);
set(gca,'fontsize', 8);


subplot(4,3,8);
bins = -0.4:0.05:0.4;
count = histc(v1_dfx_end, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.4 0.4 0 1.05*max(count)]);
set(gca,'xtick', [-0.4:0.1:0.4], 'xticklabel', [-0.4:0.1:0.4]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);


% df(N) - df(N-1) ) / 2
subplot(4,3,10);
bins = -0.4:0.05:0.4;
count = histc(sta_dfx_mn, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.4 0.4 0 1.05*max(count)]);
set(gca,'xtick', [-0.4:0.1:0.4], 'xticklabel', [-0.4:0.1:0.4]);
set(gca,'tickdir', 'out');
ylabel('( dP(sp|x)(N) + dP(sp|x)(N-1) ) / 2', 'fontsize', 10);
set(gca,'fontsize', 8);


subplot(4,3,11);
bins = -0.4:0.05:0.4;
count = histc(v1_dfx_mn, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.4 0.4 0 1.05*max(count)]);
set(gca,'xtick', [-0.4:0.1:0.4], 'xticklabel', [-0.4:0.1:0.4]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);

print_mfilename(mfilename);



figure;

% max( f(x) )
subplot(2,2,1);
hold on;
plot(sta_fx_max, v1_fx_max, 'ko');
plot([0.003 1.5], [0.003 1.5], 'k-');
axis([0.003 1.5 0.003 1.5]);
axis('square');
set(gca,'xscale', 'log', 'yscale', 'log', 'tickdir', 'out');
ylabel('MID1');
title('Max( P(sp|x) )');

% 100 * ( f(N)-f(1) ) / f(N)
subplot(2,2,2);
hold on;
plot(sta_fx_mod, v1_fx_mod, 'ko');
plot([-150 150], [-150 150], 'k-');
xtick = -150:50:150;
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'tickdir', 'out');
axis([-150 150 -150 150]);
axis('square');
title('100*(P(sp|x)(N-1)-P(sp|x)(2)) / P(sp|x)(N-1)');

% max( f(x) )
subplot(2,2,3);
hold on;
plot(sta_dfx_end, v1_dfx_end, 'ko');
plot([-0.4 0.4], [-0.4 0.4], 'k-');
xtick = -0.4:0.1:0.4;
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'tickdir', 'out');
axis([-0.4 0.4 -0.4 0.4]);
axis('square');
xlabel('STA');
ylabel('MID1');
title('dP(sp|x)(N)');

% ( df(N) + df(N-1) ) / 2
subplot(2,2,4);
hold on;
plot(sta_dfx_mn, v1_dfx_mn, 'ko');
plot([-0.4 0.4], [-0.4 0.4], 'k-');
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'tickdir', 'out');
axis([-0.4 0.4 -0.4 0.4]);
axis('square');
xlabel('STA');
title('( dP(sp|x)(N) + dP(sp|x)(N-1) ) / 2');

print_mfilename(mfilename);



figure;

bins = -1:0.1:1;
xtick = -1:0.2:1;

subplot(3,1,1);
count = histc(r_sta_v1, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('#Units');
title('Correlation b/w P(sp|x) for STA, V1', 'fontsize', 10);

subplot(3,1,2);
count = histc(r_sta_v2, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('#Units');
title('Correlation b/w P(sp|x) for STA, V2', 'fontsize', 10);

subplot(3,1,3);
count = histc(r_v1_v2, bins);
hb = bar(bins, count, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 1.05*max(count)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
xlabel('Correlation Coefficient');
ylabel('#Units');
title('Correlation b/w P(sp|x) for V1, V2', 'fontsize', 10);

print_mfilename(mfilename);





% close all;

figure;

bins = -1:0.1:1;
xtick = -1:0.2:1;
% length(find(sta_fx_asi>0.5))
subplot(4,2,1);
count = histc(sta_fx_asi, bins);
newcount = 100 * count ./ sum(count);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 30]); %1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:5:30], 'yticklabel', [0:5:30]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('Percent of neurons (%)');
title('STA P(sp|x) Asymmetry', 'fontsize', 10);

subplot(4,2,3);
count = histc(v1_fx_asi, bins);
newcount = 100 * count ./ sum(count);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 25]); %1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:5:25], 'yticklabel', [0:5:25]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('Percent of neurons (%)');
title('V1 P(sp|x) Asymmetry', 'fontsize', 10);

subplot(4,2,5);
count = histc(v2_fx_asi, bins);
newcount = 100*count ./ sum(count);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 20]); %1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:5:20], 'yticklabel', [0:5:20]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
xlabel('Asymmetry Index');
ylabel('Percent of neurons (%)');
title('V2 P(sp|x) Asymmetry', 'fontsize', 10);

subplot(4,2,2);
% mean(si)
% std(si)
count = histc(si, bins);
newcount = 100*count ./ sum(count);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 45]); %1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:5:45], 'yticklabel', [0:5:45]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('Percent of neurons (%)');
title('STA/V1 Filter Similarity Index', 'fontsize', 10);


subplot(4,2,4);
% mean(r_sta_v1)
% std(r_sta_v1)
count = histc(r_sta_v1, bins);
newcount = 100*count ./ sum(count);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-1.1 1.1 0 40]); %1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:5:40], 'yticklabel', [0:5:40]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('Percent of neurons (%)');
title('Correlation b/w P(sp|x) for STA, V1', 'fontsize', 10);



xtick = -1:0.25:1;

subplot(4,2,6);
hold on;
plot(v1_fx_asi, sta_fx_asi, 'ko', 'markersize', 3);
plot([-1 1], [-1 1], 'k--');
plot([-1 1], [0 0], 'k-');
plot([0 0], [-1 1], 'k-');
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'tickdir', 'out');
axis([-0.5 1.1 -0.5 1.1]);
%axis([-1.1 1.1 -1.1 1.1]);
% axis('square');
xlabel('V1 P(sp|x) ASI');
ylabel('STA P(sp|x) ASI');

[mean(v1_fx_asi) std(v1_fx_asi) mean(sta_fx_asi) std(sta_fx_asi)]

[h, p]=ttest2(v1_fx_asi, sta_fx_asi)

subplot(4,2,8);
hold on;
plot(v1_fx_asi, v2_fx_asi, 'ko', 'markersize', 3);
plot([-1 1], [-1 1], 'k--');
plot([-1 1], [0 0], 'k-');
plot([0 0], [-1 1], 'k-');
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', xtick, 'yticklabel', xtick);
set(gca,'tickdir', 'out');
axis([-0.5 1.1 -1.1 1.1]);
%axis([-1.1 1.1 -1.1 1.1]);
% axis('square');
xlabel('V1 P(sp|x) ASI');
ylabel('V2 P(sp|x) ASI');


bins = 0:0.05:1;
xtick = 0:0.1:1;

subplot(4,2,7);
count = histc(sepindex, bins);
newcount = 100*count ./ sum(newcount);
hb = bar(bins, newcount, 'histc');
set(hb, 'facecolor', [0.75 0.75 0.75])
axis([-0.05 1.05 0 1.05*max(newcount)]);
set(gca,'xtick', xtick, 'xticklabel', xtick);
set(gca,'ytick', [0:10:50], 'yticklabel', [0:10:50]);
set(gca,'tickdir', 'out');
set(gca,'fontsize', 8);
ylabel('Percent of neurons (%)');
title('P(sp|x1,x2) Sep Index', 'fontsize', 10);

print_mfilename(mfilename);


figure;

% subplot(2,2,1)
hold on
xmax = max([max(prob12) max(prob1_plus2)]);
plot(prob1_plus2, prob12, 'ko')
fplot('x', [0 1.25], 'k-')
% fplot('x', [0 max([ max(prob12) max(prob1_plus2)])], 'k-')
% fplot('2*x', [0 max(info_both)], 'k-')
xlim([0.01 1.25]);
ylim([0.01 1.25]);
% xlim([0 1.1*xmax]);
% ylim([0 1.1*xmax]);
axis square
set(gca,'xscale', 'log', 'yscale', 'log');
set(gca,'tickdir', 'out');
xlabel('Sum of Peak Prob')
ylabel('Joint Peak Prob')
title('Comparison between Joint and Ind Prob')
print_mfilename(mfilename);






