
function mid_fsi_from_proj(proj)
   

% % Next we want p(projection onto mid1 & mid2 | random spikes), so use the original filter but this time
% % with 50,000 random spikes 
% [x1x2rand, spindexrand] = get_mid1_mid2_projection(mid1, mid2, sprfile, 0.005, 0.095, spetrand, trigger, fs, spl, mdb);
% 
% 
% % First we want to find the p(proj onto mid1 & mid2 | spike), so use the original spike
% % train and the original filter
% [x1x2spk, spindex] = get_mid1_mid2_projection(mid1, mid2, sprfile, 0.005, 0.095, spet, trigger, fs, spl, mdb);


x1x2rand = proj(1).x1x2rand_mean;

x1x2spk = proj(1).x1x2spk_mean;


x1rand = x1x2rand(:,1);
x2rand = x1x2rand(:,2);

std1 = std(x1rand);
std2 = std(x2rand);

x1spk = x1x2spk(:,1);
x2spk = x1x2spk(:,2);

x1rand = x1rand ./ std1;
x2rand = x2rand ./ std2;

x1spk = x1spk ./ std1;
x2spk = x2spk ./ std2;

close all;

subplot(2,3,1);
[n1, x1] = hist(x1rand, 50);
n1 = n1 ./ sum(n1);
hb = bar(x1, n1, 1);
set(hb, 'facecolor', 0.6*[1 1 1]);
title('P(x1)');

subplot(2,3,2);
[n1s, x1s] = hist(x1spk, 50);
n1s = n1s ./ sum(n1s);
hb = bar(x1s, n1s, 1);
set(hb, 'facecolor', 0.6*[1 1 1]);
title('P(x1|spk)');

subplot(2,3,3);
hold on;
plot(x1, n1, 'k-');
plot(x1s, n1s, 'r-');


subplot(2,3,4);
[n2, x2] = hist(x2rand, 50);
n2 = n2 ./ sum(n2);
hb = bar(x2, n2, 1);
set(hb, 'facecolor', 0.6*[1 1 1]);
title('P(x2)');

subplot(2,3,5);
[n2s, x2s] = hist(x2spk, 50);
n2s = n2s ./ sum(n2s);
hb = bar(x2s, n2s, 1);
set(hb, 'facecolor', 0.6*[1 1 1]);
title('P(x2|spk)');

subplot(2,3,6);
hold on;
plot(x2, n2, 'k-');
plot(x2s, n2s, 'r-');


set(gcf,'position', [250 160 900 700]);
print_mfilename(mfilename);




[x1cdf, y1cdf] = cum_dist_func(x1spk);
[x1rcdf, y1rcdf] = cum_dist_func(x1rand);

[y1cdf, x1cdf] = cdfcalc(x1spk);
[y1rcdf, x1rcdf] = cdfcalc(x1rand);


min([length(x1cdf) length(y1cdf)]);

x1cdf = x1cdf( 1:min([length(x1cdf) length(y1cdf)]) );
y1cdf = y1cdf( 1:min([length(x1cdf) length(y1cdf)]) );

x1rcdf = x1rcdf( 1:min([length(x1rcdf) length(y1rcdf)]) );
y1rcdf = y1rcdf( 1:min([length(x1rcdf) length(y1rcdf)]) );


figure;

hold on;
plot(x1rcdf, y1rcdf, 'k-');
plot(x1cdf, y1cdf, 'r-');


[fsi1, fsi2, fsi3] = fsi_caa(x1spk, x1rand)

set(gcf,'position', [250 160 900 700]);
print_mfilename(mfilename);


return;



function [xcdf, ycdf] = cum_dist_func(data)

% Now we'll plot the CDF of the data

[yy, xx, nn] = cdfcalc(data);
k = length(xx);
nn = reshape(repmat(1:k, 2, 1), 2*k, 1);
xcdf = [-Inf; xx(nn); Inf];
ycdf = [0; 0; yy(1+nn)];

return;


function [FSI1,FSI2,FSI3] = fsi_caa(p,pr)


mxp = max( [max(p) max(abs(p))] );

mxpr = max( [max(pr) max(abs(pr))] );


%Generating Cummulative Distribution Function
[N,X] = hist(p, linspace(-mxp, mxp, 64) );
N = N ./ sum(N);
[Nr,Xr] = hist(pr, linspace(-mxpr, mxpr, 64) );
Nr = Nr ./ sum(Nr);
CDF = (intfft(N)-min(intfft(N))) / (max(intfft(N))-min(intfft(N)));
CDFr = (intfft(Nr)-min(intfft(Nr))) / (max(intfft(Nr))-min(intfft(Nr)));

%Generating Cummulative Distribution Function - Ideal Feature Detector 
CDFi = zeros(1,64);
CDFi(64) = 1;

%Feature Selectivity Index
FSI1 = sum(CDFr-CDF) / sum(CDFr-CDFi);
FSI2 = (mean(p)-mean(pr)) / (1-mean(pr));
FSI3 = (median(p)-median(pr)) / (1-median(pr));
















