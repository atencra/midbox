function get_pspkx_data(fio_params, mid)
%
% Figure 4 can be save using the following command:
%
% exportfig(gcf,'sta_v1_v2_asymmetry.eps', 'color', 'gray', 'height', 8, 'width', 6, 'fontmode', 'fixed', 'fontsize', 8);

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
   simindex(i) = c(1,2);

   fx1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   fx2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   fx1x2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;
   [u,s,v] = svd(fx1x2);
   singvals = sum(s);
   eigvals = singvals .^ 2;
   sepindex(i) = 1 - eigvals(1) / (sum(eigvals)+eps);

end % (for)


mnsta = mean(sta_fx_asi);
stdsta = std(sta_fx_asi);
numsta = length(sta_fx_asi);

mnv1 = mean(v1_fx_asi);
stdv1 = std(v1_fx_asi);
numv1 = length(v1_fx_asi);

mnv2 = mean(v2_fx_asi);
stdv2 = std(v2_fx_asi);
numv2 = length(v2_fx_asi);

mnsim = mean(simindex);
stdsim = std(simindex);
numsim = length(simindex);

mnsep = mean(sepindex);
stdsep = std(sepindex);
numsep = length(sepindex);


fprintf('\nSTA: mean asi=%.3f, std = %.3f, N=%.0f\n',  mnsta, stdsta, numsta);

fprintf('MID1: mean asi=%.3f, std = %.3f, N=%.0f\n',  mnv1, stdv1, numv1);

fprintf('MID2: mean asi=%.3f, std = %.3f, N=%.0f\n',  mnv2, stdv2, numv2);

fprintf('SimIndex: mean = %.3f, std = %.3f, N=%.0f\n',  mnsim, stdsim, numsim); 

fprintf('SepIndex: mean = %.3f, std = %.3f, N=%.0f\n\n',  mnsep, stdsep, numsep); 



