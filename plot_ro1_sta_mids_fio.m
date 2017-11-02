function plot_ro1_sta_mids_fio(mid, spk, trigger)
%plot_sta_pspx_mid1_pspx_mid2_pspx_pspx1x2 - Plots tuning curve data
%   obtained using the Michigan 16 channel silicon probe.
%
%

% 2004-1-14 site 16: 21, 19, 12, 7, 
% 2004-1-14 site 16: 2, 6, 8, 13
% 2003-3-5 site32: 10, 11, 15, 16
% 2002-8-27 site17: 7, 8, 9, 10, 17, 19 

ntot = 1;
fignum = 1;
figure;
nr = 6;
max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];

numplots = ceil(length(mid)/5);

for i = 1:length(mid)

   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
%    sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   fx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   pspx_sta = fx_sta;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxsp_sta = fx_sta ./ pspike .* px_sta;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
%    v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   pspx_v1 = fx_v1;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxsp_v1 = fx_v1 ./ pspike .* px_v1;

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
%    v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
   max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
   max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];

   % Plot the STA
   %------------------------------
   subplot(nr,7,(fignum-1)*7+1);
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
%    imagesc(smoothmat(sta));
   plot_strf_symmetric_colormap(smoothmat(sta));
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   tickpref;
%    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( isfield(mid,'location') )
      ylabel(sprintf('%.0f : %.0f - %.0f', i, mid(i).location, mid(i).unit));
   else
      ylabel(sprintf('#%.0f %s %.0f\nc%.0f m%.0f %.0fum', ...
         i, mid(i).exp, mid(i).site, mid(i).chan, mid(i).model, mid(i).position));
   end
   if ( fignum == 1 )
      title('STA', 'fontsize', 8);
   end;

   % Plot the STA nonlinearity
   %------------------------------
   subplot(nr,7,(fignum-1)*7+2);
   hold on;
   plot(x_sta, fx_sta, 'ko-', 'markerfacecolor', 'k', 'markersize', 2, 'linewidth', 1);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f',1.1*max(fx_sta)));
   ylim([0 1.1*max(fx_sta)]);
   set(gca,'ytick', [0 max(fx_sta)/2 max(fx_sta)], 'yticklabel', [0 maxfx/2 maxfx]);
   tickpref;
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('STA P(spk|proj)', 'fontsize', 8);
   end;


   % Plot the MID1
   %------------------------------
   subplot(nr,7,(fignum-1)*7+3);
   v1(:,end) = 0.9 * (rand(length(v1(:,end)),1)-0.5);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
%    imagesc(smoothmat(v1));
   plot_strf_symmetric_colormap(smoothmat(v1));
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   tickpref;
%    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1', 'fontsize', 8);
   end;

   % Plot the MID1 nonlinearity
   %------------------------------
   subplot(nr,7,(fignum-1)*7+4);
   hold on;
   plot(x_v1, fx_v1, 'ko-', 'markerfacecolor', 'k', 'markersize', 2, 'linewidth', 1);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f',1.1*max(fx_v1)));
   ylim([0 1.1*max(fx_v1)]);
   set(gca,'ytick', [0 max(fx_v1)/2 max(fx_v1)], 'yticklabel', [0 maxfx/2 maxfx]);
   tickpref;
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1 P(spk|proj)', 'fontsize', 8);
   end;


   % Plot the MID2
   %------------------------------
   subplot(nr,7,(fignum-1)*7+5);
   v2(:,end) = 0.25 * (rand(length(v2(:,end)),1)-0.5);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(smoothmat(v2));
   plot_strf_symmetric_colormap(smoothmat(v2));
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   tickpref;
%    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_2', 'fontsize', 8);
   end;

   % Plot the MID2 nonlinearity
   %------------------------------
   subplot(nr,7,(fignum-1)*7+6);
   hold on;
   plot(x_v2, fx_v2, 'ko-', 'markerfacecolor', 'k', 'markersize', 2, 'linewidth', 1);
   plot([min(x_sta) max(x_sta)], [pspike pspike], 'k:');
   xlim([-7.5 7.5]);
   set(gca,'xtick', [-7.5 -5 -2.5 0 2.5 5 7.5], 'xticklabel', [-7.5 -5 -2.5 0 2.5 5 7.5]);
   maxfx = str2num(sprintf('%.3f',1.1*max(fx_v2)));
   ylim([0 1.1*max(fx_v2)]);
   set(gca,'ytick', [0 max(fx_v2)/2 max(fx_v2)], 'yticklabel', [0 maxfx/2 maxfx]);
   tickpref;
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_2 P(spk|proj)', 'fontsize', 8);
   end;


   % Plot the MID1 and MID2 nonlinearity
   %------------------------------------
   subplot(nr,7,(fignum-1)*7+7);
   %imagesc(x12, y12, fx_v1v2);
   %get(hi)
   %set(hi, 'zscale', 'log')
   h = imagesc(x12, y12, log10(fx_v1v2+0.001));
%    get(h,'xdata')
%    pause
   %cm = colormap(gray);
   %colormap(gca, flipud(cm));
   set(gca,'clim', [-3 0]);
   tickpref;
   if ( fignum == 5 )
      xlabel('MID1 Proj','fontsize',8);
      ylabel('MID2 Proj', 'fontsize',8);
   end
   axis('xy');
   %axis('square');
%    hc = colorbar;
   set(gca,'fontsize', 8);
%    set(hc,'ytick', [-3 -2 -1 0], 'yticklabel', [0.001 0.01 0.1 1]);
   if ( fignum == 1 )
      title('P(spk|x1,x2)', 'fontsize', 8);
   end


   if ( mod(fignum,6) )
      fignum = fignum + 1;
   else
      orient landscape;
      print_mfilename(mfilename);
      set(gcf,'position',[50 100 700 800]);
      fignum = 1;
      if ( i ~= numplots ), figure; end;
   end

end % (for i)



close all;

ntot = 1;
fignum = 1;
figure;
nr = 6;
nc = 2;
max_min_sta = [];
max_min_v1 = [];
max_min_v2 = [];

numplots = ceil(length(mid)/5);

for i = 1:length(mid)

   tbins = mid(i).tbins;
   fbins = mid(i).fbins;

   taxis = (0:tbins-1) * 5;
   faxis = logspace(log10(500), log10(20000), fbins);

   pspike = mid(i).rpx1pxpxt_sta.mean_firing;

   sta = mid(i).rpsta.filter;
   sta = reshape(sta, fbins, tbins);
%    sta(abs(sta)<0.75) = 0;
   x_sta = mid(i).rpx1pxpxt_sta.x;
   fx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   pspx_sta = fx_sta;
   px_sta = mid(i).rpx1pxpxt_sta.px_mean;
   pxsp_sta = fx_sta ./ pspike .* px_sta;


   v1 = mid(i).rpdtest2_v1.filter;
   v1 = reshape(v1, fbins, tbins);
%    v1(abs(v1)<0.75) = 0;
   x_v1 = mid(i).rpdx1x2px_pxt_2.x1;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   pspx_v1 = fx_v1;
   px_v1 = mid(i).rpdx1x2px_pxt_2.px1_mean;
   pxsp_v1 = fx_v1 ./ pspike .* px_v1;

   v2 = mid(i).rpdtest2_v2.filter;
   v2 = reshape(v2, fbins, tbins);
%    v2(abs(v2)<0.75) = 0;
   x_v2 = mid(i).rpdx1x2px_pxt_2.x2;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;
   pspx_v2 = fx_v2;
   px_v2 = mid(i).rpdx1x2px_pxt_2.px2_mean;
   pxsp_v2 = fx_v2 ./ pspike .* px_v2;

   x12 = mid(i).rpdx1x2px_pxt_2.x12;
   y12 = mid(i).rpdx1x2px_pxt_2.y12;
   fx_v1v2 = mid(i).rpdx1x2px_pxt_2.ior12_mean;

   max_min_sta = [max_min_sta; min(x_sta) max(x_sta)];
   max_min_v1 = [max_min_v1; min(x_v1) max(x_v1)];
   max_min_v2 = [max_min_v2; min(x_v2) max(x_v2)];


   % Plot the MID1
   %------------------------------
   subplot(nr,nc,(fignum-1)*nc+1);
   v1(:,end) = 0.9 * (rand(length(v1(:,end)),1)-0.5);
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(smoothmat(v1));
%    plot_strf_symmetric_colormap(smoothmat(v1));
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   tickpref;
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_1', 'fontsize', 8);
   end;


   % Plot the MID2
   %------------------------------
   subplot(nr,nc,(fignum-1)*nc+2);
   v2(:,end) = 0.25 * (rand(length(v2(:,end)),1)-0.5);
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   imagesc(smoothmat(v2));
%    plot_strf_symmetric_colormap(smoothmat(v2));
   axis('xy');
   %colorbar;
   %set(gca,'ydir', 'normal');
   tickpref;
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'fontsize', 8);
   if ( fignum == 1 )
      title('MID_2', 'fontsize', 8);
   end;


   if ( mod(fignum,6) )
      fignum = fignum + 1;
   else
      orient landscape;
      print_mfilename(mfilename);
      set(gcf,'position',[50 100 700 800]);
      fignum = 1;
      if ( i ~= numplots ), figure; end;
   end

end % (for i)

















if ( nargin == 3 )
   % Assign beginning, ending, and duration of stimulus
   onset = trigger(1); % Find beginning of stimulus, in sample number
   offset = trigger(end); % End of stimulus
   fsspk = spk(1).fs;
   stimdur = 900000;
   fsd = 1000; % sampling rate for spike train; 1000 = 1 ms resolution
   cc = calc_cross_correlation(spk, fsspk, fsd, stimdur, onset, offset);
end

return;





function cc = calc_cross_correlation(spk, fsspk, fsd, stimdur, onset, offset)
% Calculate the cross-correlations

% Define a struct to hold the results
cc = struct(...
'exp',        [], ...
'site',       [], ...
'chan',       [], ...
'model',      [], ...
'depth',      [], ...
'position',   [], ...
'stim',       [], ...
'atten',      [], ...
'fsd',        [], ...
'stimdur',    [], ...
'na',         [], ...
'nb',         [], ...
'delay',      [], ...
'rab',        []);


cmb = nchoosek(1:length(spk),2); % determine all possible pairwise 
                                 % combinations for every recording
[nr,nc] = size(cmb);

for i = 1:nr

   c = cmb(i,:);

   spk1 = spk(c(1));
   spk2 = spk(c(2));

   chan1 = spk1.chan;
   model1 = spk1.model;

   chan2 = spk2.chan;
   model2 = spk2.model;
      
   pos1 = spk1.position;
   pos2 = spk2.position;


   fprintf('Processing channel %.0f, model %.0f and channel %.0f, model %.0f\n', ...
           chan1, model1(1), chan2, model2(1));

   % Assign some basic parameters
   cc(i).exp = spk1.exp;
   cc(i).site = spk1.site;
   cc(i).chan = {chan1 chan2};
   cc(i).model = {model1 model2};
   cc(i).depth = spk1.depth;
   cc(i).position = {pos1 pos2};
   cc(i).stim = spk1.stim;
   cc(i).atten = spk1.atten;
   cc(i).fsd = fsd;
   cc(i).stimdur = stimdur;

   duration = (offset - onset); % duration in sample number

   spiketimes1 = spk1.spiketimes; % spikes in units of milliseconds
   spiketimes2 = spk2.spiketimes;

   maxspiketime = max([ max(spiketimes1) max(spiketimes2) ] );
   spiketimes1 = [spiketimes1 maxspiketime];
   spiketimes2 = [spiketimes2 maxspiketime];


   spet1 = spiketimes1 ./ 1000 * fsspk; % now in sample number
   spet1 = spet1 - onset;
   spet1 = spet1( spet1 > 0 & spet1 < duration );

   spet2 = spiketimes2 ./ 1000 * fsspk;
   spet2 = spet2 - onset;
   spet2 = spet2( spet2 > 0 & spet2 < duration );

   n1 = length(spet1);
   n2 = length(spet2);

   train1 = spet2train(spet1, fsspk, fsd); 
   train2 = spet2train(spet2, fsspk, fsd); 

   fsd = 1000; % Hz = 1 ms resolution
   dt = 1 / fsd * 1000; % in ms


   % Compute cab = E[xa(t+u)xb(t)]
   [r12, delay] = xcorr(train1, train2, 100); % make max lag = 100 * 1 ms = 100 ms
   r12 = real(r12);
   r12(find(r12 < 0)) = 0;
   delay = delay .* dt; % delay in ms, not bin number


   % save the data to the struct array
   ccstr(i).fsd = fsd;
   ccstr(i).dt = dt;
   ccstr(i).n1 = n1; % the number of spikes
   ccstr(i).n2 = n2; % the number of spikes
   ccstr(i).delay = delay;
   ccstr(i).r12 = r12;

   [q12, conf_limit, stimdur] = get_cross_covariance(r12, n1, n2, dt);

   [pd, cent, ca, hw, sigfeature] = ...
      get_pairs_cross_covariance_params(delay, q12, r12, conf_limit);

   rho = corrstrength(delay, r12, n1, n2, stimdur); % correlation coefficient


   figure;
   subplot(2,1,1);
   bar(delay, r12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
   tickpref;
   xlim([min(delay) max(delay)]);
   title(sprintf('c%.0f m%.0f %.0fum / c%.0f m%.0f %.0fum', ...
      chan1, model1, pos1, chan2, model2, pos2));

   subplot(2,1,2);
   upper95qab = conf_limit;
   lower95qab = -conf_limit;
   hold on;
   ymin = min([min(q12) lower95qab]);
   ymax = max([max(q12) upper95qab]);
   range = ymax-ymin;
   bar(delay, q12, 'k'); %, 'markerfacecolor', 'k', 'markersize', 2);
   plot([0 0], [ymin-0.1*range ymax+0.1*range], 'k-');
   plot([min(delay) max(delay)], [0 0], 'k-');
   plot([min(delay) max(delay)], [upper95qab upper95qab], 'r-');
   plot([min(delay) max(delay)], [lower95qab lower95qab], 'r-');
   tickpref;
   xlim([min(delay) max(delay)]);


   pause(1);

end % (for i)


return;





