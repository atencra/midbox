function midptparams = get_mid_puretone_params(filtstr)
%get_mid_puretone_params   CF, BW, Latency from MIDs. 
%
%   ptparams = get_mid_puretone_params(filtstr)
%
%   filtstr is a struct array holding the receptive field data
%
%   caa 11/20/10


midptparams = filtstr;
fields = {'time', 'freq', 'v_sta', 'mtx_sta', 'v1', 'mtx_v1', 'v2', 'mtx_v2'};
midptparams = rmfield(midptparams, fields);

close all;

figure;

for i = 1:length(filtstr)

    % Assign the initial data parameters 
        % These have already been assigned earlier
    
    
    % plot the filters
   time = 1000 .* filtstr(i).time; % change to ms
   freq = filtstr(i).freq / 1000;

   ytick = [1 floor(length(freq)/2) length(freq)];
   yticklabel = round(freq(ytick)*100)/100;
   
   xtick = [1 floor(length(time)/2)+1 length(time)];
   xticklabel = round( time(xtick) );

   sta = filtstr(i).v_sta;
   sta = fliplr(sta);
   v1 = filtstr(i).v1;
   v1 = fliplr(v1);


%    h = fspecial('gaussian',5,5);
%    rfsig = imfilter(rfsig, h);

   clf;
   subplot(1,2,1);
   hold on;
   imagesc( sta );
% %    imagesc(time, freq, sta );
   axis xy;
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    cmap = color_brewer_colormap('rdbu');
   colormap(jet);
%    colormap(cmap);
   title(sprintf('STA: %.0f of %.0f', i, length(filtstr) ));
   print_mfilename(mfilename);



   subplot(1,2,2);
   hold on;
   imagesc( v1 );
% %    imagesc(time, freq, sta );
   axis xy;
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    cmap = color_brewer_colormap('rdbu');
   colormap(jet);
%    colormap(cmap);
   title(sprintf('MID1' ));

   process = input('\nProcess STA and MID1? 1=yes, 0=no: ');

   if ( process )

      % Get frequency parameters
      sta(sta<0) = 0;
      temp = sum(sta,2);
      xi = linspace(min(freq),max(freq),1000);
      yi = interp1(freq, temp, xi, 'spline');
      maxmax = max(temp);
      thr = 0.25 * maxmax;

      clf;
      hold on;
      plot(freq,temp,'ko');
      plot(xi, yi, 'r-');

      [x,y] = ginput(2);
      index = find( xi > min(x) & xi < max(x) );
      xtemp = xi(index);
      ytemp = yi(index);
      index = find(ytemp > thr);
      flow = xtemp(min(index));
      fhi = xtemp(max(index));
      index = find( ytemp == max(ytemp) );
      cf = xtemp(index);
      bw = fhi - flow;
      q = cf / bw;

      plot([cf cf], [0 1.25 * maxmax], 'b-');
      plot([flow fhi], [thr thr], 'b-');
      xlabel('Frequency [kHz]');
pause

      % Get latency
      temp = sum(sta,1);
      xi = linspace(min(time),max(time),1000);
      yi = interp1(time, temp, xi, 'spline');
      maxmax = max(temp);

      clf;
      hold on;
      plot(time,temp,'ko');
      plot(xi, yi, 'r-');

      [x,y] = ginput(2);
      index = find( xi > min(x) & xi < max(x) );
      xtemp = xi(index);
      ytemp = yi(index);
      index = find( ytemp == max(ytemp) );
      latency = xtemp(index);
      plot([latency latency], [0 1.25 * maxmax], 'b-');
      xlabel('Time [ms]');
pause

        midptparams(i).sta_latency = latency;
        midptparams(i).sta_cf = cf;
        midptparams(i).sta_flow = flow;
        midptparams(i).sta_fhi = fhi;
        midptparams(i).sta_bw = bw;
        midptparams(i).sta_q = q;

        clear latency cf flow fhi bw q

      % MID1: Get frequency parameters
      v1(v1<0) = 0;
      temp = sum(v1,2);
      xi = linspace(min(freq),max(freq),1000);
      yi = interp1(freq, temp, xi, 'spline');
      maxmax = max(temp);
      thr = 0.25*maxmax;

      clf;
      hold on;
      plot(freq,temp,'ko');
      plot(xi, yi, 'r-');

      [x,y] = ginput(2);
      index = find( xi > min(x) & xi < max(x) );
      xtemp = xi(index);
      ytemp = yi(index);
      index = find(ytemp > thr);
      flow = xtemp(min(index));
      fhi = xtemp(max(index));
      index = find( ytemp == max(ytemp) );
      cf = xtemp(index);
      bw = fhi - flow;
      q = cf / bw;

      plot([cf cf], [0 1.25 * maxmax], 'b-');
      plot([flow fhi], [thr thr], 'b-');
      xlabel('Frequency [kHz]');
pause

      % MID1: Get latency
      temp = sum(v1,1);
      xi = linspace(min(time),max(time),1000);
      yi = interp1(time, temp, xi, 'spline');
      maxmax = max(temp);

      clf;
      hold on;
      plot(time,temp,'ko');
      plot(xi, yi, 'r-');

      [x,y] = ginput(2);
      index = find( xi > min(x) & xi < max(x) );
      xtemp = xi(index);
      ytemp = yi(index);
      index = find( ytemp == max(ytemp) );
      latency = xtemp(index);
      plot([latency latency], [0 1.25 * maxmax], 'b-');
      xlabel('Time [ms]');
pause

        midptparams(i).mid1_latency = latency;
        midptparams(i).mid1_cf = cf;
        midptparams(i).mid1_flow = flow;
        midptparams(i).mid1_fhi = fhi;
        midptparams(i).mid1_bw = bw;
        midptparams(i).mid1_q = q;

        clear latency cf flow fhi bw q

   end 

end


return;
