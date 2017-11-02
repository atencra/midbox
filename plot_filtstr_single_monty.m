function plot_filtstr_single_monty(filtstr)
%plot_filtstr_single_monty - STA, MID1, and MID2 filters for one neuron
%
% plot_filtstr_single_monty(filtstr)
% ---------------------------------------------------------
%
% filtstr : struct array holding filters from mid analysis.
%
% If length(filtstr) > 1, then only the first element of filtstr is
% plotted.
%
%   caa 1/23/10

if ( nargin ~= 1 )
   error('You need 1 input arg: filtstr');
end

if ( length(filtstr) ~= 1 )
   filtstr = filtstr(1);
end


time = 1000 .* filtstr(1).time; % change to ms
freq = filtstr(1).freq / 1000;

ytick = [1 13 25];
yticklabel = round(freq(ytick)*10)/10;

xtick = [1 11 20];
xticklabel = round( time(xtick) );

close all;

figure;

i = 1;

location = filtstr(i).location;
unit = filtstr(i).unit;

sta = filtstr(i).v_sta;
sta = fliplr(sta);

v1 = filtstr(i).v1;
v1 = fliplr(v1);

if ( isfield(filtstr(i), 'v2') )
   v2 = filtstr(i).v2;
   v2 = fliplr(v2);
else
	v2 = [];
end


if ( ~isempty(v2) )

   subplot(2,2,1);
   imagesc( sta );
   axis xy;
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
% get(gca)
%    hc = colorbar;
%    set(hc, 'position', [0.93 0.72 0.02 0.2]);
%    set(hc, 'position', [0.48 0.8 0.02 0.1]);
   ylabel('Frequency (kHz)');
   title('STA');


   subplot(2,2,3);
   imagesc( v1 );
   axis xy;
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    hc = colorbar;
%    set(hc, 'position', [0.93 0.42 0.02 0.2]);
	xlabel('Time [ms]');
   ylabel('Frequency (kHz)');
   title('MID1');


   subplot(2,2,4);
   imagesc( v2 );
   axis xy;
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    cmap = color_brewer_colormap('rdbu');
%    hc = colorbar;
%    colormap(cmap);
	colormap(jet);
% get(hc)
%    set(hc, 'tickdir', 'out');   
%    set(hc, 'position', [0.93 0.12 0.02 0.2]);
	xlabel('Time [ms]');
   title('MID2');

else

   subplot(2,1,1);
   imagesc( sta );
   axis xy;
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    hc = colorbar;
%    set(hc, 'position', [0.93 0.72 0.02 0.2]);
   ylabel(sprintf('%.0f - %.0f', location, unit));
   title('STA');
   set(gca, 'xticklabel', []);


   subplot(2,1,2);
   imagesc( v1 );
   axis xy;
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    hc = colorbar;
%    set(hc, 'position', [0.93 0.42 0.02 0.2]);
   ylabel('Frequency (kHz)');
   title('MID1');
   set(gca, 'xticklabel', []);

end

   orient tall;
%    print_mfilename(mfilename);
   set(gcf,'position',[560   100   750   800]);

if ( filtstr.location == 23 )
	unit = 'CATICC07PQA1Block10Tetrode2u3';
	suptitle(unit);
elseif ( filtstr.location == 35 )
	unit = 'CATICC07PFK4 3Block35Tetrode3u5';
	suptitle(unit);
elseif ( filtstr.location == 43 )
	unit = 'CATICC07PFK4Block11Tetrode4u3';
	suptitle(unit);
elseif ( filtstr.location == 44 )
	unit = 'CATICC08JKP4cBlock33Tetrode4u4';
	suptitle(unit);
end



return;


























