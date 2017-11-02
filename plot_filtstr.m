function plot_filtstr(filtstr)
%plot_filtstr STA, MID1, and MID2 filters
%
% plot_filtstr(filtstr)
% ---------------------------------------------------------
%
% filtstr : struct array holding filters from mid analysis.
%
% For use with exportfig.m, the following option give good results:
%
% opts = struct('FontMode','fixed', 'color', 'rgb', 'height', 6.2, 'width',
% 4.8);
%
%   caa 1/23/10

mapscheme = 'rdbu';

if ( nargin ~= 1 )
   error('You need 1 input arg.');
end

gaussian = fspecial('gaussian', [2 2], 1);


time = 1000 .* filtstr(1).time; % change to ms
freq = filtstr(1).freq / 1000;

ytick = [1 floor(length(freq)/2) length(freq)];
yticklabel = round(freq(ytick)*10)/10;

xtick = [1 floor(length(time)/2)+1 length(time)];
xticklabel = round( time(xtick) );

close all;

ntot = 1;
len = length(filtstr);
fignum = 1;
figure;

for i = 1:length(filtstr)

   location = filtstr(i).location;
   unit = filtstr(i).unit;
   site = filtstr(i).site;
   exp = filtstr(i).exp;
   stim = filtstr(i).stim;

   sta = filtstr(i).v_sta;
	sta = imfilter(sta,gaussian);

%    sta = fliplr(sta);
   v1 = filtstr(i).v1;
	v1 = imfilter(v1,gaussian);
%    v1 = fliplr(v1);
   v2 = filtstr(i).v2;
   v2 = imfilter(v2,gaussian);
%    v2 = fliplr(v2);

   subplot(5,3,fignum);
   imagesc( sta );
%    imagesc(time, freq, sta );
   axis xy;
   minmin = min(min(sta));
   maxmax = max(max(sta));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    cmap = color_brewer_colormap('rdbu');
   cmap = cschemes(mapscheme,15);
%    colormap(jet);
%    colormap(cmap);
   ylabel(sprintf('%.0f - %.0f - %.0f', location, site, unit));
   if ( fignum == 1 )
      title('STA');
   end
   if ( fignum ~= 13 )
      set(gca, 'xticklabel', []);
   end


   subplot(5,3,fignum+1);
   imagesc( v1 );
   axis xy;
   minmin = min(min(v1));
   maxmax = max(max(v1));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
%    colormap(jet);
%    colorbar;
   colormap(cmap);
   if ( fignum == 1 )
      title('MID1');
   end
   if ( fignum ~= 13 )
      set(gca, 'xticklabel', []);
   end


   subplot(5,3,fignum+2);
   imagesc( v2 );
   axis xy;
   minmin = min(min(v2));
   maxmax = max(max(v2));
   boundary = max([abs(minmin) abs(maxmax)]);
   set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   set(gca,'ytick', ytick, 'yticklabel', yticklabel);
   set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
   set(gca, 'yticklabel', []);
   cmap = cschemes(mapscheme,15);
%    colormap(jet);
%    colorbar;
   colormap(cmap);
   if ( fignum == 1 )
      title('MID2');
   end
   if ( fignum ~= 13 )
      set(gca, 'xticklabel', []);
   end


   if ( mod(fignum,13) )
      fignum = fignum + 3;
      figuretitle = 0;
   else
      orient tall;
      print_mfilename(mfilename);
%       subplotspace('vertical', 20);
      set(gcf,'position',[50 100 450 750]);
      figuretitle = 1;

      if ( i ~= len )
         fignum = 1;
         figure;
      end
   end

   ntot = ntot + 1;

end

if ( figuretitle == 0 )
   orient tall;
%    suptitle(sprintf('ICC Data'));
   print_mfilename(mfilename);
   set(gcf,'position',[50 100 450 750]);
end

   suptitle(sprintf('%s site%.0f %s', exp, site, stim));


return;


























