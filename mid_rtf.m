function [tmf, xmf, rtf] = mid_rtf(v)
% MID_RTF - Calculate the ripple transfer function for an MID
%
% [tmf, xmf, rtf] = mid_rtf(v)
%
% v : mid filter. May be either MID1 or MID2. This function assumes that
% the time bins were 5 ms and the spectral resolution was 
% 6 carriers / octave.
%
% tmf : temporal modulation frequency axis (cycles / s).
% xmf : spectral modulaiton frequency axis (cycles / octave).
%
% rtf : ripple transfer function. A matrix where spectral modulation
%    varies along the rows and temporal modulation varies along the
%    columns of the matrix.
%
% caa



%    dt = diff(t);
%    dt = dt(1); % strf temporal resolution
%    dx = diff(x);
%    dx = dx(1); % strf spectral resolution

   dt = 0.005;
   dx = 1/6;

   maxtmf = ceil( 1 / dt );
   maxsmf = ceil( 1 / dx );

   ntbins = maxtmf; % this will give us 1 Hz tmtf resolution
   nfbins = ceil(maxsmf / 0.1); % make resolution 0.1 cycles / octave

   dtfreq = maxtmf / ntbins; %size(rtftemp,2); % fft temporal frequency resolution
   dxfreq = maxsmf / nfbins; %size(rtftemp,1); % fft spectral frequency resolution



   % get the ripple transfer function/s
   rfk = fft2(v, nfbins, ntbins);
   rtftemp = fftshift(abs(rfk));

   % Get the tm frequency vector - will be in cycles per second
   if ( mod(size(rtftemp,2),2) )
      tmf = [-(size(rtftemp,2)-1)/2:(size(rtftemp,2)-1)/2]*dtfreq;
   else
      tmf = [-size(rtftemp,2)/2:(size(rtftemp,2)/2-1)]*dtfreq;
   end

   itmf0 = find(tmf==0);
   itmfp40 = find(tmf<=40);
   ditmf = max(itmfp40)-itmf0;
   itmf = [itmf0-ditmf:itmf0+ditmf];
   tmf = tmf(itmf);

   % Get the sm frequency vector - will be in cycles per octave
   if ( mod(size(rtftemp,1),2) )
      xmf = [-(size(rtftemp,1)-1)/2:(size(rtftemp,1)-1)/2]*dxfreq;
   else
      xmf = [-size(rtftemp,1)/2:(size(rtftemp,1)/2-1)]*dxfreq;
   end

   ixmf0 = find(xmf == 0);
   ixmf4 = find(xmf <= 4);
   ixmf = ixmf0:max(ixmf4);
   xmf = xmf(ixmf);

   rtf = rtftemp(ixmf, itmf);


% figure;
% subplot(2,1,1);
% imagesc(v)
% axis('xy');
% title('Filter');
% 
% subplot(2,1,2);
% imagesc(rtf);
% title('RTF');
% axis('xy');
% pause

return;

