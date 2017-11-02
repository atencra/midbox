function plot_filter_files(location, global_cell)
% plot_filter_files - display sta, mids from directory files
%
% plot_filter_files(location, global_cell)
%
% location : value such as 707, 708, 2171, etc.
%
% global_cell : cell number used in analysis, example, 1, 3, 5, etc.
% 
% %Ncells = 27; % Number of cells in the penetration
% 
% pathname = 'C:\MATLAB65\work\tatyana\Filters\site532';
% location = 532; % one of the sites I gave Tatyana
% nlags = 20; % number of time delay bins in receptive fields
% 
% if (location<600)
%     Nv = 25; % Number of frequency bands
% else
%     Nv = 30; % Number of frequency bands
% end
% 
% cy = 2; % used to identify file name for MID1 and MID2
% optlevel = 2;
% 
% %for curr_cell=[2,3,4,6:Ncells]
% 
% for curr_cell=[1:10]
%    hfig = plot_one_cell_tatyana(pathname, location, curr_cell, nlags, Nv, cy, optlevel);
%    %print(hfig,'-dpsc',sprintf('%u_%u_dir2', location, curr_cell))
%    pause
% end
%
%
% % Do we have four files for each cell because Tatyana
% % used 4 different random starting points for the 
% % optimization routine?
% %
% % rpsta_532_10_1x25x20_1_1.dat
% % rpsta_532_10_1x25x20_1_2.dat
% % rpsta_532_10_1x25x20_1_3.dat
% % rpsta_532_10_1x25x20_1_4.dat
% % 
% % rpx1pxpxt_sta_532_10_1x25x20_1_1.dat
% % rpx1pxpxt_sta_532_10_1x25x20_1_2.dat
% % rpx1pxpxt_sta_532_10_1x25x20_1_3.dat
% % rpx1pxpxt_sta_532_10_1x25x20_1_4.dat
%
% return;
% 
% 
% 
% 
% 
% function hfig = plot_one_cell_tatyana(pathname, location, cell, Nh, Nv, cy, optlevel)
% % cell : cell number
% % Nh : number of horizontal elements in RF matrix, i.e. number of delay
% %        elements
% % Nv : number of vertical elements in RF matrix, i.e. number of frequency
% %        bands
% % cy : I have no idea
% % optlevel : I have no idea
% 
% % Nh=20, Nv=30
% 
% Nwinx = 3; % number of subplot columns
% Nwiny = 4; % number of subplot rows
% between_column = 1; % for subplot command
% row = 1; % for subplot command
% start_column = 1; % for subplot command
% 
% 
% nlag_start = 1; % helps give size of RF estimate
% nlag_end = 1;   % if 1 then the RF is one image, not
% nlags = 1;      % multiple images as in vision work
% 
% cx = 1; % not used later
% 
% Nbins = 21; % used for input/output function but not RF
% Nbins_medium = 15; % used for input/output function but not RF
% Nbins_short = 14; % used for input/output function but not RF
% 
% Nparts = 4; % used to identify files, not sure what the significance
%             % is, though
% 
% clf;
% hfig = figure(cell); % one figure window for each cell


% The following code plots the spike triggered average

% close all;


location = 707;
cell = 1;
Nv = 25;
Nh = 20;
cx = 1;
cy = 1;
nlags = 1;

Nbins = 21;
Nbins_medium = 15;
Nbins_short = 14;
Nparts = 4;


% rpsta_707_1_1x25x20_1_*.dat

%--------------------------------------------------------------------
%   Data Location
%--------------------------------------------------------------------
prefix = 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\Cell_707_1_1234b\';


%--------------------------------------------------------------------
%   STA Nonlinearities
%--------------------------------------------------------------------
file_sta = sprintf('%srpsta_%u_%u_1x%ux%u_1', prefix, location, cell, Nv, Nh);

[v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_filter(file_sta, Nh, Nv, Nparts);


%--------------------------------------------------------------------
%   STA Nonlinearities
%--------------------------------------------------------------------
file_sta_fio{1} = sprintf('%srpx1pxpxt_sta_%u_%u_1x%ux%u_1_1.dat', prefix, location, cell, Nv, Nh);
file_sta_fio{2} = sprintf('%srpx1pxpxt_sta_%u_%u_1x%ux%u_1_2.dat', prefix, location, cell, Nv, Nh);
file_sta_fio{3} = sprintf('%srpx1pxpxt_sta_%u_%u_1x%ux%u_1_3.dat', prefix, location, cell, Nv, Nh);
file_sta_fio{4} = sprintf('%srpx1pxpxt_sta_%u_%u_1x%ux%u_1_4.dat', prefix, location, cell, Nv, Nh);

[fiosta] = get_dat_sta_fio(file_sta_fio, coeff_sta);



%--------------------------------------------------------------------
%   MID1
%--------------------------------------------------------------------

% rpdtest2_prelim_v1_707_1_1x25x20_1_*1.dat
optlevel = 2;

file_v1 = sprintf('%srpdtest%u_prelim_v1_%u_%u_1x%ux%u_%u', prefix, optlevel, location, cell, Nv, Nh, cy);

[v1, coeff1, projection1, mtx1] = get_auditory_filter(file_v1, Nh, Nv, Nparts);


%--------------------------------------------------------------------
%   MID2
%--------------------------------------------------------------------

% rpdtest2_prelim_v2_707_1_1x25x20_1_*1.dat
optlevel = 2;

file_v2 = sprintf('%srpdtest%u_prelim_v2_%u_%u_1x%ux%u_%u', prefix, optlevel, location, cell, Nv, Nh, cy);

[v2, coeff2, projection2, mtx2] = get_auditory_filter(file_v2, Nh, Nv, Nparts);


%--------------------------------------------------------------------
%   MID1, MID2, MID12 Nonlinearities
%--------------------------------------------------------------------
files{1} = sprintf('%srpdx1x2px_pxt_1_%u_%u_1x%ux%u_1_1.dat', prefix, location, cell, Nv, Nh);
files{2} = sprintf('%srpdx1x2px_pxt_1_%u_%u_1x%ux%u_1_2.dat', prefix, location, cell, Nv, Nh);
files{3} = sprintf('%srpdx1x2px_pxt_1_%u_%u_1x%ux%u_1_3.dat', prefix, location, cell, Nv, Nh);
files{4} = sprintf('%srpdx1x2px_pxt_1_%u_%u_1x%ux%u_1_4.dat', prefix, location, cell, Nv, Nh);

files{1} = sprintf('%srpd_prelim_x1x2px_pxt_2_%u_%u_1x%ux%u_1_1.dat', prefix, location, cell, Nv, Nh);
files{2} = sprintf('%srpd_prelim_x1x2px_pxt_2_%u_%u_1x%ux%u_1_2.dat', prefix, location, cell, Nv, Nh);
files{3} = sprintf('%srpd_prelim_x1x2px_pxt_2_%u_%u_1x%ux%u_1_3.dat', prefix, location, cell, Nv, Nh);
files{4} = sprintf('%srpd_prelim_x1x2px_pxt_2_%u_%u_1x%ux%u_1_4.dat', prefix, location, cell, Nv, Nh);

Nbins = 21;
Nparts = 4;

[fio1] = get_dat_mid1_fio(files, coeff1, coeff2, Nbins, Nparts);

[fio2] = get_dat_mid2_fio(files, coeff1, coeff2, Nbins, Nparts);

[fio12] = get_dat_mid12_fio(files, coeff1, coeff2, Nbins, Nparts);

fiosta;
fio1;
fio2;
fio12;


figure;

subplot(3,2,1);
cm = max([abs(min(min(v_sta))),max(max(v_sta))]);
vprocess = reshape(v_sta, Nv, Nh);
imagesc(vprocess)
axis tight
caxis([-cm cm])
axis tight
colorbar;
colormap('hot')
axis tight;
ylabel('STA');
title(sprintf('cell %u - %u',location, cell));

subplot(3,2,2);
pspk = fiosta.pspk_mean;
xsta = fiosta.x_mean;
pspkxsta = fiosta.ior_mean;
hold on;
plot(xsta, pspkxsta, 'ko-', 'markerfacecolor', 'k');
plot([min(xsta) max(xsta)], [pspk pspk], 'k--');



subplot(3,2,3);
cm = max([abs(min(min(v1))),max(max(v1))]);
vprocess = reshape(v1, Nv, Nh);
imagesc(vprocess)
axis tight
caxis([-cm cm])
axis tight
colorbar;
colormap('hot')
axis tight;
ylabel('MID_1');

subplot(3,2,4);
x1 = fio1.x1_mean;
pspkx1 = fio1.pspkx1_mean;
hold on;
plot(x1, pspkx1, 'ko-', 'markerfacecolor', 'k');
plot([min(x1) max(x1)], [pspk pspk], 'k--');



subplot(3,2,5);
cm = max([abs(min(min(v2))),max(max(v2))]);
vprocess = reshape(v2,Nv,Nh);
imagesc(vprocess)
axis tight
caxis([-cm cm])
axis tight
colorbar;
colormap('hot')
axis tight;
ylabel('MID_2');


subplot(3,2,6);
x2 = fio2.x2_mean;
pspkx2 = fio2.pspkx2_mean;
hold on;
plot(x2, pspkx2, 'ko-', 'markerfacecolor', 'k');
plot([min(x2) max(x2)], [pspk pspk], 'k--');

set( gcf, 'position', [168 134 916 794] );

return;








% function [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = ...
% 	plot_an_ior_improved2(fname_first, coeff, varargin)
% % summry of changes made to accomodate ripple stimuli(May 31 2005):
% % downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
% %larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
% x_mtx=[];    
% px_mtx=[];    
% ior_mtx=[];
% Nbins_short=14;
% mean_firing=0;
% 
% if isempty(varargin)
%     Nparts=8;
%     Nbins=21;
%     minpx=0;
% elseif length(varargin)==1
%     Nparts=8;
%     Nbins=varargin{1};
%     minpx=0;
% elseif length(varargin)==2
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     minpx=0;
% elseif length(varargin)==3
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
%     minpx=0;
% elseif length(varargin)==4
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
%     minpx=varargin{4}
% end
% 
% if (length(coeff)~=Nparts) 
%     length(coeff)
%     coeff
%     Nparts
%     display('wrong length of coeff vector');
%     return
% end
% 
% for i=1:Nparts
% 
%     fname = sprintf('%s_%u.dat',fname_first,i);
%     fp = fopen(fname,'r');
% 
%     if (fp==-1)
%         display('error opening file');
%         display(fname);
%         pause
%         return
%     end
% 
%     x = fread(fp,Nbins,'double');
%     px = fread(fp,Nbins,'double');
%     ind0 = find(px<minpx);
%     xm = sum(x.*px);
%     x = x-xm;
%     pxt = fread(fp,Nbins,'double');
%     pxt(ind0)=0;
% 
%     if (abs(sum(px)-1)>0.001) 
%         display(sprintf('sum of px=%f',sum(px))); 
%         return;
%     end
% 
%     rbar = fread(fp,1,'double');
%     mean_firing = mean_firing+rbar/Nparts;
% %    pause
%     Nrep_eff=fread(fp,1,'double');
%     fclose(fp);
%     x_r=x;
% %    pxt=pxt/rbar;   
%  
%     if (coeff(i)==1)
%         ior=pxt./(px+eps);
%         
%     else 
%         for j=1:Nbins
%             ior(Nbins+1-j)=pxt(j)/(px(j)+eps);
%             %x_r(Nbins+1-j)=x(j);
%         end
%         x_r=-x;
%     end
%     ior=reshape(ior,Nbins,1);
%     px=reshape(px,Nbins,1);
%     maxior=max(ior);
%     minior=min(ior);
% 
%     x_r=reshape(x_r,Nbins,1); %
%     ior=reshape(ior,Nbins,1); %
%     ior_mtx=[ior_mtx,ior];
%     px_mtx=[px_mtx,px];
%     x_mtx=[x_mtx,x_r];
%     
% end
% 
% 
% if (isempty(ior_mtx) )
%     display('empty ior_mtx in plot_an_ior');
%     return;
% end
% 
% 
% ior_mtx = reshape(ior_mtx,Nbins,Nparts);
% 
% x0 = min(min(x_mtx));
% x1 = max(max(x_mtx));
% maxmax = max([abs(x0) x1]);
% edges = linspace(-maxmax, maxmax, Nbins_short);
% 
% %edges=linspace(x0,x1,Nbins_short);
% 
% ior_rescaled = zeros(Nbins_short-1,Nparts);
% px_rescaled = zeros(Nbins_short-1,Nparts);
% npoints = zeros(Nbins_short-1,Nparts);
% 
% for i=1:Nbins_short-1
% 
%     for j=1:Nparts
% 
%         ind=find((x_mtx(:,j)>=edges(i))&(x_mtx(:,j)<edges(i+1)));
% 
%         if ~isempty(ind)
%             %      npoints(i,j)=npoints(i,j)+length(ior_mtx(ind,j));
%             npoints(i,j) = npoints(i,j)+length(ind);
%             px_rescaled(i,j) = px_rescaled(i,j)+sum(px_mtx(ind,j));
%             ior_rescaled(i,j) = ior_rescaled(i,j)+sum(ior_mtx(ind,j).*px_mtx(ind,j));
%             %            ior_rescaled(i,j)=ior_rescaled(i,j)+sum(ior_mtx(ind,j));
%         end
% 
%     end
% 
%     x_mean(i)=0.5*(edges(i)+edges(i+1));
% 
% end
% 
% %ior_rescaled=ior_rescaled./(npoints+eps);
% %px_rescaled=px_rescaled./(npoints+eps);
% %sum(px_rescaled)
% 
% ior_rescaled = ior_rescaled ./ (px_rescaled+eps);
% px_mean = mean(px_rescaled');
% 
% %sum(px_mean)
% %pause
% %px_mean=mean(px_mtx,2);
% 
% ior_mean = mean(ior_rescaled');
% 
% ior_std = sqrt(var(ior_rescaled') / Nparts);
% px_std = sqrt(var(px_rescaled') / Nparts);    
% 
% %px_std=sqrt(var(px_mtx')/Nparts);    
% %     ior_mean=mean(ior_mtx');
% %     ior_std=sqrt(var(ior_mtx')/Nparts);
% %     x_mtx=reshape(x_mtx,Nbins,Nparts);
% %     x_mean=mean(x_mtx');
% 
% subplot(Nwiny,Nwinx,Nwinx*(row-1)+column)
% 
% %px_mean=px_mean/sum(px_mean);
% 
% dx = (max(x_mean)-min(x_mean))/(Nbins_short-1);
% %dx=(max(x_mean)-min(x_mean))/(Nbins);
% 
% px_mean = px_mean ./ dx;
% px_std = px_std ./ dx;
% 
% errorbar(x_mean, ior_mean, ior_std)
% hold on
% errorbar(x_mean, px_mean, px_std, 'm')
% 
% % errorbar(x_mean, px_mean/dx, px_std/dx, 'm')
% %errorbar(linspace(x_mean(1),x_mean(Nbins_short-1),Nbins),px_mean/dx,px_std,'m')
% 
% hold off
% 
% xmax = max(x_mean)*1.1;
% xmin = min(x_mean);
% if (xmin<0) 
%     xmin=xmin*1.1;
% else
%     xmin=xmin*0.9;
% end
% xlim([xmin xmax])
% %ylim([0 max(0.05,max(ior_mean+ior_std)*1.2)])
% ylim([0 max(0.05,max(px_mean+px_std)*1.05)])
% 
% return;
% 










% function [v_mean,coeff,projection,mtx] = get_auditory_filter(fname_first, Nh, Nv, nlags, varargin)
% 
% mtx = [];
% if ( ~isempty(varargin) )
%     Nparts = varargin{1};
% else 
%     Nparts = 8; 
% end
% 
% fsize = Nv * Nh;
% Nn = fsize * nlags; % total number of elements in the STRF
% 
% for i = 1:Nparts
% 	fname = sprintf('%s_%u.dat', fname_first, i);
%    fp = fopen(fname,'rb');
% 	if ~(fp== -1)
% 		[v] = fread(fp,'double');
% %         size(v)
% %         i
% %         pause
% 		v = reshape(v,Nn,1);
% %         size(v)
% %         pause
%         %  v(Nn)=mean(v);
% 		mtx = [mtx,v];
% 		fclose(fp);
% %         figure(2)
% %         colormap('gray');
% %         v=reshape(v,fsize,nlags);
% %         for j=1:nlags
% %             subplot(Nparts,nlags,j+(i-1)*nlags);
% %             imagesc(reshape(v(:,j),Nv,Nh)')
% %             caxis([min(min(v)) max(max(v))])
% %         end
% %         colorbar;
%         %  pause
% 	else
% 		v = [];
% 	end
%     
% end
% 
% %pause
% 
% if (isempty(mtx) )
%     fname
%     error('empty mtx in plot_a_vector');
%     %return;
% end
% 
% mtx = reshape(mtx, Nn, Nparts);
% coeff(1) = 1;
% 
% 
% for i = 2:Nparts
%     %     sum((mtx(:,i)-mean(mtx(:,i))).*(mtx(:,1)-mean(mtx(:,1))))
%     %     mean(mtx(:,i))
%     %pause
%     %coeff(i)=sign(sum((mtx(:,i)-mean(mtx(:,i))).*(mtx(:,1)-mean(mtx(:,1)))));
%     coeff(i)=sign(sum(mtx(:,i).*mtx(:,1)));
%     %pause
%     if ( coeff(i) == -1 )
%         mtx(:,i)=mtx(:,i)*coeff(i);
%     end
% end
% 
% v_mean = mean(mtx')';
% v_mean = v_mean./sqrt(sum(sum(v_mean.*v_mean)));
% v_std = sqrt(sum(var(mtx'))*(Nparts-1)/Nn);
% v_mean = v_mean/v_std;
% v_mean = reshape(v_mean,fsize,nlags);
% % for testrep=1:Nparts
% %     temp_vector=mtx(:,testrep);
% %     temp_vector=temp_vector/sqrt(sum(temp_vector.*temp_vector));
% %     temp_vector=Nparts*v_mean-(Nparts-1)*temp_vector;%ATTENTION
% %     mtx(:,testrep)=temp_vector/sqrt(sum(temp_vector.*temp_vector));
% %     if (testrep >1)
% %         mtx(:,testrep)=mtx(:,testrep)*sign(sum(mtx(:,testrep).*mtx(:,1)));
% %     end
% %     
% % end
% % 
% % v_mean=mean(mtx');
% % v_mean=v_mean./sqrt(sum(sum(v_mean.*v_mean)));
% 
% % if ( nlag_end > nlags )
% %     for i=1:Nv*Nh
% %         v_mean(i,nlag_end)=0;
% %     end    
% % end
% 
% cm = max([abs(min(min(v_mean))),max(max(v_mean))]);
% 
% % nlag_start
% % nlag_end
% % size(v_mean)
% % pause
% 
% vprocess = reshape(v_mean(:,1),Nv,Nh);
% 
% 
% k = 1;
% for i = 1:Nparts
%     for j = i+1:Nparts
%         projection(k) = sum(mtx(:,i).*mtx(:,j));
%         k = k+1;
%     end
% end
% 
% return;


