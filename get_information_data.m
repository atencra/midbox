function info = get_information_data(filestruct)
% get_information_data - get information values for MID analysis
%
% info = get_information_data(file)
%
% The folder C:\MATLAB65\work\tatyana\Filters contains all the data from
% Tatyana's maximally informative dimensions analysis. This function goes
% through that folder, extracts all the files for one location, and then
% creates a struct array of the file information, where each element of the
% struct array corresponds to one single unit in the location.
%
% location : a scalar. May be 500, 515, 516, 517, 519, 532, 534, 537, 602,
% 603, 604, 608, 609, 611, 614, or 616.
%
%
% caa 4/11/06

if ( nargin~=1 )
   error('You need one input arg.');
end

current_directory = pwd;

if ( isempty(findstr(current_directory, 'C:\MATLAB65\work\tatyana\Filters')) )
   error('You need to run this function in C:\MATLAB65\work\tatyana\Filters');
end


% Parameters that Tatyana uses in her programs
nlag_start = 1; % helps give size of RF estimate
nlag_end = 1;   % if 1 then the RF is one image, not multiple images
nlags = 1;      % multiple images as in vision work
Nbins = 21; % used for input/output function but not RF
Nbins_medium = 15; % used for input/output function but not RF
Nbins_short = 14; % used for input/output function but not RF
Nparts = 4; % used to identify files, not sure what the significance
            % is, though

% Number of time and frequency bins in the receptive fields
Nh = filestruct(1).tbins;
Nv = filestruct(1).fbins;


for i = 1:length(filestruct)
i
   info(i).location = filestruct(i).location;
   info(i).unit = filestruct(i).unit;
   info(i).tbins = filestruct(i).tbins;
   info(i).fbins = filestruct(i).fbins;


   % The following code plot the spike triggered average
   %====================================================================
   [v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_strf(filestruct(i).rpsta, Nh, Nv, nlags);


   % Now plot the nonlinearity that goes along with the STA
   %====================================================================
   [x_mtx, px_mtx, pxspk_mtx, pspk_mtx, pspkx_mtx, info_sta] = ...
      get_sta_information(filestruct(i).rpx1pxpxt_sta, coeff_sta, Nbins, Nparts, Nbins_short);

   sta.x = x_mtx;
   sta.px = px_mtx;
   sta.pxspk = pxspk_mtx;
   sta.pspk = pspk_mtx;
   sta.pspkx = pspkx_mtx;
   sta.information = info_sta;

   info(i).sta = sta;



   % Process the rpdtest2 files
   %====================================================================

   % Filter for projection 1 and projection 2
   [v1, coeff1, projection1, mtx1] = get_auditory_strf(filestruct(i).rpdtest2_v1, Nh, Nv, nlags);
   [v2, coeff2, projection2, mtx2] = get_auditory_strf(filestruct(i).rpdtest2_v2, Nh, Nv, nlags);

   [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, info_mid1] = ...
      get_mid1_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, Nbins_medium, Nparts, Nbins_short);

   mid1.x = x1_mtx;
   mid1.px = px1_mtx;
   mid1.pxspk = px1spk_mtx;
   mid1.pspk = pspk_mtx;
   mid1.pspkx = pspkx1_mtx;
   mid1.information = info_mid1;

   info(i).mid1 = mid1;


   [x2_mtx, px2_mtx, px2spk_mtx, pspk_mtx, pspkx2_mtx, info_mid2] = ...
      get_mid2_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, Nbins_medium, Nparts, Nbins_short);

   mid2.x = x2_mtx;
   mid2.px = px2_mtx;
   mid2.pxspk = px2spk_mtx;
   mid2.pspk = pspk_mtx;
   mid2.pspkx = pspkx2_mtx;
   mid2.information = info_mid2;

   info(i).mid2 = mid2;



   [x1_mtx, x2_mtx, px1x2_mtx, px1x2spk_mtx, pspk_mtx, pspkx1x2_mtx, info_mid12] = ...
      get_mid1_mid2_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, Nparts, Nbins_medium);

   mid12.x1 = x1_mtx;
   mid12.x2 = x2_mtx;
   mid12.px1x2 = px1x2_mtx;
   mid12.px1x2spk = px1x2spk_mtx;
   mid12.pspk = pspk_mtx;
   mid12.pspkx1x2 = pspkx1x2_mtx;
   mid12.information = info_mid12;

   info(i).mid12 = mid12;

end % (for i)

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v_mean, coeff, projection, mtx] = get_auditory_strf(files, Nh, Nv, nlags)

mtx = [];

fsize = Nv*Nh;
Nn = fsize*nlags; % total number of elements in the STRF
Nparts = length(files);

for i=1:Nparts

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'rb');

    if ~(fp== -1)
        [v]=fread(fp,'double');
        v=reshape(v,Nn,1);
        mtx=[mtx,v];
        fclose(fp);
    else v=[];
    end
    
end

if (isempty(mtx) )
    files{1}
    display('empty mtx in plot_a_vector');
    return;
end

mtx = reshape(mtx,Nn,Nparts);
coeff(1) = 1;


% This just corrects a possible pos/neg symmetry issue. It's possible some
% of the filters could be similar, though with opposite sign. We want to correct
% for that and make sure all the filters overlap correctly
for i=2:Nparts
    coeff(i)=sign(sum(mtx(:,i).*mtx(:,1)));
    if (coeff(i) == -1)
        mtx(:,i)=mtx(:,i)*coeff(i);
    end
end

v_mean = mean(mtx')';
v_mean = v_mean ./ sqrt(sum(sum(v_mean.*v_mean)));
v_std = sqrt(sum(var(mtx'))*(Nparts-1)/Nn);
v_mean = v_mean ./ v_std;
v_mean = reshape(v_mean,fsize,nlags);

k=1;

for i=1:Nparts
    for j=i+1:Nparts
        projection(k) = sum(mtx(:,i).*mtx(:,j)); % find every possible correlation between the 
        k=k+1;                                   % 4 estimated STAs?
    end
end

return;




function [x_mtx, px_mtx, pxspk_mtx, pspk_mtx, pspkx_mtx, information] = get_sta_information(files, coeff, varargin)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density


x_mtx = [];    
px_mtx = [];
pxspk_mtx = [];
pspkx_mtx = [];
pspk_mtx = [];

Nbins_short = 14;

if isempty(varargin)
    Nparts=8;
    Nbins=21;
    minpx=0;
elseif length(varargin)==1
    Nparts=8;
    Nbins=varargin{1};
    minpx=0;
elseif length(varargin)==2
    Nbins=varargin{1};
    Nparts=varargin{2};
    minpx=0;
elseif length(varargin)==3
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3};
    minpx=0;
elseif length(varargin)==4
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3};
    minpx=varargin{4}
end

if (length(coeff)~=Nparts) 
    length(coeff)
    coeff
    Nparts
    display('wrong length of coeff vector');
    return
end

for i=1:Nparts

   fp = fopen(files{i},'r');

   if (fp==-1)
       display('error opening file');
       display(files{i});
       pause
       return
   end

   x = fread(fp,Nbins,'double');
   px = fread(fp,Nbins,'double');
   pxspk_pspk = fread(fp,Nbins,'double');
   pspk = fread(fp,1,'double');
   Nrep_eff = fread(fp,1,'double');
   fclose(fp);

   ind0 = find(px < minpx);
   pxspk_pspk(ind0)=0;
   pxspk = pxspk_pspk ./ pspk; % prob of projection given a spike

    xm = sum(x.*px);
    x = x - xm;

    if (abs(sum(px)-1)>0.001) 
        display(sprintf('sum of px=%f',sum(px))); 
        return;
    end

% sum(px)
% sum(pxspk)
% plot(x,px,'r-', x, pxspk, 'k-');
% pause

    x_r = x;
 
   if (coeff(i)==1)
      pspkx = pxspk_pspk ./ (px+eps);
   else 
      for j=1:Nbins
         pspkx(Nbins+1-j) = pxspk_pspk(j) ./ (px(j)+eps);
      end
      x_r = -x;
   end

   x_r = reshape(x_r,Nbins,1);
   px = reshape(px,Nbins,1);
   pxspk = reshape(pxspk,Nbins,1);
   pspkx = reshape(pspkx,Nbins,1);

   x_mtx = [x_mtx x_r];
   px_mtx = [px_mtx px];
   pxspk_mtx = [pxspk_mtx pxspk];
   pspk_mtx = [pspk_mtx pspk];
   pspkx_mtx = [pspkx_mtx pspkx];

end

if ( size(pspkx_mtx,2) ~= Nparts )
    error('empty pxspk_pspk');
end

% size(x_mtx)
% size(px_mtx)
% size(pxspk_mtx)
% size(pspkx_mtx)
% size(pspk_mtx)
% pause


% Calculate information values using a different technique:
pxs1 = pxspk_mtx(:,1);
pxs2 = pxspk_mtx(:,2);
pxs3 = pxspk_mtx(:,3);
pxs4 = pxspk_mtx(:,4);

px1 = px_mtx(:,1);
px2 = px_mtx(:,2);
px3 = px_mtx(:,3);
px4 = px_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];


% pspkx_mtx = reshape(pspkx_mtx,Nbins,Nparts);
% 
% xmin = min(min(x_mtx));
% xmax = max(max(x_mtx));
% maxmax = max([abs(xmin) xmax]);
% edges = linspace(-maxmax, maxmax, Nbins_short);
% 
% %edges=linspace(x0,x1,Nbins_short);
% 
% pspkx_rescaled = zeros(Nbins_short-1,Nparts);
% px_rescaled = zeros(Nbins_short-1,Nparts);
% pxspk_rescaled = zeros(Nbins_short-1,Nparts);
% npoints = zeros(Nbins_short-1,Nparts);
% 
% for i=1:Nbins_short-1
% 
%     for j=1:Nparts
% 
%         ind = find((x_mtx(:,j)>=edges(i))&(x_mtx(:,j)<edges(i+1)));
% 
%         if ~isempty(ind)
%             npoints(i,j) = npoints(i,j)+length(ind);
%             px_rescaled(i,j) = px_rescaled(i,j) + sum(px_mtx(ind,j));
%             pxspk_rescaled(i,j) = pxspk_rescaled(i,j) + sum(pxspk_mtx(ind,j));
%             pspkx_rescaled(i,j) = pspkx_rescaled(i,j) + sum(pspkx_mtx(ind,j).*px_mtx(ind,j));
%         end
% 
%     end
% 
%     x_mean(i)=0.5*(edges(i)+edges(i+1));
% 
% end
% 
% pspkx_rescaled = pspkx_rescaled ./ (px_rescaled+eps);
% 
% pspkx_mean = mean(pspkx_rescaled');
% pspkx_std = sqrt(var(pspkx_rescaled') / Nparts);
% 
% px_mean = mean(px_rescaled');
% px_std = sqrt(var(px_rescaled') / Nparts);    
% 
% dx = (max(x_mean)-min(x_mean))/(Nbins_short-1);
% 
% px_mean = px_mean ./ dx;
% px_std = px_std ./ dx;
% 
% 
% % The following code will also calculate the nonlinearity by 
% % using a loess procedure to find the appropriate shape
% % to the nonlinearity
% projection = -10:0.25:10;
% xmin = min(x_mtx(:))
% xmax = max(x_mtx(:))
% 
% ipos = find(projection < xmax);
% ineg = find(projection > xmin);
% 
% pmax = projection(max(ipos))
% pmin = projection(min(ineg))
% pause
% 
% 
% newx = pmin:0.25:pmax;
% alpha = .1; % minimal smoothing
% 
% lambda = 1; % linear local curve fitting
% newylin = loess(x_mtx(:), pspkx_mtx(:), newx, alpha, lambda);
% newylin(newylin<0) = 0;
% 
% 
% lambda = 2; % quadratic local curve fitting
% newyquad = loess(x_mtx(:), pspkx_mtx(:), newx, alpha, lambda);
% newyquad(newyquad<0) = 0;
% 
% % the linear local curve fitting approach is the best
% 
% figure;
% 
% subplot(2,1,1);
% hold on;
% plot(x_mtx(:,1), pspkx_mtx(:,1), 'ko');
% plot(x_mtx(:,2), pspkx_mtx(:,2), 'ko');
% plot(x_mtx(:,3), pspkx_mtx(:,3), 'ko');
% plot(x_mtx(:,4), pspkx_mtx(:,4), 'ko');
% plot(newx, newylin, 'r-', 'linewidth', 2);
% plot(newx, newyquad, 'g-', 'linewidth', 2);
% xlim([-6 6])
% 
% subplot(2,1,2);
% hold on;
% plot(newx, newylin, 'ro-', 'linewidth', 2);
% plot(newx, newyquad, 'go-', 'linewidth', 2);
% plot(x_mean, pspkx_mean, 'ko-');
% xlim([-6 6])
% 
% pause




return;




function [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, information] = ...
   get_mid1_information(files, coeff1, coeff2, varargin)

%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

% Variable name definitions and what they mean:
% 
% pspx1 = P(spike|x1)
% psp = P(spike)
% px1sp = P(x1|spike)
% px1 = P(x1)
% px1x2 = P(x1,x2)
% px1x2sp = P(x1,x2|spike)
% px1sp_psp = P(x1|sp) * P(spike)
%
% px1spk_mean 
% px1_mean
% pspkx1_mean
% pspk_mean

x1_mtx = [];
pspkx1_mtx = [];
px1_mtx = [];
px1spk_mtx = [];
px1spk_pspk_mtx = [];

pspkx1x2_mtx=[];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];



%Nparts=8;
%Nbins=21;

minpx=0;

if isempty(varargin)
    Nparts=8;
    Nbins=21;
    Nbins_short=19;
elseif length(varargin)==1
    Nparts=8;
    Nbins_short=19;
    Nbins=varargin{1};
elseif length(varargin)==2
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=19;
elseif length(varargin)==3
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3};
elseif length(varargin)==4
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3};
    minpx=varargin{4};
end

Nbins_short;
Nbins;
Nparts;


if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   display('wrong length of coeff vector');
   return
end

Nok=0;

for i=1:Nparts

   %fname = sprintf('%s_%u.dat',fname_first,i);
   fp = fopen(files{i},'r');

   if (fp==-1)
      display('error opening file');
      display(files{i});
      return
   end

   x1 = fread(fp, Nbins, 'double'); % first filter projection values
   x2 = fread(fp, Nbins, 'double'); % second filter projection values
   px1x2 = fread(fp, Nbins*Nbins, 'double'); % prior distr
   px1x2spk_pspk = fread(fp, Nbins*Nbins, 'double'); % really rbar * pxt
   pspk = fread(fp, 1, 'double');
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike
   Neff = fread(fp, 1, 'double');

   fclose(fp);

   % Get rid of machine error values
   ind0 = find(px1x2+eps<minpx);
   px1x2spk(ind0) = 0;

   % Reshape into a matrix so that marginals can be calculated
   px1x2 = reshape(px1x2, Nbins, Nbins);  %p(x1,x2)
   px1x2spk = reshape(px1x2spk, Nbins, Nbins);
   px1x2spk_pspk = reshape(px1x2spk_pspk, Nbins, Nbins);

   if (coeff1(i)==1)
       px1x2_r = px1x2;
       px1x2spk_r = px1x2spk;
       x1_r = x1;
   else
       x1_r = -x1;
   end

   if (coeff2(i) == -1)
       x2_r = -x2;
   end

   if (coeff2(i)==1) x2_r=x2; end

   % assign the complete probability distributions
   px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk,Nbins*Nbins,1)]; % P(x1,x2|spike)
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
   pspk_mtx = [pspk_mtx pspk]; % P(spike)
   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)

   % assign the x1 marginal probability distributions

   % P(x1|spike) * P(spike)
   px1spk_pspk = sum(px1x2spk_pspk,2); % sum across columns to get the first marginal
   px1spk_pspk = px1spk_pspk(:);
   px1spk_pspk_mtx = [px1spk_pspk_mtx px1spk_pspk];

% sum(px1sp_psp)

   % P(x1|spike)
   px1spk = sum(px1x2spk,2);
   px1spk = px1spk(:);
   px1spk_mtx = [px1spk_mtx px1spk]; 

   % P(x1)
   px1 = sum(px1x2,2);
   px1 = px1(:);
   px1_mtx = [px1_mtx px1]; 

   % P(spike|x1)
   pspkx1 = px1spk_pspk ./ (px1+eps);
   pspkx1 = pspkx1(:);
   pspkx1_mtx = [pspkx1_mtx pspkx1]; 

   xm1 = sum(x1_r .* px1); % expected value: x * p(x) of the projection
   x1_r = x1_r - xm1;
   x1_r = x1_r(:);
   x1_mtx = [x1_mtx x1_r];

end

if (isempty(pspkx1x2_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end

% sum(pspkx1_mtx)
% sum(px1spk_mtx)
% sum(px1_mtx)
% pause

% Calculate information values using a different technique:
pxs1 = px1spk_mtx(:,1);
pxs2 = px1spk_mtx(:,2);
pxs3 = px1spk_mtx(:,3);
pxs4 = px1spk_mtx(:,4);

px1 = px1_mtx(:,1);
px2 = px1_mtx(:,2);
px3 = px1_mtx(:,3);
px4 = px1_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

return;



function [x2_mtx, px2_mtx, px2spk_mtx, pspk_mtx, pspkx2_mtx, information] = ...
   get_mid2_information(files, coeff1, coeff2, varargin)
%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

% Variable name definitions and what they mean:
% 
% pspx1 = P(spike|x1)
% psp = P(spike)
% px1sp = P(x1|spike)
% px1 = P(x1)
% px1x2 = P(x1,x2)
% px1x2sp = P(x1,x2|spike)
% px1sp_psp = P(x1|sp) * P(spike)
%
% px1spk_mean 
% px1_mean
% pspkx1_mean
% pspk_mean

x2_mtx = [];
pspkx2_mtx = [];
px2_mtx = [];
px2spk_mtx = [];
px2spk_pspk_mtx = [];

pspkx1x2_mtx=[];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];

minpx=0;

if isempty(varargin)
    Nparts=8;
    Nbins=21;
    Nbins_short=19;
elseif length(varargin)==1
    Nparts=8;
    Nbins_short=19;
    Nbins=varargin{1};
elseif length(varargin)==2
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=19;
elseif length(varargin)==3
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3};
elseif length(varargin)==4
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3};
    minpx=varargin{4};
end

Nbins_short;
Nbins;
Nparts;


if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   error('wrong length of coeff vector');
end

Nok=0;

for i=1:Nparts

   fp = fopen(files{i},'r');

   if (fp==-1)
      error(sprintf('error opening file %s', files{i}));
   end

   x1 = fread(fp, Nbins, 'double'); % first filter projection values
   x2 = fread(fp, Nbins, 'double'); % second filter projection values
   px1x2 = fread(fp, Nbins*Nbins, 'double'); % prior distr
   px1x2spk_pspk = fread(fp, Nbins*Nbins, 'double'); % really rbar * pxt
   pspk = fread(fp, 1, 'double');
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike
   Neff = fread(fp, 1, 'double');

   fclose(fp);

   % Get rid of machine error values
   ind0 = find(px1x2+eps<minpx);
   px1x2spk(ind0) = 0;

   % Reshape into a matrix so that marginals can be calculated
   px1x2 = reshape(px1x2, Nbins, Nbins);  %p(x1,x2)
   px1x2spk = reshape(px1x2spk, Nbins, Nbins);
   px1x2spk_pspk = reshape(px1x2spk_pspk, Nbins, Nbins);

   if ( coeff1(i) == -1 )
       x1_r = -x1;
   else
       x1_r = x1;
   end

   if ( coeff2(i) == -1 )
       x2_r = -x2;
   else
       x2_r = x2;
   end

   if (coeff2(i)==1) x2_r=x2; end

   % assign the complete probability distributions
   px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk,Nbins*Nbins,1)]; % P(x1,x2|spike)
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
   pspk_mtx = [pspk_mtx pspk]; % P(spike)
   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)

   % assign the x2 marginal probability distributions

   % P(x2|spike) * P(spike)
   px2spk_pspk = sum(px1x2spk_pspk,1); % sum across columns to get the first marginal
   px2spk_pspk = px2spk_pspk(:);
   px2spk_pspk_mtx = [px2spk_pspk_mtx px2spk_pspk];


   % P(x2|spike)
   px2spk = sum(px1x2spk,1);
   px2spk = px2spk(:);
   px2spk_mtx = [px2spk_mtx px2spk];

   % P(x2)
   px2 = sum(px1x2,1);
   px2 = px2(:);
   px2_mtx = [px2_mtx px2];

   % P(spike|x2)
   pspkx2 = px2spk_pspk ./ (px2+eps);
   pspkx2 = pspkx2(:);
   pspkx2_mtx = [pspkx2_mtx pspkx2];

   xm2 = sum( x2_r .* px2 ); % expected value: x * p(x) of the projection
   x2_r = x2_r - xm2;
   x2_r = x2_r(:);
   x2_mtx = [x2_mtx x2_r];

end


% Calculate information values:
pxs1 = px2spk_mtx(:,1);
pxs2 = px2spk_mtx(:,2);
pxs3 = px2spk_mtx(:,3);
pxs4 = px2spk_mtx(:,4);

px1 = px2_mtx(:,1);
px2 = px2_mtx(:,2);
px3 = px2_mtx(:,3);
px4 = px2_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

 return;




function [x1_mtx, x2_mtx, px1x2_mtx, px1x2spk_mtx, pspk_mtx, pspkx1x2_mtx, information] = ...
      get_mid1_mid2_information(files, coeff1, coeff2, Nparts, Nbins, varargin)
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within 
% the larger bin, and dividing the two numbers, also introduced dx, and probability distribution 
% is now shown as a probability density1_mtx=[];
% also properly rotated the two-dimensional matrix, and its labels

x1_mtx = [];
x2_mtx = [];

pspkx1x2_mtx = [];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];

if isempty(varargin)
    alpha1=0;
    alpha2=0;
    beta1=0;
    beat2=0;
else
    alpha1=varargin{1};
    alpha2=varargin{2};
    beta1=varargin{3};
    beta2=varargin{4};
    non_ortho=varargin{5}
end

if (length(coeff1)~=Nparts)
    length(coeff1)
    display('wrong length of coeff vector');
    return
end

for i=1:Nparts

    fp = fopen(files{i},'r');

    if (fp==-1)
        error(sprintf('error opening file %s', files{i}));
        return
    end

   x1 = fread(fp,Nbins,'double');
   x2 = fread(fp,Nbins,'double');
   px1x2 = fread(fp,Nbins*Nbins,'double');
   px1x2spk_pspk = fread(fp,Nbins*Nbins,'double');
   pspk = fread(fp, 1, 'double');

   px1x2 = reshape(px1x2,Nbins,Nbins); % make it a matrix
   px1x2spk_pspk = reshape(px1x2spk_pspk,Nbins,Nbins); % make it a matrix
   px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike

   fclose(fp);

    if (coeff1(i)==1)
        px1x2_r = px1x2;
        px1x2spk_r = px1x2spk;
        px1x2spk_pspk_r = px1x2spk_pspk;
        x1_r = x1;
    else
        for j=1:Nbins
            px1x2_r(Nbins+1-j,:) = px1x2(j,:);
            px1x2spk_pspk_r(Nbins+1-j,:) = px1x2spk_pspk(j,:);
            px1x2spk_r(Nbins+1-j,:) = px1x2spk(j,:);
            x1_r = -x1;
        end
    end

   if (coeff2(i)==-1)
      for j=1:Nbins
         px1x2_r(:,Nbins+1-j) = px1x2(:,j);
         px1x2spk_pspk_r(:,Nbins+1-j) = px1x2spk_pspk(:,j);
         px1x2spk_r(:,Nbins+1-j) = px1x2spk(:,j);
         x2_r = -x2;
      end
   end

   if (coeff2(i)==1) x2_r = x2; end

   pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk_r./(px1x2_r+eps),Nbins*Nbins,1)];

   px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk_r,Nbins*Nbins,1)]; % make it a column vector
   px1x2_mtx = [px1x2_mtx reshape(px1x2_r,Nbins*Nbins,1)]; % make it a column vector
   pspk_mtx = [pspk_mtx pspk];
   px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk_r,Nbins*Nbins,1)]; % make it a column vector

   x1_r = reshape(x1_r,Nbins,1); % stolbets
   x2_r = reshape(x2_r,Nbins,1); % stolbets

   x1_mtx = [x1_mtx x1_r];
   x2_mtx = [x2_mtx x2_r];

end


% Calculate information values:
pxs1 = px1x2spk_mtx(:,1);
pxs2 = px1x2spk_mtx(:,2);
pxs3 = px1x2spk_mtx(:,3);
pxs4 = px1x2spk_mtx(:,4);

px1 = px1x2_mtx(:,1);
px2 = px1x2_mtx(:,2);
px3 = px1x2_mtx(:,3);
px4 = px1x2_mtx(:,4);

ind1 = find(pxs1>0 & px1>0);
ind2 = find(pxs2>0 & px2>0);
ind3 = find(pxs3>0 & px3>0);
ind4 = find(pxs4>0 & px4>0);

info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );

information = [info1 info2 info3 info4];

return;






% 
% function [x1_mean, px1_mean, px1_std] = ...
%    get_information_for_mid_filter1(files, coeff1, coeff2, varargin)
% %plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% %% summry of changes made to accomodate ripple stimuli(May 31 2005):
% %% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
% %larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
% 
% % Variable name definitions and what they mean:
% % 
% % pspx1 = P(spike|x1)
% % psp = P(spike)
% % px1sp = P(x1|spike)
% % px1 = P(x1)
% % px1x2 = P(x1,x2)
% % px1x2sp = P(x1,x2|spike)
% % px1sp_psp = P(x1|sp) * P(spike)
% %
% % px1spk_mean 
% % px1_mean
% % pspkx1_mean
% % pspk_mean
% 
% x1_mtx = [];
% ior1_mtx = [];
% pspx1_mtx = [];
% px1_mtx = [];
% px1sp_mtx = [];
% px1sp_psp_mtx = [];
% 
% ior_mtx=[];
% pspx1x2_mtx=[];
% psp_mtx = [];
% px1x2_mtx = [];
% px1x2sp_mtx = [];
% px1x2sp_psp_mtx = [];
% 
% 
% 
% %Nparts=8;
% %Nbins=21;
% 
% minpx=0;
% 
% if isempty(varargin)
%     Nparts=8;
%     Nbins=21;
%     Nbins_short=19;
% elseif length(varargin)==1
%     Nparts=8;
%     Nbins_short=19;
%     Nbins=varargin{1};
% elseif length(varargin)==2
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=19;
% elseif length(varargin)==3
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
% elseif length(varargin)==4
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
%     minpx=varargin{4};
% end
% 
% Nbins_short;
% Nbins;
% Nparts;
% 
% 
% if ( length(coeff1) ~= Nparts )
%    length(coeff1)
%    Nparts
%    display('wrong length of coeff vector');
%    return
% end
% 
% Nok=0;
% 
% for i=1:Nparts
% 
%    %fname = sprintf('%s_%u.dat',fname_first,i);
%    fp = fopen(files{i},'r');
% 
%    if (fp==-1)
%       display('error opening file');
%       display(files{i});
%       return
%    end
% 
%    x1 = fread(fp, Nbins, 'double'); % first filter projection values
%    x2 = fread(fp, Nbins, 'double'); % second filter projection values
%    px1x2 = fread(fp, Nbins*Nbins, 'double'); % prior distr
%    px1x2sp_psp = fread(fp, Nbins*Nbins, 'double'); % really rbar * pxt
%    psp = fread(fp, 1, 'double');
%    px1x2sp = px1x2sp_psp ./ psp; % prob of projection given a spike
%    Neff = fread(fp, 1, 'double');
% 
%    fclose(fp);
% 
%    % Get rid of machine error values
%    ind0 = find(px1x2+eps<minpx);
%    px1x2sp(ind0) = 0;
% 
%    % Reshape into a matrix so that marginals can be calculated
%    px1x2 = reshape(px1x2, Nbins, Nbins);  %p(x1,x2)
%    px1x2sp = reshape(px1x2sp, Nbins, Nbins);
%    px1x2sp_psp = reshape(px1x2sp_psp, Nbins, Nbins);
% 
%    if (coeff1(i)==1)
%        px1x2_r = px1x2;
%        px1x2sp_r = px1x2sp;
%        x1_r = x1;
%    else
%        x1_r = -x1;
%    end
% 
%    if (coeff2(i) == -1)
%        x2_r = -x2;
%    end
% 
%    if (coeff2(i)==1) x2_r=x2; end
% 
%    % assign the complete probability distributions
%    px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
%    px1x2sp_mtx = [px1x2sp_mtx reshape(px1x2sp,Nbins*Nbins,1)]; % P(x1,x2|spike)
%    px1x2sp_psp_mtx = [px1x2sp_psp_mtx reshape(px1x2sp_psp,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
%    psp_mtx = [psp_mtx psp]; % P(spike)
% %    ior_mtx = [ior_mtx reshape(px1x2sp_psp ./ (px1x2+eps), Nbins*Nbins, 1)];
%    pspx1x2_mtx = [pspx1x2_mtx reshape(px1x2sp_psp ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)
% 
%    % assign the x1 marginal probability distributions
% 
%    % P(x1|spike) * P(spike)
%    px1sp_psp = sum(px1x2sp_psp,2); % sum across columns to get the first marginal
%    px1sp_psp = px1sp_psp(:);
%    px1sp_psp_mtx = [px1sp_psp_mtx px1sp_psp]; 
% 
%    % P(x1|spike)
%    px1sp = sum(px1x2sp,2);
%    px1sp = px1sp(:);
%    px1sp_mtx = [px1sp_mtx px1sp]; 
% 
%    % P(x1)
%    px1 = sum(px1x2,2);
%    px1 = px1(:);
%    px1_mtx = [px1_mtx px1]; 
% 
% %    ior1 = px1sp_psp ./ (px1+eps);
% %    ior1 = ior1(:);
% %    ior1_mtx = [ior1_mtx ior1];
% 
%    % P(spike|x1)
%    pspx1 = px1sp_psp ./ (px1+eps);
%    pspx1 = pspx1(:);
%    pspx1_mtx = [pspx1_mtx pspx1]; 
% 
%    xm1 = sum(x1_r .* px1); % expected value: x * p(x) of the projection
%    x1_r = x1_r - xm1;
%    x1_r = x1_r(:);
%    x1_mtx = [x1_mtx x1_r];
% 
% end
% 
% if (isempty(pspx1x2_mtx) )
%     display('empty ior_mtx in plot_an_ior');
%     return;
% end
% 
% % if (isempty(ior_mtx) )
% %     display('empty ior_mtx in plot_an_ior');
% %     return;
% % end
% 
% xmin = min(min(x1_mtx));
% xmax = max(max(x1_mtx));
% maxmax = max([abs(xmin) xmax]);
% edges = linspace(-maxmax, maxmax, Nbins_short);
% 
% % Possible other bins and edges that could be used
% % bincenters = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
% % bincenters = bincenters(:);
% % binedges = [-100 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5 100];
% % binedges = binedges(:);
% 
% % ior1_rescaled = zeros(Nbins_short-1,Nparts);
% pspx1_rescaled = zeros(Nbins_short-1,Nparts);
% px1_rescaled = zeros(Nbins_short-1,Nparts);
% px1sp_rescaled = zeros(Nbins_short-1,Nparts);
% npoints = zeros(Nbins_short-1,Nparts);
% 
% for i=1:Nbins_short-1
%    for j=1:Nparts
%       ind = find((x1_mtx(:,j)>=edges(i))&(x1_mtx(:,j)<edges(i+1)));
%       if ~isempty(ind)
%          npoints(i,j)=1; %npoints(i,j)+length(ior1_mtx(ind,j));
% %          ior1_rescaled(i,j)= ior1_rescaled(i,j) + sum(ior1_mtx(ind,j) .* px1_mtx(ind,j));
%          pspx1_rescaled(i,j)= pspx1_rescaled(i,j) + sum(pspx1_mtx(ind,j) .* px1_mtx(ind,j));
%          px1_rescaled(i,j)= px1_rescaled(i,j) + sum(px1_mtx(ind,j));
%          px1sp_rescaled(i,j) = px1sp_rescaled(i,j) + sum(px1sp_mtx(ind,j));
%       end
%    end
%    x1_mean(i) = 0.5*(edges(i)+edges(i+1));%mean(x1_mtx');
% end
% 
% nsamples = sum(npoints');
% 
% % P(spike|x1)
% pspx1_rescaled = pspx1_rescaled./(px1_rescaled+eps);
% pspx1_mean = sum(pspx1_rescaled')./(nsamples+eps);
% pspx1_std = sqrt((sum(pspx1_rescaled'.^2)./(nsamples+eps)-pspx1_mean.^2)./(nsamples-1).*(nsamples+eps));
% 
% % ior1_rescaled = ior1_rescaled./(px1_rescaled+eps);
% % ior1_mean = sum(ior1_rescaled')./(nsamples+eps);
% % ior1_std = sqrt((sum(ior1_rescaled'.^2)./(nsamples+eps)-ior1_mean.^2)./(nsamples-1).*(nsamples+eps));
% 
% % P(x1)
% px1_mean = sum(px1_rescaled')./(nsamples+eps);
% px1_std = sqrt((sum(px1_rescaled'.^2)./(nsamples+eps)-px1_mean.^2)./(nsamples-1));
% sumpx1 = sum(px1_mean);
% px1_mean = px1_mean./sum(px1_mean);
% dx = (max(x1_mean)-min(x1_mean))/(Nbins_short-1);
% 
% %px1_mean = px1_mean ./ dx;
% %px1_std = px1_std ./ dx;
% 
% % P(x1|sp)
% px1sp_mean = sum(px1sp_rescaled')./(nsamples+eps);
% px1sp_std = sqrt((sum(px1sp_rescaled'.^2)./(nsamples+eps)-px1_mean.^2)./(nsamples-1));
% sumpx1sp = sum(px1sp_mean);
% px1sp_mean = px1sp_mean./sum(px1sp_mean);
% 
% psp_mean = mean(psp_mtx);
% 
% % pspx1_mean = ior1_mean;
% fio = pspx1_mean / psp_mean; % F = P(sp|x1) / P(sp) = P(x1|xp) / P(x1)
% 
% index = find(fio>0);
% information = sum( px1_mean(index) .* fio(index) .* log2(fio(index)) )
% 
% index = find(px1sp_mean>0 & px1_mean>0);
% information = sum( px1sp_mean(index) .* log2( px1sp_mean(index) ./ px1_mean(index) ) )
% 
% disp('pause')
% pause
% 
% % Calculate information values using a different technique:
% pxs1 = px1sp_mtx(:,1);
% pxs2 = px1sp_mtx(:,2);
% pxs3 = px1sp_mtx(:,3);
% pxs4 = px1sp_mtx(:,4);
% 
% px1 = px1_mtx(:,1);
% px2 = px1_mtx(:,2);
% px3 = px1_mtx(:,3);
% px4 = px1_mtx(:,4);
% 
% ind1 = find(pxs1>0 & px1>0);
% ind2 = find(pxs2>0 & px2>0);
% ind3 = find(pxs3>0 & px3>0);
% ind4 = find(pxs4>0 & px4>0);
% 
% info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
% info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
% info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
% info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );
% 
% info = [info1 info2 info3 info4]
% mean(info)
% 
% disp('mid1 info')
% pause
% 
% 
% return;
% 
% 
% function [x2_mean, px2_mean, px2_std] = ...
%    get_information_for_mid_filter2(files, coeff1, coeff2, varargin)
% %plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% %% summry of changes made to accomodate ripple stimuli(May 31 2005):
% %% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
% %larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
% 
% % Variable name definitions and what they mean:
% % 
% % pspx1 = P(spike|x1)
% % psp = P(spike)
% % px1sp = P(x1|spike)
% % px1 = P(x1)
% % px1x2 = P(x1,x2)
% % px1x2sp = P(x1,x2|spike)
% % px1sp_psp = P(x1|sp) * P(spike)
% %
% % px1spk_mean 
% % px1_mean
% % pspkx1_mean
% % pspk_mean
% 
% x2_mtx = [];
% pspx2_mtx = [];
% px2_mtx = [];
% px2sp_mtx = [];
% px2sp_psp_mtx = [];
% 
% pspx1x2_mtx=[];
% psp_mtx = [];
% px1x2_mtx = [];
% px1x2sp_mtx = [];
% px1x2sp_psp_mtx = [];
% 
% minpx=0;
% 
% if isempty(varargin)
%     Nparts=8;
%     Nbins=21;
%     Nbins_short=19;
% elseif length(varargin)==1
%     Nparts=8;
%     Nbins_short=19;
%     Nbins=varargin{1};
% elseif length(varargin)==2
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=19;
% elseif length(varargin)==3
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
% elseif length(varargin)==4
%     Nbins=varargin{1};
%     Nparts=varargin{2};
%     Nbins_short=varargin{3};
%     minpx=varargin{4};
% end
% 
% Nbins_short;
% Nbins;
% Nparts;
% 
% 
% if ( length(coeff1) ~= Nparts )
%    length(coeff1)
%    Nparts
%    error('wrong length of coeff vector');
% end
% 
% Nok=0;
% 
% for i=1:Nparts
% 
%    fp = fopen(files{i},'r');
% 
%    if (fp==-1)
%       error(sprintf('error opening file %s', files{i}));
%    end
% 
%    x1 = fread(fp, Nbins, 'double'); % first filter projection values
%    x2 = fread(fp, Nbins, 'double'); % second filter projection values
%    px1x2 = fread(fp, Nbins*Nbins, 'double'); % prior distr
%    px1x2sp_psp = fread(fp, Nbins*Nbins, 'double'); % really rbar * pxt
%    psp = fread(fp, 1, 'double');
%    px1x2sp = px1x2sp_psp ./ psp; % prob of projection given a spike
%    Neff = fread(fp, 1, 'double');
% 
%    fclose(fp);
% 
%    % Get rid of machine error values
%    ind0 = find(px1x2+eps<minpx);
%    px1x2sp(ind0) = 0;
% 
%    % Reshape into a matrix so that marginals can be calculated
%    px1x2 = reshape(px1x2, Nbins, Nbins);  %p(x1,x2)
%    px1x2sp = reshape(px1x2sp, Nbins, Nbins);
%    px1x2sp_psp = reshape(px1x2sp_psp, Nbins, Nbins);
% 
%    if ( coeff1(i) == -1 )
%        x1_r = -x1;
%    else
%        x1_r = x1;
%    end
% 
%    if ( coeff2(i) == -1 )
%        x2_r = -x2;
%    else
%        x2_r = x2;
%    end
% 
%    if (coeff2(i)==1) x2_r=x2; end
% 
%    % assign the complete probability distributions
%    px1x2_mtx = [px1x2_mtx reshape(px1x2,Nbins*Nbins,1)]; % P(x1,x2)
%    px1x2sp_mtx = [px1x2sp_mtx reshape(px1x2sp,Nbins*Nbins,1)]; % P(x1,x2|spike)
%    px1x2sp_psp_mtx = [px1x2sp_psp_mtx reshape(px1x2sp_psp,Nbins*Nbins,1)]; % P(x1,x2|spike) * P(spike)
%    psp_mtx = [psp_mtx psp]; % P(spike)
%    pspx1x2_mtx = [pspx1x2_mtx reshape(px1x2sp_psp ./ (px1x2+eps), Nbins*Nbins, 1)]; % P(spike|x1,x2)
% 
%    % assign the x2 marginal probability distributions
% 
%    % P(x2|spike) * P(spike)
%    px2sp_psp = sum(px1x2sp_psp,1); % sum across columns to get the first marginal
%    px2sp_psp = px2sp_psp(:);
%    px2sp_psp_mtx = [px2sp_psp_mtx px2sp_psp];
% 
% 
%    % P(x2|spike)
%    px2sp = sum(px1x2sp,1);
%    px2sp = px2sp(:);
%    px2sp_mtx = [px2sp_mtx px2sp];
% 
%    % P(x2)
%    px2 = sum(px1x2,1);
%    px2 = px2(:);
%    px2_mtx = [px2_mtx px2];
% 
%    % P(spike|x2)
%    pspx2 = px2sp_psp ./ (px2+eps);
%    pspx2 = pspx2(:);
%    pspx2_mtx = [pspx2_mtx pspx2];
% 
%    xm2 = sum( x2_r .* px2 ); % expected value: x * p(x) of the projection
%    x2_r = x2_r - xm2;
%    x2_r = x2_r(:);
%    x2_mtx = [x2_mtx x2_r];
% 
% end
% 
% if (isempty(pspx1x2_mtx) )
%     error('empty pspx1x2_mtx in nonlinearity file.');
% end
% 
% xmin = min(min(x2_mtx));
% xmax = max(max(x2_mtx));
% maxmax = max([abs(xmin) xmax]);
% edges = linspace(-maxmax, maxmax, Nbins_short);
% 
% % Possible other bins and edges that could be used
% % bincenters = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
% % bincenters = bincenters(:);
% % binedges = [-100 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5 100];
% % binedges = binedges(:);
% 
% pspx2_rescaled = zeros(Nbins_short-1,Nparts);
% px2_rescaled = zeros(Nbins_short-1,Nparts);
% px2sp_rescaled = zeros(Nbins_short-1,Nparts);
% npoints = zeros(Nbins_short-1,Nparts);
% 
% for i=1:Nbins_short-1
%    for j=1:Nparts
%       ind = find((x2_mtx(:,j)>=edges(i))&(x2_mtx(:,j)<edges(i+1)));
%       if ~isempty(ind)
%          npoints(i,j)=1; 
%          pspx2_rescaled(i,j) = pspx2_rescaled(i,j) + sum( pspx2_mtx(ind,j) .* px2_mtx(ind,j) );
%          px2_rescaled(i,j) = px2_rescaled(i,j) + sum( px2_mtx(ind,j) );
%          px2sp_rescaled(i,j) = px2sp_rescaled(i,j) + sum( px2sp_mtx(ind,j) );
%       end
%    end
%    x2_mean(i) = 0.5*(edges(i)+edges(i+1));
% end
% 
% nsamples = sum(npoints');
% 
% % P(spike|x2)
% pspx2_rescaled = pspx2_rescaled ./ (px2_rescaled+eps);
% pspx2_mean = sum(pspx2_rescaled') ./ (nsamples+eps);
% pspx2_std = sqrt((sum(pspx2_rescaled'.^2) ./ (nsamples+eps)-pspx2_mean.^2) ./ (nsamples-1) .* (nsamples+eps));
% 
% % P(x2)
% px2_mean = sum(px2_rescaled')./(nsamples+eps);
% px2_std = sqrt((sum(px2_rescaled'.^2)./(nsamples+eps)-px2_mean.^2)./(nsamples-1));
% sumpx2 = sum(px2_mean);
% px2_mean = px2_mean./sum(px2_mean);
% 
% 
% % P(x2|sp)
% px2sp_mean = sum(px2sp_rescaled') ./ (nsamples+eps);
% px2sp_std = sqrt((sum(px2sp_rescaled'.^2) ./ (nsamples+eps)-px2_mean.^2) ./ (nsamples-1));
% sumpx2sp = sum(px2sp_mean);
% px2sp_mean = px2sp_mean ./ sum(px2sp_mean);
% 
% psp_mean = mean(psp_mtx);
% 
% % pspx1_mean = ior1_mean;
% fio = pspx2_mean / psp_mean; % F = P(sp|x1) / P(sp) = P(x1|xp) / P(x1)
% 
% index = find(fio>0);
% information = sum( px2_mean(index) .* fio(index) .* log2(fio(index)) )
% 
% index = find(px2sp_mean>0 & px2_mean>0);
% information = sum( px2sp_mean(index) .* log2( px2sp_mean(index) ./ px2_mean(index) ) )
% 
% disp('pause')
% pause
% 
% 
% 
% % Calculate information values using a different technique:
% pxs1 = px2sp_mtx(:,1);
% pxs2 = px2sp_mtx(:,2);
% pxs3 = px2sp_mtx(:,3);
% pxs4 = px2sp_mtx(:,4);
% 
% px1 = px2_mtx(:,1);
% px2 = px2_mtx(:,2);
% px3 = px2_mtx(:,3);
% px4 = px2_mtx(:,4);
% 
% ind1 = find(pxs1>0 & px1>0);
% ind2 = find(pxs2>0 & px2>0);
% ind3 = find(pxs3>0 & px3>0);
% ind4 = find(pxs4>0 & px4>0);
% 
% info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
% info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
% info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
% info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );
% 
% info = [info1 info2 info3 info4]
% mean(info)
% 
% disp('mid2 info')
% pause
% 
% return;
% 
% 
% 
% 
% function [x12, y12] = get_information_for_mid_filter1_and_filter2(files, coeff1, coeff2, Nparts, Nbins, varargin)
% %% summry of changes made to accomodate ripple stimuli(May 31 2005):
% %% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
% %larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density1_mtx=[];
% % also properly rotated the two-dimensional matrix, and its labels
% 
% ior1_mtx=[];
% x2_mtx=[];
% x1_mtx=[];
% ior2_mtx=[];
% ior_mtx=[];
% 
% pspkx1x2_mtx=[];
% pspk_mtx = [];
% px1x2_mtx = [];
% px1x2spk_mtx = [];
% px1x2spk_pspk_mtx = [];
% 
% if isempty(varargin)
%     alpha1=0;
%     alpha2=0;
%     beta1=0;
%     beat2=0;
% else
%     alpha1=varargin{1};
%     alpha2=varargin{2};
%     beta1=varargin{3};
%     beta2=varargin{4};
%     non_ortho=varargin{5}
% end
% 
% if (length(coeff1)~=Nparts)
%     length(coeff1)
%     display('wrong length of coeff vector');
%     return
% end
% 
% for i=1:Nparts
% 
%     fp = fopen(files{i},'r');
% 
%     if (fp==-1)
%         error(sprintf('error opening file %s', files{i}));
%         return
%     end
% 
%    x1 = fread(fp,Nbins,'double');
%    x2 = fread(fp,Nbins,'double');
%    px1x2 = fread(fp,Nbins*Nbins,'double');
%    px1x2spk_pspk = fread(fp,Nbins*Nbins,'double');
%    pspk = fread(fp, 1, 'double');
% 
%    px1x2 = reshape(px1x2,Nbins,Nbins);
%    px1x2spk_pspk = reshape(px1x2spk_pspk,Nbins,Nbins);
%    px1x2spk = px1x2spk_pspk ./ pspk; % prob of projection given a spike
% 
%    fclose(fp);
% 
%     if (coeff1(i)==1)
%         px1x2_r = px1x2;
%         px1x2spk_r = px1x2spk;
%         px1x2spk_pspk_r = px1x2spk_pspk;
%         x1_r = x1;
%     else
%         for j=1:Nbins
%             px1x2_r(Nbins+1-j,:) = px1x2(j,:);
%             px1x2spk_pspk_r(Nbins+1-j,:) = px1x2spk_pspk(j,:);
%             px1x2spk_r(Nbins+1-j,:) = px1x2spk(j,:);
%             x1_r = -x1;
%         end
%     end
% 
%    if (coeff2(i)==-1)
%       for j=1:Nbins
%          px1x2_r(:,Nbins+1-j) = px1x2(:,j);
%          px1x2spk_pspk_r(:,Nbins+1-j) = px1x2spk_pspk(:,j);
%          px1x2spk_r(:,Nbins+1-j) = px1x2spk(:,j);
%          x2_r = -x2;
%       end
%    end
% 
%    if (coeff2(i)==1) x2_r = x2; end
% 
%    pspkx1x2_mtx = [pspkx1x2_mtx reshape(px1x2spk_pspk_r./(px1x2_r+eps),Nbins*Nbins,1)];
% %    ior_mtx = [ior_mtx reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
% 
%    px1x2spk_mtx = [px1x2spk_mtx reshape(px1x2spk_r,Nbins*Nbins,1)];
%    px1x2_mtx = [px1x2_mtx reshape(px1x2_r,Nbins*Nbins,1)];
%    pspk_mtx = [pspk_mtx pspk];
%    px1x2spk_pspk_mtx = [px1x2spk_pspk_mtx reshape(px1x2spk_pspk_r,Nbins*Nbins,1)];
% 
%    x1_r = reshape(x1_r,Nbins,1); % stolbets
%    x2_r = reshape(x2_r,Nbins,1); % stolbets
% 
%    x1_mtx = [x1_mtx x1_r];
%    x2_mtx = [x2_mtx x2_r];
% 
% end
% 
% % x1_mtx
% % x2_mtx
% % px1x2_mtx
% % return;
% 
% if (isempty(px1x2spk_pspk_mtx) )
%    error('empty px1x2spk_pspk_mtx');
% end
% 
% xmin = min(min(x1_mtx));
% xmax = max(max(x1_mtx));
% maxmax = max([abs(xmin) xmax]);
% xedges = linspace(-maxmax, maxmax, Nbins);
% 
% ymin = min(min(x2_mtx));
% ymax = max(max(x2_mtx));
% maxmax = max([abs(ymin) ymax]);
% yedges = linspace(-maxmax, maxmax, Nbins);
% 
% pspkx1x2_rescaled = zeros(Nbins*Nbins,Nparts);
% px1x2_rescaled = zeros(Nbins*Nbins,Nparts);
% px1x2spk_rescaled = zeros(Nbins*Nbins,Nparts);
% 
% for i=1:Nbins-1
%     for k=1:Nbins-1
%         for j=1:Nparts
%             xind = find((x1_mtx(:,j)>=xedges(i))&(x1_mtx(:,j)<xedges(i+1)));
%             yind = find((x2_mtx(:,j)>=yedges(k))&(x2_mtx(:,j)<yedges(k+1)));
%             temp = 0;
%             sum_px1x2 = 0;
%             sum_px1x2spk = 0;
%             for ii=1:length(xind)
%                 for jj=1:length(yind)
%                     temp = temp + pspkx1x2_mtx((yind(jj)-1)*Nbins+xind(ii),j) .* px1x2_mtx((yind(jj)-1)*Nbins+xind(ii),j);
%                     sum_px1x2 = sum_px1x2 + sum(px1x2_mtx((yind(jj)-1)*Nbins+xind(ii),j));
%                     sum_px1x2spk = sum_px1x2spk + sum(px1x2spk_mtx((yind(jj)-1)*Nbins+xind(ii),j));
%                 end
%             end
%             if (temp>0)
%                 temp = temp/sum_px1x2;
%             end
%             pspkx1x2_rescaled((k-1)*Nbins+i,j) = pspkx1x2_rescaled((k-1)*Nbins+i,j) + temp;
%             px1x2_rescaled((k-1)*Nbins+i,j) = px1x2_rescaled((k-1)*Nbins+i,j) + sum_px1x2;
%             px1x2spk_rescaled((k-1)*Nbins+i,j) = px1x2spk_rescaled((k-1)*Nbins+i,j) + sum_px1x2spk;
%         end
%     end
%     x1_mean(i) = 0.5*(xedges(i)+xedges(i+1));
%     x2_mean(i) = 0.5*(yedges(i)+yedges(i+1));
% end
% 
% 
% % P(spike|x1,x2)
% pspkx1x2_mean = mean(pspkx1x2_rescaled');
% pspkx1x2_mean = reshape(pspkx1x2_mean, Nbins, Nbins)';
% pspkx1x2_mean = pspkx1x2_mean ./ sum(sum(pspkx1x2_mean));
% 
% pspkx1x2_std = sqrt(var(pspkx1x2_rescaled')/Nparts);
% pspkx1x2_std = reshape(pspkx1x2_std, Nbins, Nbins)';
% 
% 
% % P(x1,x2)
% px1x2_mean = mean(px1x2_rescaled');
% px1x2_mean = reshape(px1x2_mean, Nbins, Nbins)';
% px1x2_mean = px1x2_mean ./ sum(sum(px1x2_mean));
% 
% px1x2_std = sqrt(var(px1x2_rescaled')/Nparts);
% px1x2_std = reshape(px1x2_std, Nbins, Nbins)';
% 
% 
% % P(x1,x2|spike)
% px1x2spk_mean = mean(px1x2spk_rescaled');
% px1x2spk_mean = reshape(px1x2spk_mean, Nbins, Nbins)';
% px1x2spk_mean = px1x2spk_mean ./ sum(sum(px1x2spk_mean));
% 
% px1x2spk_std = sqrt(var(px1x2spk_rescaled')/Nparts);
% px1x2spk_std = reshape(px1x2spk_std, Nbins, Nbins)';
% 
% x12 = x1_mean;
% y12 = x2_mean;
% 
% pspk_mean = mean(pspk_mtx);
% 
% % F(x1,x2) = P(spike|x1,x2) / P(spike) = P(x1,x2|spike) / P(x1,x2)
% fx1x2 = pspkx1x2_mean ./ pspk_mean; 
% fx1x2 = px1x2spk_mean ./ px1x2_mean; 
% 
% 
% index = find(fx1x2>0);
% information = sum( sum( px1x2_mean(index) .* fx1x2(index) .* log2(fx1x2(index)) ) )
% 
% % Calculate information values using a different technique:
% pxs1 = px1x2spk_mtx(:,1);
% pxs2 = px1x2spk_mtx(:,2);
% pxs3 = px1x2spk_mtx(:,3);
% pxs4 = px1x2spk_mtx(:,4);
% 
% px1 = px1x2_mtx(:,1);
% px2 = px1x2_mtx(:,2);
% px3 = px1x2_mtx(:,3);
% px4 = px1x2_mtx(:,4);
% 
% ind1 = find(pxs1>0 & px1>0);
% ind2 = find(pxs2>0 & px2>0);
% ind3 = find(pxs3>0 & px3>0);
% ind4 = find(pxs4>0 & px4>0);
% 
% info1 = sum( pxs1(ind1) .* log2( pxs1(ind1) ./ px1(ind1) ) );
% info2 = sum( pxs2(ind2) .* log2( pxs2(ind2) ./ px2(ind2) ) );
% info3 = sum( pxs3(ind3) .* log2( pxs3(ind3) ./ px3(ind3) ) );
% info4 = sum( pxs4(ind4) .* log2( pxs4(ind4) ./ px4(ind4) ) );
% 
% info = [info1 info2 info3 info4]
% mean(info)
% 
% disp('mid1 and mid2 info')
% pause
% 
% 
% % index = find(px2sp_mean>0 & px2_mean>0);
% % information = sum( px2sp_mean(index) .* log2( px2sp_mean(index) ./ px2_mean(index) ) )
% % 
% % subplot(Nwiny,Nwinx,(row-1)*Nwinx+column);
% % imagesc(x12, y12, ior12_mean);
% % xlabel('MID1 Proj');
% % ylabel('MID2 Proj');
% % axis square;
% % colorbar;
% % ind=find(ior_mean==max(ior_mean));
% % indy=mod(ind,Nbins);
% % indx=(ind-indy)/Nbins+1;
% 
% return;




















