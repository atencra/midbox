function mid = get_sta_mid_fio_data(filestruct)
% get_sta_mid_fio_data - put MID data files into a struct array
%    so data can be analysed
%
% mid = get_sta_mid_fio_data(filestruct)
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

% if ( isempty(findstr(current_directory, 'C:\MATLAB65\work\tatyana\Filters')) )
%    error('You need to run this function in C:\MATLAB65\work\tatyana\Filters');
% end


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
   mid(i).location = filestruct(i).location;
   mid(i).unit = filestruct(i).unit;
   mid(i).tbins = filestruct(i).tbins;
   mid(i).fbins = filestruct(i).fbins;


   % The following code plot the spike triggered average
   %====================================================================
   [v_sta, coeff_sta, projection_sta, mtx_sta] = get_auditory_strf(filestruct(i).rpsta, Nh, Nv, nlags);


   % Now plot the nonlinearity that goes along with the STA
   %====================================================================
   [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = get_sta_fior(filestruct(i).rpx1pxpxt_sta, coeff_sta, Nbins, Nparts, Nbins_short);

   rpsta.filter = v_sta;
   rpsta.coeff = coeff_sta;
   rpsta.projection = projection_sta;
   rpsta.mtx = mtx_sta;

   rpx1pxpxt_sta.mean_firing = mean_firing;
   rpx1pxpxt_sta.x = x_mean;
   rpx1pxpxt_sta.ior_mean = ior_mean;
   rpx1pxpxt_sta.ior_std = ior_std;
   rpx1pxpxt_sta.px_mean = px_mean;
   rpx1pxpxt_sta.px_std = px_std;

   mid(i).rpsta = rpsta;
   mid(i).rpx1pxpxt_sta = rpx1pxpxt_sta;



   % Process the rpdtest1 files
   %====================================================================

   % Filter for projection 1 and projection 2
   [v1, coeff1, projection1, mtx1] = get_auditory_strf(filestruct(i).rpdtest1_v1, Nh, Nv, nlags);
   [v2, coeff2, projection2, mtx2] = get_auditory_strf(filestruct(i).rpdtest1_v2, Nh, Nv, nlags);


   if ( ~isempty(filestruct(i).rpdx1x2px_pxt_1) )
      % Nonlinearity For projection v1, v2, and v1&v2 simultaneously
      [x1, ior1_mean, ior1_std, px1_mean, px1_std] = get_mid_fior_for_filter1(filestruct(i).rpdx1x2px_pxt_1, coeff1, coeff2, Nbins_medium, Nparts, Nbins_short);
      [x2, ior2_mean, ior2_std, px2_mean, px2_std] = get_mid_fior_for_filter2(filestruct(i).rpdx1x2px_pxt_1, coeff1, coeff2, Nbins_medium, Nparts, Nbins_short);
      [x12, y12, ior12_mean, ior12_std] = get_mid_fior_for_filter1_and_filter2(filestruct(i).rpdx1x2px_pxt_1, coeff1, coeff2, Nparts, Nbins_medium);
   else
      x1 = [];
      ior1_mean = [];
      ior1_std = [];
      px1_mean = [];
      px1_std = [];
      x2 = [];
      ior2_mean = [];
      ior2_std = [];
      px2_mean = [];
      px2_std = [];
      x12 = [];
      y12 = [];
      ior12_mean = [];
      ior12_std = [];
   end


   % Assign the data
   rpdtest1_v1.filter = v1;
   rpdtest1_v1.coeff = coeff1;
   rpdtest1_v1.projection = projection1;
   rpdtest1_v1.mtx = mtx1;

   rpdtest1_v2.filter = v2;
   rpdtest1_v2.coeff = coeff2;
   rpdtest1_v2.projection = projection2;
   rpdtest1_v2.mtx = mtx2;

   rpdx1x2px_pxt_1.x1 = x1;
   rpdx1x2px_pxt_1.ior1_mean = ior1_mean;
   rpdx1x2px_pxt_1.ior1_std = ior1_std;
   rpdx1x2px_pxt_1.px1_mean = px1_mean;
   rpdx1x2px_pxt_1.px1_std = px1_std;

   rpdx1x2px_pxt_1.x2 = x2;
   rpdx1x2px_pxt_1.ior2_mean = ior2_mean;
   rpdx1x2px_pxt_1.ior2_std = ior2_std;
   rpdx1x2px_pxt_1.px2_mean = px2_mean;
   rpdx1x2px_pxt_1.px2_std = px2_std;

   rpdx1x2px_pxt_1.x12 = x12;
   rpdx1x2px_pxt_1.y12 = y12;
   rpdx1x2px_pxt_1.ior12_mean = ior12_mean;
   rpdx1x2px_pxt_1.ior12_std = ior12_std;


   mid(i).rpdtest1_v1 = rpdtest1_v1;
   mid(i).rpdtest1_v2 = rpdtest1_v2;
   mid(i).rpdx1x2px_pxt_1 = rpdx1x2px_pxt_1;



   % Process the rpdtest2 files
   %====================================================================

   % Filter for projection 1 and projection 2
   [v1, coeff1, projection1, mtx1] = get_auditory_strf(filestruct(i).rpdtest2_v1, Nh, Nv, nlags);
   [v2, coeff2, projection2, mtx2] = get_auditory_strf(filestruct(i).rpdtest2_v2, Nh, Nv, nlags);


   % Nonlinearity For projection v1, v2, and v1&v2 simultaneously
   [x1_mean, ior1_mean, ior1_std, px1_mean, px1_std] = get_mid_fior_for_filter1(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, Nbins_medium, Nparts, Nbins_short);
   [x2_mean, ior2_mean, ior2_std, px2_mean, px2_std] = get_mid_fior_for_filter2(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, Nbins_medium, Nparts, Nbins_short);
   [x12, y12, ior12_mean, ior12_std] = get_mid_fior_for_filter1_and_filter2(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, Nparts, Nbins_medium);


   % Assign the data
   rpdtest2_v1.filter = v1;
   rpdtest2_v1.coeff = coeff1;
   rpdtest2_v1.projection = projection1;
   rpdtest2_v1.mtx = mtx1;

   rpdtest2_v2.filter = v2;
   rpdtest2_v2.coeff = coeff2;
   rpdtest2_v2.projection = projection2;
   rpdtest2_v2.mtx = mtx2;

   rpdx1x2px_pxt_2.x1 = x1_mean;
   rpdx1x2px_pxt_2.ior1_mean = ior1_mean;
   rpdx1x2px_pxt_2.ior1_std = ior1_std;
   rpdx1x2px_pxt_2.px1_mean = px1_mean;
   rpdx1x2px_pxt_2.px1_std = px1_std;

   rpdx1x2px_pxt_2.x2 = x2_mean;
   rpdx1x2px_pxt_2.ior2_mean = ior2_mean;
   rpdx1x2px_pxt_2.ior2_std = ior2_std;
   rpdx1x2px_pxt_2.px2_mean = px2_mean;
   rpdx1x2px_pxt_2.px2_std = px2_std;

   rpdx1x2px_pxt_2.x12 = x12;
   rpdx1x2px_pxt_2.y12 = y12;
   rpdx1x2px_pxt_2.ior12_mean = ior12_mean;
   rpdx1x2px_pxt_2.ior12_std = ior12_std;


   mid(i).rpdtest2_v1 = rpdtest2_v1;
   mid(i).rpdtest2_v2 = rpdtest2_v2;
   mid(i).rpdx1x2px_pxt_2 = rpdx1x2px_pxt_2;


   % Process the rpdbest1 files
   %====================================================================

   % Filter for projection 1 and projection 2
   [v1, coeff1, projection1, mtx1] = get_auditory_strf(filestruct(i).rpdbest1_v1, Nh, Nv, nlags);
   [v2, coeff2, projection2, mtx2] = get_auditory_strf(filestruct(i).rpdbest1_v2, Nh, Nv, nlags);

   % Assign the data
   rpdbest1_v1.filter = v1;
   rpdbest1_v1.coeff = coeff1;
   rpdbest1_v1.projection = projection1;
   rpdbest1_v1.mtx = mtx1;

   rpdbest1_v2.filter = v2;
   rpdbest1_v2.coeff = coeff2;
   rpdbest1_v2.projection = projection2;
   rpdbest1_v2.mtx = mtx2;

   mid(i).rpdbest1_v1 = rpdbest1_v1;
   mid(i).rpdbest1_v2 = rpdbest1_v2;



   % Process the rpdbest2 files
   %====================================================================

   % Filter for projection 1 and projection 2
   [v1, coeff1, projection1, mtx1] = get_auditory_strf(filestruct(i).rpdbest2_v1, Nh, Nv, nlags);
   [v2, coeff2, projection2, mtx2] = get_auditory_strf(filestruct(i).rpdbest2_v2, Nh, Nv, nlags);

   % Assign the data
   rpdbest2_v1.filter = v1;
   rpdbest2_v1.coeff = coeff1;
   rpdbest2_v1.projection = projection1;
   rpdbest2_v1.mtx = mtx1;

   rpdbest2_v2.filter = v2;
   rpdbest2_v2.coeff = coeff2;
   rpdbest2_v2.projection = projection2;
   rpdbest2_v2.mtx = mtx2;


   mid(i).rpdbest2_v1 = rpdbest2_v1;
   mid(i).rpdbest2_v2 = rpdbest2_v2;

end % (for)

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [v_mean, coeff, projection, mtx] = get_auditory_strf(files, Nh, Nv, nlags)

mtx=[];


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
% of the filters could be similar, though with similar sign. We want to correct
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




function [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = get_sta_fior(files, coeff, varargin)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

x_mtx=[];    
px_mtx=[];    
ior_mtx=[];
Nbins_short=14;
mean_firing=0;

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

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'r');

    if (fp==-1)
        display('error opening file');
        display(files{i});
        pause
        return
    end

    x = fread(fp,Nbins,'double');
    px = fread(fp,Nbins,'double');
    ind0 = find(px<minpx);
    xm = sum(x.*px);
    x = x - xm;
    pxt = fread(fp,Nbins,'double');
    pxt(ind0)=0;

    if (abs(sum(px)-1)>0.001) 
        display(sprintf('sum of px=%f',sum(px))); 
        return;
    end

    rbar = fread(fp,1,'double');
    mean_firing = mean_firing+rbar/Nparts;
    Nrep_eff = fread(fp,1,'double');
    fclose(fp);
    x_r = x;
 
    if (coeff(i)==1)
        ior=pxt./(px+eps);
        
    else 
        for j=1:Nbins
            ior(Nbins+1-j)=pxt(j)/(px(j)+eps);
        end
        x_r=-x;
    end
    ior=reshape(ior,Nbins,1);
    px=reshape(px,Nbins,1);
    maxior=max(ior);
    minior=min(ior);

    x_r=reshape(x_r,Nbins,1); %
    ior=reshape(ior,Nbins,1); %
    ior_mtx=[ior_mtx,ior];
    px_mtx=[px_mtx,px];
    x_mtx=[x_mtx,x_r];
    
end


if (isempty(ior_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end


ior_mtx = reshape(ior_mtx,Nbins,Nparts);

x0 = min(min(x_mtx));
x1 = max(max(x_mtx));
maxmax = max([abs(x0) x1]);
edges = linspace(-maxmax, maxmax, Nbins_short);

%edges=linspace(x0,x1,Nbins_short);

ior_rescaled = zeros(Nbins_short-1,Nparts);
px_rescaled = zeros(Nbins_short-1,Nparts);
npoints = zeros(Nbins_short-1,Nparts);

for i=1:Nbins_short-1

    for j=1:Nparts

        ind=find((x_mtx(:,j)>=edges(i))&(x_mtx(:,j)<edges(i+1)));

        if ~isempty(ind)
            npoints(i,j) = npoints(i,j)+length(ind);
            px_rescaled(i,j) = px_rescaled(i,j)+sum(px_mtx(ind,j));
            ior_rescaled(i,j) = ior_rescaled(i,j)+sum(ior_mtx(ind,j).*px_mtx(ind,j));
        end

    end

    x_mean(i)=0.5*(edges(i)+edges(i+1));

end

ior_rescaled = ior_rescaled ./ (px_rescaled+eps);
px_mean = mean(px_rescaled');

ior_mean = mean(ior_rescaled');

ior_std = sqrt(var(ior_rescaled') / Nparts);
px_std = sqrt(var(px_rescaled') / Nparts);    

dx = (max(x_mean)-min(x_mean))/(Nbins_short-1);

px_mean = px_mean ./ dx;
px_std = px_std ./ dx;

return;




function [x1_mean, ior1_mean, ior1_std, px1_mean, px1_std] = get_mid_fior_for_filter1(files, coeff1, coeff2, varargin)
%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

x1_mtx=[];
ior1_mtx=[];
px1_mtx=[];
x2_mtx=[];
ior2_mtx=[];
ior_mtx=[];

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
    x1=fread(fp,Nbins,'double');
    x2=fread(fp,Nbins,'double');
    px=fread(fp,Nbins*Nbins,'double');
    pxt=fread(fp,Nbins*Nbins,'double');
    rbar=fread(fp,1,'double');
    Neff=fread(fp,1,'double');

    fclose(fp);
    ind0=find(px+eps<minpx);
    pxt(ind0)=0;
    px=reshape(px,Nbins,Nbins);  %px(x1,x2)
    pxt=reshape(pxt,Nbins,Nbins);
    if (coeff1(i)==1)
        px_r=px;
        pxt_r=pxt;
        x1_r=x1;
    else
        x1_r=-x1;
    end
    if (coeff2(i)==-1)
        x2_r=-x2;
    end
    if (coeff2(i)==1) x2_r=x2; end
    ior_mtx=[ior_mtx,reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
    pxt1=sum(pxt_r');
    px1=sum(px_r');
    pxt2=sum(pxt_r);
    px2=sum(px_r);
    xm1=sum(x1_r.*px1');
    xm2=sum(x2_r.*px2');
    x1_r=x1_r-xm1;
    x2_r=x2_r-xm2;
    ior2=pxt2./(px2+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1=pxt1./(px1+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1=reshape(ior1,Nbins,1); % stolbets
    px1=reshape(px1,Nbins,1); % stolbets
    ior2=reshape(ior2,Nbins,1); % stolbets

    x1_r=reshape(x1_r,Nbins,1); % stolbets
    x2_r=reshape(x2_r,Nbins,1); % stolbets
    ior1_mtx=[ior1_mtx,ior1];
    ior2_mtx=[ior2_mtx,ior2];
    px1_mtx=[px1_mtx,px1];
    x1_mtx=[x1_mtx,x1_r];
    x2_mtx=[x2_mtx,x2_r];
end
if (isempty(ior_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end

x0 = min(min(x1_mtx));
x1 = max(max(x1_mtx));
maxmax = max([abs(x0) x1]);
edges = linspace(-maxmax, maxmax, Nbins_short);

ior1_rescaled = zeros(Nbins_short-1,Nparts);
px1_rescaled = zeros(Nbins_short-1,Nparts);
npoints = zeros(Nbins_short-1,Nparts);

for i=1:Nbins_short-1
    for j=1:Nparts
        ind = find((x1_mtx(:,j)>=edges(i))&(x1_mtx(:,j)<edges(i+1)));
        if ~isempty(ind)
            npoints(i,j)=1;%npoints(i,j)+length(ior1_mtx(ind,j));
            ior1_rescaled(i,j)=ior1_rescaled(i,j)+sum(ior1_mtx(ind,j).*px1_mtx(ind,j));
            px1_rescaled(i,j)=px1_rescaled(i,j)+sum(px1_mtx(ind,j));
        end
    end
    x1_mean(i) = 0.5*(edges(i)+edges(i+1));%mean(x1_mtx');
end

nsamples = sum(npoints');
ior1_rescaled = ior1_rescaled./(px1_rescaled+eps);
px1_mean = sum(px1_rescaled')./(nsamples+eps);
px1_std = sqrt((sum(px1_rescaled'.^2)./(nsamples+eps)-px1_mean.^2)./(nsamples-1));
sumpx1=sum(px1_mean);

% if abs(sumpx1-1)>0.1
%     sumpx1
%     fname
%     px1_mean
%     px1_rescaled
%     px1_mtx
%     pause
% end

px1_mean = px1_mean./sum(px1_mean);
dx = (max(x1_mean)-min(x1_mean))/(Nbins_short-1);

ior1_mean = sum(ior1_rescaled')./(nsamples+eps);
ior1_std = sqrt((sum(ior1_rescaled'.^2)./(nsamples+eps)-ior1_mean.^2)./(nsamples-1).*(nsamples+eps));

px1_mean = px1_mean ./ dx;
px1_std = px1_std ./ dx;

return;





function [x2_mean, ior2_mean, ior2_std, px2_mean, px2_std] = get_mid_fior_for_filter2(files, coeff1, coeff2, varargin)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density

x1_mtx=[];    
ior1_mtx=[];
x2_mtx=[];    
ior2_mtx=[];  
px2_mtx=[];  
ior_mtx=[];
minpx=0;

%Nparts=8;
%Nbins=21;

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
elseif length(varargin)>=4
    Nbins=varargin{1};
    Nparts=varargin{2};
    Nbins_short=varargin{3}
    minpx=varargin{4}
end

if (length(coeff1)~=Nparts) 
    length(coeff1)
    Nparts
    Nbins
    Nbins_short
    pause
    display('wrong length of coeff vector');
    return
end

for i=1:Nparts

    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'r');

    if (fp==-1)
        display('error opening file');
        display(files{i});
        return
    end
    x1=fread(fp,Nbins,'double');
    x2=fread(fp,Nbins,'double');
    px=fread(fp,Nbins*Nbins,'double');
    pxt=fread(fp,Nbins*Nbins,'double');
    px=reshape(px,Nbins,Nbins);  %px(x1,x2)  
    pxt=reshape(pxt,Nbins,Nbins);
    rbar=fread(fp,1,'double');
    Nrep_eff=fread(fp,1,'double');
    fclose(fp);
    ind0=find(px<minpx);
    pxt(ind0)=0;
    if (coeff1(i)==1)
        px_r=px;
        pxt_r=pxt;
        x1_r=x1;
    else 
        x1_r=-x1;
    end
    if (coeff2(i)==-1)
        x2_r=-x2;
    end
    
    if (coeff2(i)==1) x2_r=x2; end
    ior_mtx=[ior_mtx,reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
    pxt1=sum(pxt_r');
    px1=sum(px_r');
    pxt2=sum(pxt_r);
    px2=sum(px_r);
    ior2=pxt2./(px2+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);   
    ior1=pxt1./(px1+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1=reshape(ior1,Nbins,1); % stolbets
    ior2=reshape(ior2,Nbins,1); % stolbets
    xm1=sum(x1_r.*px1');
    xm2=sum(x2_r.*px2');
    x1_r=x1_r-xm1;
    x2_r=x2_r-xm2;
    x1_r=reshape(x1_r,Nbins,1); % stolbets
    x2_r=reshape(x2_r,Nbins,1); % stolbets
    px2=reshape(px2,Nbins,1); % stolbets
    ior1_mtx=[ior1_mtx,ior1];
    px2_mtx=[px2_mtx,px2];
    ior2_mtx=[ior2_mtx,ior2];
    x1_mtx=[x1_mtx,x1_r];
    x2_mtx=[x2_mtx,x2_r];
end
if (isempty(ior_mtx) )
    display('empty ior_mtx in plot_an_ior');
    return;
end

% x0=min(min(x2_mtx));
% x1=max(max(x2_mtx));
% edges=linspace(x0,x1,Nbins_short);

x0 = min(min(x2_mtx));
x1 = max(max(x2_mtx));
maxmax = max([abs(x0) x1]);
edges = linspace(-maxmax, maxmax, Nbins_short);

ior2_rescaled=zeros(Nbins_short-1,Nparts);
px2_rescaled=zeros(Nbins_short-1,Nparts);
npoints=zeros(Nbins_short-1,Nparts);

for i=1:Nbins_short-1
    for j=1:Nparts
        ind=find((x2_mtx(:,j)>=edges(i))&(x2_mtx(:,j)<edges(i+1)));
        if ~isempty(ind)
            npoints(i,j)=1;%npoints(i,j)+length(ior2_mtx(ind,j));
            ior2_rescaled(i,j)=ior2_rescaled(i,j)+sum(ior2_mtx(ind,j).*px2_mtx(ind,j));
            px2_rescaled(i,j)=px2_rescaled(i,j)+sum(px2_mtx(ind,j));
            %        ior2_rescaled(i,j)=ior2_rescaled(i,j)+sum(ior2_mtx(ind,j));
        end
    end
    x2_mean(i)=0.5*(edges(i)+edges(i+1));
end

ior2_rescaled=ior2_rescaled./(px2_rescaled+eps);

nsamples = sum(npoints');

px2_mean = sum(px2_rescaled')./(nsamples+eps);
px2_std = sqrt(var(px2_rescaled')*(Nparts-1)./(nsamples-1)./(nsamples+eps));

sumpx2 = sum(px2_mean);

% if abs(sumpx2-1)>0.1 
%     sumpx2
%     fname
%     px2_mean
%     px2_rescaled
%     px2_mtx
%     pause
% end

px2_mean = px2_mean./sum(px2_mean);
dx = (max(x2_mean)-min(x2_mean))/(Nbins_short-1);

ior2_mean=sum(ior2_rescaled')./(nsamples+eps);
ior2_std=sqrt((sum(ior2_rescaled'.^2)./(nsamples+eps)-ior2_mean.^2)./(nsamples-1));

px2_mean = px2_mean ./ dx;
px2_std = px2_std ./ dx;

% subplot(Nwiny,Nwinx,(row-1)*Nwinx+column)
% errorbar(x2_mean, ior2_mean, ior2_std,'g')
% hold on
% errorbar(x2_mean, px2_mean, px2_std,'k')
% xmax = max(x2_mean)*1.1;
% xmin = min(x2_mean)*1.1;
% xlim([xmin xmax])
% ylim([0 max([0.05,ior2_mean+ior2_std])*1.2])

return;





function [x12, y12, ior12_mean, ior12_std] = get_mid_fior_for_filter1_and_filter2(files, coeff1, coeff2, Nparts, Nbins, varargin)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density1_mtx=[];
% also properly rotated the two-dimensional matrix, and its labels

ior1_mtx=[];
x2_mtx=[];
x1_mtx=[];
ior2_mtx=[];
ior_mtx=[];

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
    %fname = sprintf('%s_%u.dat',fname_first,i);
    fp = fopen(files{i},'r');
    if (fp==-1)
        display('error opening file');
        display(files{i});
        return
    end

    x1=fread(fp,Nbins,'double');
    x2=fread(fp,Nbins,'double');
    px=fread(fp,Nbins*Nbins,'double');
    pxt=fread(fp,Nbins*Nbins,'double');
    px=reshape(px,Nbins,Nbins);
    pxt=reshape(pxt,Nbins,Nbins);

    if (coeff1(i)==1)
        px_r = px;
        pxt_r = pxt;
        x1_r = x1;
    else
        for j=1:Nbins
            px_r(Nbins+1-j,:) = px(j,:);
            pxt_r(Nbins+1-j,:) = pxt(j,:);
            x1_r = -x1;
        end
    end
    if (coeff2(i)==-1)
        for j=1:Nbins
            px_r(:,Nbins+1-j) = px(:,j);
            pxt_r(:,Nbins+1-j) = pxt(:,j);
            x2_r = -x2;
        end
    end
    if (coeff2(i)==1) x2_r = x2; end
    ior_mtx=[ior_mtx,reshape(pxt_r./(px_r+eps),Nbins*Nbins,1)];
    px_mtx=[ior_mtx,reshape(px_r,Nbins*Nbins,1)];
    pxt1=sum(pxt_r,2);
    px1=sum(px_r,2);
    pxt2=sum(pxt_r,1);
    px2=sum(px_r,1);
    ior2=pxt2./(px2+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1=pxt1./(px1+eps);%*Nspikes(ind)/repetitions/Ntrials(ind);
    ior1=reshape(ior1,Nbins,1); % stolbets
    ior2=reshape(ior2,Nbins,1); % stolbets
    x1_r=reshape(x1_r,Nbins,1); % stolbets
    x2_r=reshape(x2_r,Nbins,1); % stolbets
    ior1_mtx=[ior1_mtx,ior1];
    ior2_mtx=[ior2_mtx,ior2];
    x1_mtx=[x1_mtx,x1_r];
    x2_mtx=[x2_mtx,x2_r];
    if (i==1)
        rbar=fread(fp,1,'double');
        Nrep_eff=fread(fp,1,'double');
    end
    fclose(fp);
end

if (isempty(ior_mtx) )
   display('empty ior_mtx in plot_an_ior');
   return;
end

x0 = min(min(x1_mtx));
x1 = max(max(x1_mtx));
maxmax = max([abs(x0) x1]);
xedges = linspace(-maxmax, maxmax, Nbins);

y0 = min(min(x2_mtx));
y1 = max(max(x2_mtx));
maxmax = max([abs(y0) y1]);
yedges = linspace(-maxmax, maxmax, Nbins);

ior_rescaled=zeros(Nbins*Nbins,Nparts);
for i=1:Nbins-1
    for k=1:Nbins-1
        for j=1:Nparts
            xind = find((x1_mtx(:,j)>=xedges(i))&(x1_mtx(:,j)<xedges(i+1)));
            yind = find((x2_mtx(:,j)>=yedges(k))&(x2_mtx(:,j)<yedges(k+1)));
            temp = 0;
            sum_px = 0;
            for ii=1:length(xind)
                for jj=1:length(yind)
                    temp = temp+ior_mtx((yind(jj)-1)*Nbins+xind(ii),j).*px_mtx((yind(jj)-1)*Nbins+xind(ii),j);
                    sum_px = sum_px+sum(px_mtx((yind(jj)-1)*Nbins+xind(ii),j));
                end
            end
            if (temp>0)
                temp = temp/sum_px;
            end
            ior_rescaled((k-1)*Nbins+i,j) = ior_rescaled((k-1)*Nbins+i,j)+temp;
        end
    end
    x1_mean(i) = 0.5*(xedges(i)+xedges(i+1));
    x2_mean(i) = 0.5*(yedges(i)+yedges(i+1));
end


ior_mean = mean(ior_rescaled');
ior12_mean = reshape(ior_mean, Nbins, Nbins)';

ior_std = sqrt(var(ior_rescaled')/Nparts);
ior12_std = reshape(ior_std, Nbins, Nbins)';

x12 = x1_mean;
y12 = x2_mean;

% subplot(Nwiny,Nwinx,(row-1)*Nwinx+column);
% imagesc(x12, y12, ior12_mean);
% xlabel('MID1 Proj');
% ylabel('MID2 Proj');
% axis square;
% colorbar;
% ind=find(ior_mean==max(ior_mean));
% indy=mod(ind,Nbins);
% indx=(ind-indy)/Nbins+1;
% 
% hold on
% 
% if ~isempty(varargin)
% 
% disp('~isempty(varargin)');
% 
%     xstart=mean(x1_mean)
%     ystart=0
%     aspect_ratio=(max(x2_mean)-min(x2_mean))/(max(x1_mean)-min(x1_mean));
%     x_extend=max(abs([x1_mean,x2_mean]))/aspect_ratio;
%     xend=x_extend*(alpha1+alpha2*non_ortho)*aspect_ratio;
%     yend=x_extend*(beta1+beta2*non_ortho)*aspect_ratio;
%     h=line([xstart,xstart+xend],[ystart,yend+ystart]);
%     set(h,'Linewidth',5,'Color','w');
%     h=text(1.1*(xend+xstart),1.1*(yend+ystart),'1');
%     set(h,'Color','w');
%     
%     xend=x_extend*(alpha2+alpha1*non_ortho)*aspect_ratio;
%     yend=x_extend*(beta2+beta1*non_ortho)*aspect_ratio;
%     %    h=line([xstart,ystart],[xend+xstart,yend+ystart]);
%     %    h=line([ystart,xstart],[yend+ystart,xend+xstart]);
%     h=line([xstart,xstart+xend],[ystart,yend+ystart]);
%     set(h,'Linewidth',3,'Color','w');
%     h=text(1.1*(xend+xstart),(yend+ystart)*1.1,'2');
%     set(h,'Color','w');
%     % pause
% end

return;























