function info = get_rippledir_detailed_filters_nonlinearities(filestruct)
% get_information_data - get information values for MID analysis
%
% info = get_rippledir_detailed_filters_nonlinearities(file)
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

% current_directory = pwd;
% 
% if ( isempty(findstr(current_directory, 'C:\MATLAB65\work\tatyana\Filters')) )
%    error('You need to run this function in C:\MATLAB65\work\tatyana\Filters');
% end


% Parameters that Tatyana uses in her programs
nlag_start = 1; % helps give size of RF estimate
nlag_end = 1;   % if 1 then the RF is one image, not multiple images
% nlags = 1;      % multiple images as in vision work
Nbins = 21; % used for input/output function but not RF
Nbins_medium = 15; % used for input/output function but not RF
Nbins_short = 14; % used for input/output function but not RF
Nparts = 4; % used to identify files, not sure what the significance
            % is, though

% Number of time and frequency bins in the receptive fields
Nh = filestruct(1).nh;
Nv = filestruct(1).nv;
nlags = filestruct(1).nlags;


for i = 1:length(filestruct)

   info(i).location = filestruct(i).location;
   info(i).unit = filestruct(i).unit;

   info(i).x0 = filestruct(i).x0;
   info(i).nh = filestruct(i).nh;
   info(i).nv = filestruct(i).nv;
   info(i).nlags = filestruct(i).nlags;

   info(i).tbins = filestruct(i).tbins;
   info(i).fbins = filestruct(i).fbins;

   if ( ~isempty(filestruct(i).rpsta) )

      % get the spike triggered average
      %====================================================================
      [v_sta, coeff_sta, projection_sta, mtx_sta] = ...
         get_auditory_strf(filestruct(i).rpsta, Nh, Nv, nlags);


      % get the sta nonlinearity data and information 
      %====================================================================
      [x_mtx, px_mtx, pxspk_mtx, pspk_mtx, pspkx_mtx, info_sta] = ...
         get_sta_information(filestruct(i).rpx1pxpxt_sta, coeff_sta, Nbins, Nparts);

      sta.x = x_mtx;
      sta.px = px_mtx;
      sta.pxspk = pxspk_mtx;
      sta.pspk = pspk_mtx;
      sta.pspkx = pspkx_mtx;
      sta.information = info_sta;

   else

      sta.x = [];
      sta.px = [];
      sta.pxspk = [];
      sta.pspk = [];
      sta.pspkx = [];
      sta.information = [];

   end

   info(i).sta = sta;




   % Process the rpdtest2 files
   %====================================================================


   if ( ~isempty(filestruct(i).rpdx1x2px_pxt_2) ) % make sure the 2d nonlinearity exists

      % mid1 filter
      [v1, coeff1, projection1, mtx1] = ...
         get_auditory_strf(filestruct(i).rpdtest2_v1, Nh, Nv, nlags);



      % mid2 filter
      [v2, coeff2, projection2, mtx2] = ...
         get_auditory_strf(filestruct(i).rpdtest2_v2, Nh, Nv, nlags);

% figure;
% subplot(1,3,1);
% imagesc(reshape(v_sta,Nh,nlags))
% set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% subplot(1,3,2);
% imagesc(reshape(v1,Nh,nlags))
% set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% subplot(1,3,3);
% imagesc(reshape(v2,Nh,nlags))
% set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% set(gcf,'position', [567   493   425   105]);
% pause


      % mid1 information
      [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, info_mid1] = ...
         get_mid1_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, ...
         Nbins, Nparts);

      mid1.x = x1_mtx;
      mid1.px = px1_mtx;
      mid1.pxspk = px1spk_mtx;
      mid1.pspk = pspk_mtx;
      mid1.pspkx = pspkx1_mtx;
      mid1.information = info_mid1;


      % mid2 information
      [x2_mtx, px2_mtx, px2spk_mtx, pspk_mtx, pspkx2_mtx, info_mid2] = ...
         get_mid2_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, ...
         Nbins, Nparts);


      mid2.x = x2_mtx;
      mid2.px = px2_mtx;
      mid2.pxspk = px2spk_mtx;
      mid2.pspk = pspk_mtx;
      mid2.pspkx = pspkx2_mtx;
      mid2.information = info_mid2;


      % mid1 and mid2 information
      [x1_mtx, x2_mtx, px1x2_mtx, px1x2spk_mtx, pspk_mtx, pspkx1x2_mtx, info_mid12] = ...
         get_mid1_mid2_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, ...
         coeff2, Nbins, Nparts);


      mid12.x1 = x1_mtx;
      mid12.x2 = x2_mtx;
      mid12.px1x2 = px1x2_mtx;
      mid12.px1x2spk = px1x2spk_mtx;
      mid12.pspk = pspk_mtx;
      mid12.pspkx1x2 = pspkx1x2_mtx;
      mid12.information = info_mid12;


   else % can't process files, so save as empty data

      % mid1 information
      mid1.x = [];
      mid1.px = [];
      mid1.pxspk = [];
      mid1.pspk = [];
      mid1.pspkx = [];
      mid1.information = [];

      % mid2 information
      mid2.x = [];
      mid2.px = [];
      mid2.pxspk = [];
      mid2.pspk = [];
      mid2.pspkx = [];
      mid2.information = [];

      % mid1 and mid2 information
      mid12.x1 = [];
      mid12.x2 = [];
      mid12.px1x2 = [];
      mid12.px1x2spk = [];
      mid12.pspk = [];
      mid12.pspkx1x2 = [];
      mid12.information = [];

   end

   % assign the data
   info(i).mid1 = mid1;
   info(i).mid2 = mid2;
   info(i).mid12 = mid12;

end % (for i)

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v_mean, coeff, projection, mtx] = get_auditory_strf(files, Nh, Nv, nlags)

mtx = [];

fsize = Nh*nlags;
Nn = fsize*Nv; % total number of elements in the STRF
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
v_mean = reshape(v_mean,fsize,1); % make it a column vector

k=1;

for i=1:Nparts
    for j=i+1:Nparts
        projection(k) = sum(mtx(:,i).*mtx(:,j)); % find every possible correlation between the 
        k=k+1;                                   % 4 estimated STAs?
    end
end

return;




function [x_mtx, px_mtx, pxspk_mtx, pspk_mtx, pspkx_mtx, information] = ...
   get_sta_information(files, coeff, Nbins, Nparts)
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density


if ( nargin < 2 | nargin > 4 )
   error('You need 2 to 4 input args.');
end

if ( nargin == 2 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 3 )
   Nparts = 4;
end

minpx = 0;


x_mtx = [];    
px_mtx = [];
pxspk_mtx = [];
pspkx_mtx = [];
pspk_mtx = [];



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


if (length(coeff)~=Nparts) 
    length(coeff)
    coeff
    Nparts
    display('wrong length of coeff vector');
    return
end

for i = 1:Nparts

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

return;




function [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, information] = ...
   get_mid1_information(files, coeff1, coeff2, Nbins, Nparts)

%plot_an_1d_from2dior_improved(fname,Nwiny,Nwinx,row,Nwinx,coeff_wnbest1,coeff_wnbest2,Nbins,Nparts,Nbins_short);
% summry of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the
%larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density


%       % mid1 information
%       [x1_mtx, px1_mtx, px1spk_mtx, pspk_mtx, pspkx1_mtx, info_mid1] = ...
%          get_mid1_information(filestruct(i).rpdx1x2px_pxt_2, coeff1, coeff2, ...
%          Nbins_medium, Nparts, Nbins_short);




% Variable name definitions and what they mean:
% 
% pspx1 = P(spike|x1)
% psp = P(spike)
% px1sp = P(x1|spike)
% px1 = P(x1)
% px1x2 = P(x1,x2)
% px1x2sp = P(x1,x2|spike)
% px1sp_psp = P(x1|sp) * P(spike)


if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx = 0;



% initialize variables
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




% Don't really this stuff
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
   get_mid2_information(files, coeff1, coeff2, Nbins, Nparts)
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

if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx=0;



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



if ( length(coeff1) ~= Nparts )
   length(coeff1)
   Nparts
   error('wrong length of coeff vector');
end

Nok = 0;

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

   if (coeff2(i)==1) 
      x2_r = x2; 
   end

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
      get_mid1_mid2_information(files, coeff1, coeff2, Nbins, Nparts)
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within 
% the larger bin, and dividing the two numbers, also introduced dx, and probability distribution 
% is now shown as a probability density1_mtx=[];
% also properly rotated the two-dimensional matrix, and its labels


if ( nargin < 3 | nargin > 5 )
   error('You need 3 to 5 input args.');
end

if ( nargin == 3 )
   Nbins = 21;
   Nparts = 4;
end

if ( nargin == 4 )
   Nparts = 4;
end

minpx = 0;

x1_mtx = [];
x2_mtx = [];

pspkx1x2_mtx = [];
pspk_mtx = [];
px1x2_mtx = [];
px1x2spk_mtx = [];
px1x2spk_pspk_mtx = [];

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


disp('mid2')


if ( length(coeff1) ~= Nparts )
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

pspk

    if ( coeff1(i)==1 )
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

   if ( coeff2(i) == -1 )
      for j=1:Nbins
         px1x2_r(:,Nbins+1-j) = px1x2(:,j);
         px1x2spk_pspk_r(:,Nbins+1-j) = px1x2spk_pspk(:,j);
         px1x2spk_r(:,Nbins+1-j) = px1x2spk(:,j);
         x2_r = -x2;
      end
   end

   if ( coeff2(i)==1 ) 
      x2_r = x2; 
   end

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





