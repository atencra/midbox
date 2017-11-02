function spkloc = get_spkloc_x0_dimx_nf_nlags_cy(spkloc)
% get_spkloc_from_isk_spk  Add locator to spk struct array 
% 
%    spkloc = get_spkloc_from_isk_spk(spk) reads in .isk files for
%    the corresponding elements in spk and returns a struct away
%    with two new fields: iskfile and locator. 
% 
%    iskfile is a binary file that holds the spike train of the 
%    neuron
% 
%    locator is a vector of 0s and 1s, corresponding to the spiketimes
% 
%    The data that is coordinated was used by Tatyana Sharpee to estimate
%    MIDs in Atencio et al (2008) and Atencio et al (2009).
% 
%    spkloc now has the locator vector. To estimate an STA the
%    stimulus matrix may be read in and correlated with the
%    spike train
% 
%    Craig Atencio
%    4/10/12

error(nargchk(1,1,nargin));

exp = spkloc(1).exp;

if ( strcmp(exp, '2002-8-26') || strcmp(exp, '2002-8-27') || ...
     strcmp(exp, '2003-3-5') || strcmp(exp, '2003-4-8') )

   stimfile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min_6_carriers_per_octave-matrix.mat';
   dimx = 33; 
   nf = 25;
   nlags = 20;
   x0 = dimx - nf + 1;
   cy = 2;

elseif ( strcmp(exp, '2003-11-12') || strcmp(exp, '2004-1-14') )

   stimfile = 'dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min_6_carriers_per_octave-matrix.mat';
   dimx = 40;
   nf = 30;
   nlags = 20;
   x0 = dimx - nf + 1;
   cy = 2;

else

   error('Experiment in spkloc does not match mid experiment list.');

end % (if -elseif)


for i = 1:length(spkloc)
   spkloc(i).stimfile = stimfile;
   spkloc(i).dimx = dimx;
   spkloc(i).nf = nf;
   spkloc(i).nlags = nlags;
   spkloc(i).x0 = x0;
   spkloc(i).cy = cy;
end


return;







% 
% 
%   switch (exp){
%   case 516:
%     sprintf(spk_file,"2002-8-26-site16-2414um-20db-dmr1-fs17857-%u.isk",global_cell);
%     //    Nh=12; dimx=33;    x0=dimx-Nh+1;   nlags=32;
%     signl=501;
%     break;
%   case 534:
%     sprintf(spk_file,"2003-3-5-site34-2371um-30db-dmr1-fs24038-%u.isk",global_cell);
%     signl=501;
%     break;
%   case 537:
%     sprintf(spk_file,"2003-3-5-site37-2371um-30db-dmr1-fs24038-%u.isk",global_cell);
%     signl=501;
%     break;
%   case 515:
%     sprintf(spk_file,"2003-4-8-site15-2332um-20db-dmr1-fs27173-%u.isk",global_cell);
%     signl=501;
%     cout<<spk_file<<endl;
%     break;
%   case 519:
%     sprintf(spk_file,"2002-8-27-site19-2400um-20db-dmr1-fs17857-%u.isk",global_cell);
%     signl=501;
%     break;
% 
% 
% 
% 
%   case 609:
%     sprintf(spk_file,"2004-1-14-site9-2389um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 602:
%     sprintf(spk_file,"2004-1-14-site2-2352um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 603:
%     sprintf(spk_file,"2004-1-14-site3-2388um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 608:
%     sprintf(spk_file,"2004-1-14-site8-2395um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 614:
%     sprintf(spk_file,"2004-1-14-site14-2392um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 616:
%     sprintf(spk_file,"2004-1-14-site16-2390um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 611:
%     sprintf(spk_file,"2004-1-14-site11-2391um-30db-dmr1-fs18115-%u.isk",global_cell);
%     signl=502;
%     break;
%   case 604:
%     sprintf(spk_file,"2003-11-12-site4-2355um-30db-dmr1-fs20833-%u.isk",global_cell);
%     signl=502;
%     break;
% 
%   case 532:  
%     sprintf(spk_file,"2003-3-5-site32-2357um-30db-dmr2-fs24038-%u.isk",global_cell);
%     //nlags=24;Nh=15;    dimx=33;    x0=dimx-Nh+1;
%     signl=501;
%     break;
%   case 517: 
%     sprintf(spk_file,"2002-8-27-site17-2398um-20db-dmr1-fs17857-%u.isk",global_cell);
%     //Nh=12;    x0=18;    dimx=33;nlags=32;
%     signl=501;
%     break;
% 
%   
% 
%   if (signl==501) {
%     dimx=33; 
%     Nh=25;
%     cy=2;
%     x0=dimx-Nh+1;
%     Movie_length=320*1240/cy;
%     Movie_length=(Movie_length>>2)*4;
%   }
% 
% 
% 
%   if ( signl==502) {
%     dimx=40;
%     Nh=30;
%     x0=dimx-Nh+1;
%     cy=2;
%     Movie_length=3779*(134/cy);
%     Movie_length=(Movie_length>>2)*4;
%   }
% 
% 
%   nlags=20;
% 
























