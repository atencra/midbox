function write_locator_to_isk_file(sprfile, spk, trigger)

% For the first try we can use the data for site 7, neuron #1


% load 2003-3-5-site32-2357um-30db-dmr2-fs24038-spkcmb.mat
% spkfile_part1='2003-3-5-site32-2357um-30db-dmr2-fs24038-spkcmb.mat';
% site_number=32;

%load 2002-8-26-site16-2414um-20db-dmr1-fs17857-spkcmb.mat
%load 2002-8-27-site17-2398um-20db-dmr1-fs17857-spkcmb.mat
%load 2002-8-27-site17-2398um-20db-dmr1-fs17857-spkcmb.mat


% load 2002-8-27-site19-2400um-20db-dmr1-fs17857-spkcmb.mat
% spkfile_part1='2002-8-27-site19-2400um-20db-dmr1-fs17857-spk';
% site_number=19;

% load 2003-3-5-site34-2371um-30db-dmr1-fs24038-spkcmb.mat
% spkfile_part1='2003-3-5-site34-2371um-30db-dmr1-fs24038-spk';
% site_number=34;
%
% load 2003-3-5-site37-2337um-30db-dmr1-fs24038-spkcmb.mat
% spkfile_part1='2003-3-5-site37-2371um-30db-dmr1-fs24038-spk';
% site_number=37;

% load 2003-4-8-site15-2332um-20db-dmr1-fs27173-spkcmb.mat
% spkfile_part1='2003-4-8-site15-2332um-20db-dmr1-fs27173-spk';
% site_number=15;

% load 2003-11-12-site4-2355um-30db-dmr1-fs20833-spkcmb.mat
% spkfile_part1='2003-11-12-site4-2355um-30db-dmr1-fs20833-spk';
% site_number=4;

% load 2004-1-14-site11-2391um-30db-dmr1-fs18115-spkcmb.mat
% spkfile_part1='2004-1-14-site11-2391um-30db-dmr1-fs18115-spkcmb.mat';
% site_number=11;

% load 2004-1-14-site14-2392um-30db-dmr1-fs18115-spkcmb.mat
% spkfile_part1='2004-1-14-site14-2392um-30db-dmr1-fs18115-spkcmb.mat';
% site_number=14;

% load 2004-1-14-site16-2390um-30db-dmr1-fs18115-spkcmb.mat
% spkfile_part1='2004-1-14-site16-2390um-30db-dmr1-fs18115-spkcmb.mat';
% site_number=16;
%
% load 2004-1-14-site2-2352um-30db-dmr1-fs18115-spkcmb.mat
% spkfile_part1='2004-1-14-site2-2352um-30db-dmr1-fs18115-spkcmb.mat';
% site_number=2;

% load 2004-1-14-site3-2388um-30db-dmr1-fs18115-spkcmb.mat
% spkfile_part1='2004-1-14-site3-2352um-30db-dmr1-fs18115-spkcmb.mat';
% site_number=3;

% load 2004-1-14-site8-2395um-30db-dmr1-fs18115-spkcmb.mat
% spkfile_part1='2004-1-14-site8-2395um-30db-dmr1-fs18115-spkcmb.mat';
% site_number=8;

%load 2004-1-14-site9-2389um-30db-dmr1-fs18115-spkcmb.mat
%spkfile_part1='2004-1-14-site9-2389um-30db-dmr1-fs18115-spkcmb.mat';
%site_number=9;
%fss=18115;
%fss=20833;
%fss=24038;
%fss=17857;
% fss=24038;
%fss=27173;
%first stimulus file:
% specfile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min-dnsmp.spr'
%second stimulu file:
%specfile='dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min_12_carriers_per_octave.spr'
%specfile='dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min_6_carriers_per_octave.spr'
%fid = fopen(specfile);
%[a]=fread(fid,'inf','float');
%fclose(fid);
%fid=fopen('dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min-dnsmp.raw','wb');
%fwrite(fid,a,'double');
%fclose(fid);



% sprfile = 'dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8.spr';


t1 = 0;
t2 = 0.1;


for i = 1:length(spk) %Ncells

	fsspk = spk(i).fs;
	spet = spk(i).spiketimes / 1000 * fsspk;

	[t, f, locator, ns, ar] = get_locator_for_mid_analysis(sprfile, t1, t2, spet, trigger, fsspk);



% 	[tax,fax,rf,n0,averate,locator_Craig]=get_downsampled_strf(specfile,t1,t2,spet,trigger,fss) ;
% 	if (i>20)
% 		figure(2)
% 		subplot(5,4,i-20)
% 	else
% 		figure (1)
% 		subplot(5,4,i)
% 	end
% 
%    imagesc(tax,fax,rf)
%    title(sprintf('cell %u site %u Nspikes=%u',i,site_number,n0))
%    %  imagesc(rf)
%    [taxis, faxis, locator, numspikes, averate]=get_locator(specfile, t1, t2, spet, trigger, fss);
%    %n0
%    ts=size(taxis)
%    fs=size(faxis)


	exp = spk(i).exp;
	site = spk(i).site;
	depth = spk(i).depth;
	atten = spk(i).atten;
	stim = spk(i).stim;


%     spkfile=sprintf('%s-%u.isk',spkfile_part1,i);

	spkfile = sprintf('%s-site%.0f-%.0fum-%.0fdb-%s-fs%.0f-spk-%u.isk', ...
		exp, site, depth, atten, stim, fsspk, i)

%     spkfile = sprintf('2002-8-26-site16-2414um-20db-dmr1-fs17857-spk-%u.isk',i);

    %   spkfile=sprintf('2003-3-5-site32-2357um-30db-dmr2-fs24038-spk-%u.isk',i);


% pause

	fid = fopen( spkfile, 'w' );
	count = fprintf( fid, '%u \n', locator );
	fclose( fid );

   fprintf( '\nNtrials, or Movie_length, = %.0f\n', length(locator) );

end



%numer of spikes
% interasting 11,12,15,16 - lots of spikes and very little structure
%[2536;9663;3869;1836;4144;1486;2722;2433;6525;21962;20839;18857;2310;1661;
%15088;30671;8045;969;546]



