function [projection, location] = get_filter_projection(strf, SpecFile, T1, T2, spet, Trig, Fss, SPL, MdB)
%
%       FILE NAME       : get_filter_projection
%       DESCRIPTION     : Find projection of ripple stimulus segments onto
%                          the strf
%
%  strf : filter, with low-high frequencies in first-last rows
%         and time from high-low in first-last columns
%	SpecFile	: Spectral Profile File (*.spr file)
%	T1, T2		: Evaluation delay interval for WSTRF(T,F)
%			  T E [- T1 , T2 ], Note that T1 and T2 > 0
%	spet		: Array of spike event times in sample number
%	Trig		: Array of Trigger Times
%	Fss		: Sampling Rate for TRIGGER and SPET
%	SPL		: Signal RMS Sound Pressure Level
%	MdB		: Signal Modulation Index in dB
%
%
%	RETURNED VALUES 
%
%	projection : Vector of projection values, 1 value for each spike
%	spindex	: Spet indecees for channel 1
%
% To get things to match up with Tatyana's code you use the following line
% to compute the STA:
%
% [projection, spindex] = get_filter_projection(strf, file, 0.005, 0.095, spet, trigger, fs, SPL, MdB);
%
% 
% caa 12/18/06


Sound = 'MR';
ModType = 'dB';
sprtype='float';

NTrig = length(Trig);

%Loading Parameter Data
index = findstr(SpecFile,'.spr');
ParamFile = [SpecFile(1:index(1)-1) '_param.mat'];
f = ['load ' ParamFile];
eval(f);

clear App MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 f2 
clear Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn

%Fliping Trig and Spet for channel 2 and channel 1 STRFs
MinTime = min([Trig spet]);
MaxTime = max([Trig spet]);
spet = spet-MinTime+1;
Trig = Trig-MinTime+1;

%Converting Temporal Delays to Sample Numbers
N1 = round(T1*Fs/DF);
N2 = round(T2*Fs/DF);

%Opening Spectral Profile File
fid = fopen(SpecFile);

%Initializing Some Variables
TrigCount = 2;
tbins = size(strf,2);
fbins = size(strf,1);
spindex = []; % Spike Index
projection = [];
location = [];

%Fiding Mean Spectral Profile and RMS Power
if strcmp(Sound,'RN')
	RMSP=-MdB/2;		% RMS value of normalized Spectral Profile
	PP=MdB^2/12;		% Modulation Depth Variance 
elseif strcmp(Sound,'MR')
	RMSP=-MdB/2;		% RMS value of normalized Spectral Profile
	PP=MdB^2/8;		% Modulation Depth Variance 
end

%Initializing First and Second Spectral Profile Segments
frewind(fid);
if strcmp(sprtype,'float')
	S1=fread(fid,NT*NF,'float');
	S2=fread(fid,NT*NF,'float');
	S3=fread(fid,NT*NF,'float');
else
	S1=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
	S2=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
	S3=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
end

S1 = reshape(S1,NF,NT);
S2 = reshape(S2,NF,NT);
S3 = reshape(S3,NF,NT);


%Computing Spectro Temporal Receptive Field Variability - 'dB'
while ( ~feof(fid) & TrigCount<length(Trig)-1 )

	%Finding SPET in between triggers
	index = find(spet>=Trig(TrigCount) & spet<Trig(TrigCount+1));

	%Spike Indeces Used for Computing Projection Values
	spindex = [spindex index];

	%Resampling spet relative to the Spectral Profile samples
	spettrig = ceil( (spet(index)-Trig(TrigCount)+1) * Fs / Fss /DF );

	%Finding Receptive Field Variability for Channel 1
	epsilon = 10^(-MdB/20);

   for k = 1:length(spettrig)

      %Setting Spike Time and STRF length
      M = size(S1,2); % # of time bins
      L = spettrig(k);

      if ( L <= M )  %Condition to circumvent rounding errors 

         % Finding Pre-Event Spectral Profiles
			if ( L < N2 )
			   segment = MdB*[S1(:,M-(N2-L-1):M) S2(:,1:L+N1)] - RMSP;
            segment = segment(end-(fbins-1):end,:);
 			   %segment1 = segment(end-(fbins-1):end,end-(tbins-1):end);
			elseif ( L+N1 > M )
            segment = MdB*[S2(:,L-N2+1:M) S3(:,1:N1-M+L)] - RMSP;
            segment = segment(end-(fbins-1):end,:);
            %segment1 = segment(end-(fbins-1):end,end-(tbins-1):end);
			else
            segment = MdB*[S2(:,L-N2+1:L+N1)] - RMSP;
            segment = segment(end-(fbins-1):end,:);
            %segment1 = segment(end-(fbins-1):end,end-(tbins-1):end);
         end % (if)


         % Now downsample the segment using Tatyana's method
         flipsegment = fliplr(segment); % flip for easier indexing, t=0 is now in column 1
         stimdnsmp = zeros(fbins, tbins);
         for i = 1:tbins
            stimdnsmp(:,i) = ( flipsegment(:,2*i-1) + flipsegment(:,2*i) ) / 2;
         end
         stimdnsmp = fliplr(stimdnsmp); % flip back to normal orientation, t=0 is now in last column


% size(segment)
% size(stimdnsmp)
% size(strf)
% 
% subplot(2,1,1);
% imagesc(segment(:,end-2*tbins+1:end));
% subplot(2,1,2);
% imagesc(stimdnsmp);
% 
% pause;

         if ( size(stimdnsmp,1)==size(strf,1) & size(stimdnsmp,2)==size(strf,2) )
            stimvec = reshape(stimdnsmp, 1, size(stimdnsmp,1)*size(stimdnsmp,2)); % 1 x N
            strfvec = reshape(strf, size(strf,1)*size(strf,2), 1); % N x 1
            projection = [projection; stimvec*strfvec];
            location = [location; TrigCount]; % Keep track of which trigger corresponds to which projection.
                                              % Will be useful for later information analysis.
         else
            error('stim and strf dimensions do not match.');
         end

%          size(stimvec)
%          size(strfvec)
%          pause

      end % (if L<=M)

   end % (for k)

	%Reading Spectral Profile Data File
	S1 = S2;
	S2 = S3;

	if ( strcmp(sprtype,'float') )
		S3 = fread(fid,NT*NF,'float');
	else
		S3 = fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
	end

	if ( ~feof(fid) )
		S3 = reshape(S3, NF, NT);
	end

	%Updating Trigger Counter
	TrigCount = TrigCount + 1;
   if ( mod(TrigCount,250) == 0 )
      disp(['Block Number ' num2str(TrigCount) ' of ' num2str(NTrig)])
   end

%    %Displaying Output
%    if TrigCount/NBlocks==round(TrigCount/NBlocks)
%       if length(p1)>1
%          subplot(211)
%          hist(p1,-1:.1:1)
% 			xlabel('Correlation Coefficient')
% 			ylabel('Counts')
%          pause(0)
%       end
% 
%       if length(p2)>1
%          subplot(212)
%          hist(p2,-1:.1:1)
% 			xlabel('Correlation Coefficient')
% 			ylabel('Counts')
%          pause(0)
%       end
%    end

end % (while)

%Closing all opened files
fclose('all');




