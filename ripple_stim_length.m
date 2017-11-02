%
% [n, NT, NF] = ripple_stim_length(specfile)
%
% Determine length of ripple stimulus - How many time trials were there? 
%
%  INPUT VALUES
%
%  specfile       :  Spectral Profile File ( *.spr file )
%
%
%  RETURNED VALUES 
%
%  n : Total number of envelope samples. To get number of trials, or time
%      bins, divide n by NF, n / NF.
%  NT : Number of time samples in a ripple envelope segment
%  NF : Number of freq samples in a ripple segment
%
%
function [n, NT, NF] = ripple_stim_length(specfile)


% Loading Parameter Data
index = findstr(specfile,'.spr');
paramfile = [specfile(1:index(1)-1) '_param.mat'];
f = ['load ' paramfile];
eval(f);

% Clear possible useless variables in the param file
clear App MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 
clear f2 Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn


fid = fopen(specfile); % open ripple envelop file
frewind(fid); % make sure the file is set to the beginning


n = 0; 

while ( ~feof(fid) )

      [env, count] = fread(fid, NT * NF, 'float');

      n = n + count;

end % ( while ~feof(fid) )


% Closing all opened files

fclose('all');





