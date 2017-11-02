function locator = locator_from_isk_file(iskfile)
% locator_from_isk_file - Takes an iskfile, returns a spike train vector.
%
% locator = locator_from_isk_file(iskfile)
% ------------------------------------------------------
%
% The vector is composed of either 0s or 1s. Each element corresponds to
% a particular trial, or time bin. 
%
% The iskfile is created using the function
%
% get_locator_for_mid_analysis.m
%
% The isk files are used in the MID analysis.
%
% caa 3/11/09

fid = fopen(iskfile, 'r');

locator = fscanf(fid, '%d');

fclose(fid);

[nr,nc] = size(locator);

if ( nr > nc )
   locator = locator';
end

return