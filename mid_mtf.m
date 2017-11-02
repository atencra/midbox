function [tmf1s, tmtf, xmf1s, xmtf, rtf1s] = mid_mtf(tmf, xmf, rtf)
%MID_MTF - Calculate MTFs from a ripple transfer function
%
% [tmf1s, tmtf, xmf1s, xmtf] = mid_mtf(tmf, xmf, rtf)
%
% tmf : temporal modulation frequency of the rtf
% xmf : spectral modulation frequency of the rtf
%
% rtf : ripple transfer function. Goes from 0 to Max(xmf) on the spectral
% modulation axis and from negative to positive tmf values.
%
% tmf1s : temporal modulation frequency from 0 to max(tmf)
% tmtf : temporal modulation transfer function
%
% xmf1s : spectral modulation freq's from 0 to max(xmf)
% xmtf : spectral modulation transfer function
%
% caa

% temporal and spectral modulation frequency axes

tmf = tmf(:)';
rtf_tmf = tmf;
ind0 = find(tmf==0);
tmf1s = tmf(ind0:end);

xmf1s = xmf(:)';


% the ripple transfer function and its folded version
rtf_right = rtf(:,ind0:end);
rtf_left = rtf(:,1:ind0-1);
rtfnew = zeros(size(rtf_right));
rtfnew = rtfnew + rtf_right;
rtfnew(:,2:end) = rtfnew(:,2:end) + fliplr(rtf_left);
rtfnew(:,1) = 2*rtfnew(:,1);
rtf1s = rtfnew;


% temporal modulation transfer function
tmtf = sum(rtf,1);
tmtf = tmtf(:)';
tmtf_left = tmtf(1:ind0-1);
tmtf_right = tmtf(ind0:end);
tmtf = tmtf_right + [0 fliplr(tmtf_left)];
tmtf(1) = 2 * tmtf(1);
tmtf = tmtf ./ max(tmtf);


% spectral modulation transfer function
xmtf = sum(rtf,2)';
xmtf = xmtf(:)';
xmtf = xmtf ./ max(xmtf);


return;

