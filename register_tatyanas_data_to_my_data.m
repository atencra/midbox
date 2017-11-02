function [midnew, fionew] = register_tatyanas_data_to_my_data(mid, fio_params, strf)


if ( length(mid) ~= length(fio_params) | mid(1).location ~= fio_params(1).location )
   error('mid and fio_params must match.');
end


exp = strf(1).exp;
site = strf(1).site;
stim = strf(1).stim;

good_data = 0;

if ( mid(1).location==515 )
   if ( strcmp(exp, '2003-4-8') & site==15 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==516 )
   if ( strcmp(exp, '2002-8-26') & site==16 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==517 )
   if ( strcmp(exp, '2002-8-27') & site==17 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==519 )
   if ( strcmp(exp, '2002-8-27') & site==19 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==532 )
   if ( strcmp(exp, '2003-3-5') & site==32 & strcmp(stim,'dmr2') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==534 )
   if ( strcmp(exp, '2003-3-5') & site==34 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==537 )
   if ( strcmp(exp, '2003-3-5') & site==37 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==602 )
   if ( strcmp(exp, '2004-1-14') & site==2 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==603 )
   if ( strcmp(exp, '2004-1-14') & site==3 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==604 )
   if ( strcmp(exp, '2003-11-12') & site==4 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==608 )
   if ( strcmp(exp, '2004-1-14') & site==8 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==609 )
   if ( strcmp(exp, '2004-1-14') & site==9 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==611 )
   if ( strcmp(exp, '2004-1-14') & site==11 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==614 )
   if ( strcmp(exp, '2004-1-14') & site==14 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( mid(1).location==616 )
   if ( strcmp(exp, '2004-1-14') & site==16 & strcmp(stim,'dmr1') & length(mid)==length(strf) )
      good_data = 1;
   end
end

if ( ~good_data )
   error('MID and strf do not match up.');
end


for i = 1:length(mid)

   midnew(i).exp = strf(i).exp;
   midnew(i).site = strf(i).site;
   midnew(i).chan = strf(i).chan;
   midnew(i).model = strf(i).model;
   midnew(i).depth = strf(i).depth;
   midnew(i).position = strf(i).position;
   midnew(i).stim = strf(i).stim;
   midnew(i).atten = strf(i).atten;
   midnew(i).spl = strf(i).spl;
   midnew(i).sm = strf(i).sm;
   midnew(i).tm = strf(i).tm;
   midnew(i).mdb = strf(i).mdb;

   midnew(i).tbins = mid(i).tbins;
   midnew(i).fbins = mid(i).fbins;
   midnew(i).rpsta = mid(i).rpsta;
   midnew(i).rpx1pxpxt_sta = mid(i).rpx1pxpxt_sta;
   midnew(i).rpdtest1_v1 = mid(i).rpdtest1_v1;
   midnew(i).rpdtest1_v2 = mid(i).rpdtest1_v2;
   midnew(i).rpdx1x2px_pxt_1 = mid(i).rpdx1x2px_pxt_1;
   midnew(i).rpdtest2_v1 = mid(i).rpdtest2_v1;
   midnew(i).rpdtest2_v2 = mid(i).rpdtest2_v2;
   midnew(i).rpdx1x2px_pxt_2 = mid(i).rpdx1x2px_pxt_2;
   midnew(i).rpdbest1_v1 = mid(i).rpdbest1_v1;
   midnew(i).rpdbest1_v2 = mid(i).rpdbest1_v2;
   midnew(i).rpdbest2_v1 = mid(i).rpdbest2_v1;
   midnew(i).rpdbest2_v2 = mid(i).rpdbest2_v2;



   fionew(i).exp = strf(i).exp;
   fionew(i).site = strf(i).site;
   fionew(i).chan = strf(i).chan;
   fionew(i).model = strf(i).model;
   fionew(i).depth = strf(i).depth;
   fionew(i).position = strf(i).position;
   fionew(i).stim = strf(i).stim;
   fionew(i).atten = strf(i).atten;
   fionew(i).spl = strf(i).spl;
   fionew(i).sm = strf(i).sm;
   fionew(i).tm = strf(i).tm;
   fionew(i).mdb = strf(i).mdb;

   fionew(i).sta = fio_params(i).sta;
   fionew(i).v1 = fio_params(i).v1;
   fionew(i).v2 = fio_params(i).v2;
   fionew(i).r_sta_v1 = fio_params(i).r_sta_v1;
   fionew(i).p_sta_v1 = fio_params(i).p_sta_v1;
   fionew(i).r_sta_v2 = fio_params(i).r_sta_v2;
   fionew(i).p_sta_v2 = fio_params(i).p_sta_v2;
   fionew(i).r_v1_v2 = fio_params(i).r_v1_v2;
   fionew(i).p_v1_v2 = fio_params(i).p_v1_v2;


end % (for i)




