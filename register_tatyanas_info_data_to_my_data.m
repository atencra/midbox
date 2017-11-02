function [infonew, fionew] = register_tatyanas_data_to_my_data(info, strf)



exp = strf(1).exp;
site = strf(1).site;
stim = strf(1).stim;

good_data = 0;

if ( info(1).location==515 )
   if ( strcmp(exp, '2003-4-8') & site==15 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==516 )
   if ( strcmp(exp, '2002-8-26') & site==16 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==517 )
   if ( strcmp(exp, '2002-8-27') & site==17 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==519 )
   if ( strcmp(exp, '2002-8-27') & site==19 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==532 )
   if ( strcmp(exp, '2003-3-5') & site==32 & strcmp(stim,'dmr2') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==534 )
   if ( strcmp(exp, '2003-3-5') & site==34 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==537 )
   if ( strcmp(exp, '2003-3-5') & site==37 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==602 )
   if ( strcmp(exp, '2004-1-14') & site==2 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==603 )
   if ( strcmp(exp, '2004-1-14') & site==3 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==604 )
   if ( strcmp(exp, '2003-11-12') & site==4 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==608 )
   if ( strcmp(exp, '2004-1-14') & site==8 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==609 )
   if ( strcmp(exp, '2004-1-14') & site==9 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==611 )
   if ( strcmp(exp, '2004-1-14') & site==11 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==614 )
   if ( strcmp(exp, '2004-1-14') & site==14 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( info(1).location==616 )
   if ( strcmp(exp, '2004-1-14') & site==16 & strcmp(stim,'dmr1') & length(info)==length(strf) )
      good_data = 1;
   end
end

if ( ~good_data )
   error('info and strf do not match up.');
end


for i = 1:length(info)

   infonew(i).exp = strf(i).exp;
   infonew(i).site = strf(i).site;
   infonew(i).chan = strf(i).chan;
   infonew(i).model = strf(i).model;
   infonew(i).depth = strf(i).depth;
   infonew(i).position = strf(i).position;
   infonew(i).stim = strf(i).stim;
   infonew(i).atten = strf(i).atten;
   infonew(i).spl = strf(i).spl;
   infonew(i).sm = strf(i).sm;
   infonew(i).tm = strf(i).tm;
   infonew(i).mdb = strf(i).mdb;

   infonew(i).tbins = info(i).tbins;
   infonew(i).fbins = info(i).fbins;
   infonew(i).sta = info(i).sta;
   infonew(i).mid1 = info(i).mid1;
   infonew(i).mid2 = info(i).mid2;
   infonew(i).mid12 = info(i).mid12;


end % (for i)




