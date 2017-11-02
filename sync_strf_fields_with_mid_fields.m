function newstr = sync_strf_fields_with_mid_fields(strf, middata, structtype)
%sync_strf_fields_with_mid_fields - align mid struct fields with strf fields
%
% newstr = sync_strf_fields_with_mid_fields(strf, middata, structtype)
%
% strf : strf data from which the original middata were originally extracted
%
% middata : struct array of data obtained from *.dat files from mid analysis
%
% structtype : the type of middata: can be: 'filestruct', 'filtstr', 'fio', 
% or 'proj'
%
% caa 1/27/10



if ( ~strcmp(structtype, 'filestruct') & ~strcmp(structtype, 'filtstr') & ...
~strcmp(structtype, 'fio') & ~strcmp(structtype, 'proj') )
   error('Unknown structure type.');
end


names = fieldnames(middata);

for i = 1:length(middata)

   unit = middata(i).unit;

   middata(i).exp = strf(unit).exp;
   middata(i).site = strf(unit).site;
   middata(i).chan = strf(unit).chan;
   middata(i).model = strf(unit).model;
   middata(i).depth = strf(unit).depth;
   middata(i).position = strf(unit).position;
   middata(i).stim = strf(unit).stim;
   middata(i).atten = strf(unit).atten;
   middata(i).spl = strf(unit).spl;
   middata(i).sm = strf(unit).sm;
   middata(i).tm = strf(unit).tm;
   middata(i).mdb = strf(unit).mdb;

end % (for i)



f = ['exp'; 'site'; 'chan'; 'model'; 'depth'; 'position'; 'stim'; ...
     'atten'; 'spl'; 'sm'; 'tm'; 'mdb'; names];

newstr = orderfields(middata, f);




return;