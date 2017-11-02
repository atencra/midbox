function infodata = plot_proj_info(info)
%
%
% caa 12/20/06


% Get information for the STA:
datalength = info(1).datalength;
numbins = info(1).numbins;

for k = 1:length(info)

   cell_num = k;
   exp = info(k).exp;
   site = info(k).site;
   chan = info(k).chan;
   model = info(k).model;

   % The data will be organized into matrices. Each row of the matrix
   % corresponds to a bin size. Each column corresponds to a
   % datafraction

   % STA - Part 1:
   ista_part1_test = zeros( length(numbins), length(datalength) );
   ista_part1_train = zeros( length(numbins), length(datalength) );

   % STA - Part 2:
   ista_part2_test = zeros( length(numbins), length(datalength) );
   ista_part2_train = zeros( length(numbins), length(datalength) );

   % STA - Part 3:
   ista_part3_test = zeros( length(numbins), length(datalength) );
   ista_part3_train = zeros( length(numbins), length(datalength) );

   % STA - Part 4:
   ista_part4_test = zeros( length(numbins), length(datalength) );
   ista_part4_train = zeros( length(numbins), length(datalength) );


   % MID - Part 1:
   imid1_part1_test = zeros( length(numbins), length(datalength) );
   imid2_part1_test = zeros( length(numbins), length(datalength) );
   imid12_part1_test = zeros( length(numbins), length(datalength) );

   imid1_part1_train = zeros( length(numbins), length(datalength) );
   imid2_part1_train = zeros( length(numbins), length(datalength) );
   imid12_part1_train = zeros( length(numbins), length(datalength) );


   % MID - Part 2:
   imid1_part2_test = zeros( length(numbins), length(datalength) );
   imid2_part2_test = zeros( length(numbins), length(datalength) );
   imid12_part2_test = zeros( length(numbins), length(datalength) );

   imid1_part2_train = zeros( length(numbins), length(datalength) );
   imid2_part2_train = zeros( length(numbins), length(datalength) );
   imid12_part2_train = zeros( length(numbins), length(datalength) );


   % MID - Part 3:
   imid1_part3_test = zeros( length(numbins), length(datalength) );
   imid2_part3_test = zeros( length(numbins), length(datalength) );
   imid12_part3_test = zeros( length(numbins), length(datalength) );

   imid1_part3_train = zeros( length(numbins), length(datalength) );
   imid2_part3_train = zeros( length(numbins), length(datalength) );
   imid12_part3_train = zeros( length(numbins), length(datalength) );


   % MID - Part 4:
   imid1_part4_test = zeros( length(numbins), length(datalength) );
   imid2_part4_test = zeros( length(numbins), length(datalength) );
   imid12_part4_test = zeros( length(numbins), length(datalength) );

   imid1_part4_train = zeros( length(numbins), length(datalength) );
   imid2_part4_train = zeros( length(numbins), length(datalength) );
   imid12_part4_train = zeros( length(numbins), length(datalength) );


   for i = 1:length(numbins)

      for j = 1:length(datalength)

         % Assign STA info Part 1-4:
         ista_part1_test(i, j)  = info(k).sta_part1_test{j}{i};
         ista_part1_train(i, j) = info(k).sta_part1_train{j}{i};

         ista_part2_test(i, j)  = info(k).sta_part2_test{j}{i};
         ista_part2_train(i, j) = info(k).sta_part2_train{j}{i};

         ista_part3_test(i, j)  = info(k).sta_part3_test{j}{i};
         ista_part3_train(i, j) = info(k).sta_part3_train{j}{i};

         ista_part4_test(i, j)  = info(k).sta_part4_test{j}{i};
         ista_part4_train(i, j) = info(k).sta_part4_train{j}{i};


         % Assign MID info Part 1:
         imid12_part1_test(i, j) = info(k).mid1mid2_part1_test{j}{i}(1);
         imid1_part1_test(i, j)  = info(k).mid1mid2_part1_test{j}{i}(2);
         imid2_part1_test(i, j)  = info(k).mid1mid2_part1_test{j}{i}(3);

         imid12_part1_train(i, j) = info(k).mid1mid2_part1_train{j}{i}(1);
         imid1_part1_train(i, j)  = info(k).mid1mid2_part1_train{j}{i}(2);
         imid2_part1_train(i, j)  = info(k).mid1mid2_part1_train{j}{i}(3);


         % Assign MID info Part 2:
         imid12_part2_test(i, j) = info(k).mid1mid2_part2_test{j}{i}(1);
         imid1_part2_test(i, j)  = info(k).mid1mid2_part2_test{j}{i}(2);
         imid2_part2_test(i, j)  = info(k).mid1mid2_part2_test{j}{i}(3);

         imid12_part2_train(i, j) = info(k).mid1mid2_part2_train{j}{i}(1);
         imid1_part2_train(i, j)  = info(k).mid1mid2_part2_train{j}{i}(2);
         imid2_part2_train(i, j)  = info(k).mid1mid2_part2_train{j}{i}(3);


         % Assign MID info Part 3:
         imid12_part3_test(i, j) = info(k).mid1mid2_part3_test{j}{i}(1);
         imid1_part3_test(i, j)  = info(k).mid1mid2_part3_test{j}{i}(2);
         imid2_part3_test(i, j)  = info(k).mid1mid2_part3_test{j}{i}(3);

         imid12_part3_train(i, j) = info(k).mid1mid2_part3_train{j}{i}(1);
         imid1_part3_train(i, j)  = info(k).mid1mid2_part3_train{j}{i}(2);
         imid2_part3_train(i, j)  = info(k).mid1mid2_part3_train{j}{i}(3);


         % Assign MID info Part 4:
         imid12_part4_test(i, j) = info(k).mid1mid2_part4_test{j}{i}(1);
         imid1_part4_test(i, j)  = info(k).mid1mid2_part4_test{j}{i}(2);
         imid2_part4_test(i, j)  = info(k).mid1mid2_part4_test{j}{i}(3);

         imid12_part4_train(i, j) = info(k).mid1mid2_part4_train{j}{i}(1);
         imid1_part4_train(i, j)  = info(k).mid1mid2_part4_train{j}{i}(2);
         imid2_part4_train(i, j)  = info(k).mid1mid2_part4_train{j}{i}(3);


      end % (for j)

   end % (for i)



   % Get the information data for the training data

   [ista_extrap_part1_train, imid1_extrap_part1_train, imid2_extrap_part1_train, imid12_extrap_part1_train, ...
    ista_extrap_normr_part1_train, imid1_extrap_normr_part1_train, imid2_extrap_normr_part1_train, imid12_extrap_normr_part1_train] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part1_train, imid1_part1_train, imid2_part1_train, imid12_part1_train, datalength, numbins);


   [ista_extrap_part2_train, imid1_extrap_part2_train, imid2_extrap_part2_train, imid12_extrap_part2_train, ...
    ista_extrap_normr_part2_train, imid1_extrap_normr_part2_train, imid2_extrap_normr_part2_train, imid12_extrap_normr_part2_train] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part2_train, imid1_part2_train, imid2_part2_train, imid12_part2_train, datalength, numbins);


   [ista_extrap_part3_train, imid1_extrap_part3_train, imid2_extrap_part3_train, imid12_extrap_part3_train, ...
    ista_extrap_normr_part3_train, imid1_extrap_normr_part3_train, imid2_extrap_normr_part3_train, imid12_extrap_normr_part3_train] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part3_train, imid1_part3_train, imid2_part3_train, imid12_part3_train, datalength, numbins);


   [ista_extrap_part4_train, imid1_extrap_part4_train, imid2_extrap_part4_train, imid12_extrap_part4_train, ...
    ista_extrap_normr_part4_train, imid1_extrap_normr_part4_train, imid2_extrap_normr_part4_train, imid12_extrap_normr_part4_train] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part4_train, imid1_part4_train, imid2_part4_train, imid12_part4_train, datalength, numbins);



%    [istainf_part1_train, imid1inf_part1_train, imid2inf_part1_train, imid12inf_part1_train, ...
%     istainf_normr_part1_train, imid1inf_normr_part1_train, imid2inf_normr_part1_train, imid12inf_normr_part1_train] = ...
%       get_info_extrapolation(ista_part1_train, imid1_part1_train, imid2_part1_train, imid12_part1_train, datalength, numbins);
% 
%    [istainf_part2_train, imid1inf_part2_train, imid2inf_part2_train, imid12inf_part2_train, ...
%     istainf_normr_part2_train, imid1inf_normr_part2_train, imid2inf_normr_part2_train, imid12inf_normr_part2_train] = ...
%       get_info_extrapolation(ista_part2_train, imid1_part2_train, imid2_part2_train, imid12_part2_train, datalength, numbins);
% 
%    [istainf_part3_train, imid1inf_part3_train, imid2inf_part3_train, imid12inf_part3_train, ...
%     istainf_normr_part3_train, imid1inf_normr_part3_train, imid2inf_normr_part3_train, imid12inf_normr_part3_train] = ...
%       get_info_extrapolation(ista_part3_train, imid1_part3_train, imid2_part3_train, imid12_part3_train, datalength, numbins);
% 
%    [istainf_part4_train, imid1inf_part4_train, imid2inf_part4_train, imid12inf_part4_train, ...
%     istainf_normr_part4_train, imid1inf_normr_part4_train, imid2inf_normr_part4_train, imid12inf_normr_part4_train] = ...
%       get_info_extrapolation(ista_part4_train, imid1_part4_train, imid2_part4_train, imid12_part4_train, datalength, numbins);




   % Get the information data for the test data

   [ista_extrap_part1_test, imid1_extrap_part1_test, imid2_extrap_part1_test, imid12_extrap_part1_test, ...
    ista_extrap_normr_part1_test, imid1_extrap_normr_part1_test, imid2_extrap_normr_part1_test, imid12_extrap_normr_part1_test] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part1_test, imid1_part1_test, imid2_part1_test, imid12_part1_test, datalength, numbins);


   [ista_extrap_part2_test, imid1_extrap_part2_test, imid2_extrap_part2_test, imid12_extrap_part2_test, ...
    ista_extrap_normr_part2_test, imid1_extrap_normr_part2_test, imid2_extrap_normr_part2_test, imid12_extrap_normr_part2_test] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part2_test, imid1_part2_test, imid2_part2_test, imid12_part2_test, datalength, numbins);


   [ista_extrap_part3_test, imid1_extrap_part3_test, imid2_extrap_part3_test, imid12_extrap_part3_test, ...
    ista_extrap_normr_part3_test, imid1_extrap_normr_part3_test, imid2_extrap_normr_part3_test, imid12_extrap_normr_part3_test] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part3_test, imid1_part3_test, imid2_part3_test, imid12_part3_test, datalength, numbins);


   [ista_extrap_part4_test, imid1_extrap_part4_test, imid2_extrap_part4_test, imid12_extrap_part4_test, ...
    ista_extrap_normr_part4_test, imid1_extrap_normr_part4_test, imid2_extrap_normr_part4_test, imid12_extrap_normr_part4_test] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista_part4_test, imid1_part4_test, imid2_part4_test, imid12_part4_test, datalength, numbins);


%    [istainf_part1_test, imid1inf_part1_test, imid2inf_part1_test, imid12inf_part1_test, ...
%     istainf_normr_part1_test, imid1inf_normr_part1_test, imid2inf_normr_part1_test, imid12inf_normr_part1_test] = ...
%       get_info_extrapolation(ista_part1_test, imid1_part1_test, imid2_part1_test, imid12_part1_test, datalength, numbins);
% 
%    [istainf_part2_test, imid1inf_part2_test, imid2inf_part2_test, imid12inf_part2_test, ...
%     istainf_normr_part2_test, imid1inf_normr_part2_test, imid2inf_normr_part2_test, imid12inf_normr_part2_test] = ...
%       get_info_extrapolation(ista_part2_test, imid1_part2_test, imid2_part2_test, imid12_part2_test, datalength, numbins);
% 
%    [istainf_part3_test, imid1inf_part3_test, imid2inf_part3_test, imid12inf_part3_test, ...
%     istainf_normr_part3_test, imid1inf_normr_part3_test, imid2inf_normr_part3_test, imid12inf_normr_part3_test] = ...
%       get_info_extrapolation(ista_part3_test, imid1_part3_test, imid2_part3_test, imid12_part3_test, datalength, numbins);
% 
%    [istainf_part4_test, imid1inf_part4_test, imid2inf_part4_test, imid12inf_part4_test, ...
%     istainf_normr_part4_test, imid1inf_normr_part4_test, imid2inf_normr_part4_test, imid12inf_normr_part4_test] = ...
%       get_info_extrapolation(ista_part4_test, imid1_part4_test, imid2_part4_test, imid12_part4_test, datalength, numbins);


   % Save the data in the struct array proj
   % ------------------------------------------------
   infodata(k).exp = info(k).exp;
   infodata(k).site = info(k).site;
   infodata(k).chan = info(k).chan;
   infodata(k).model = info(k).model;
   infodata(k).depth = info(k).depth;
   infodata(k).position = info(k).position;
   infodata(k).stim = info(k).stim;
   infodata(k).atten = info(k).atten;
   infodata(k).spl = info(k).spl;
   infodata(k).sm = info(k).sm;
   infodata(k).tm = info(k).tm;
   infodata(k).mdb = info(k).mdb;
   infodata(k).datalength = info(k).datalength;
   infodata(k).numbins = info(k).numbins;
   infodata(k).numreps = info(k).numreps;


   % STA Parts 1-4:

   infodata(k).sta.part1_train = ista_part1_train;
   infodata(k).sta.part2_train = ista_part2_train;
   infodata(k).sta.part3_train = ista_part3_train;
   infodata(k).sta.part4_train = ista_part4_train;

   infodata(k).sta.part1_train_extrap = ista_extrap_part1_train;
   infodata(k).sta.part2_train_extrap = ista_extrap_part2_train;
   infodata(k).sta.part3_train_extrap = ista_extrap_part3_train;
   infodata(k).sta.part4_train_extrap = ista_extrap_part4_train;

   infodata(k).sta.part1_normr_train_extrap = ista_extrap_normr_part1_train;
   infodata(k).sta.part2_normr_train_extrap = ista_extrap_normr_part2_train;
   infodata(k).sta.part3_normr_train_extrap = ista_extrap_normr_part3_train;
   infodata(k).sta.part4_normr_train_extrap = ista_extrap_normr_part4_train;

   infodata(k).sta.part1_test = ista_part1_test;
   infodata(k).sta.part2_test = ista_part2_test;
   infodata(k).sta.part3_test = ista_part3_test;
   infodata(k).sta.part4_test = ista_part4_test;

   infodata(k).sta.part1_test_extrap = ista_extrap_part1_test;
   infodata(k).sta.part2_test_extrap = ista_extrap_part2_test;
   infodata(k).sta.part3_test_extrap = ista_extrap_part3_test;
   infodata(k).sta.part4_test_extrap = ista_extrap_part4_test;

   infodata(k).sta.part1_normr_test_extrap = ista_extrap_normr_part1_test;
   infodata(k).sta.part2_normr_test_extrap = ista_extrap_normr_part2_test;
   infodata(k).sta.part3_normr_test_extrap = ista_extrap_normr_part3_test;
   infodata(k).sta.part4_normr_test_extrap = ista_extrap_normr_part4_test;



   % First MID Parts 1-4:

   infodata(k).mid1.part1_train = imid1_part1_train;
   infodata(k).mid1.part2_train = imid1_part2_train;
   infodata(k).mid1.part3_train = imid1_part3_train;
   infodata(k).mid1.part4_train = imid1_part4_train;

   infodata(k).mid1.part1_train_extrap = imid1_extrap_part1_train;
   infodata(k).mid1.part2_train_extrap = imid1_extrap_part2_train;
   infodata(k).mid1.part3_train_extrap = imid1_extrap_part3_train;
   infodata(k).mid1.part4_train_extrap = imid1_extrap_part4_train;

   infodata(k).mid1.part1_normr_train_extrap = imid1_extrap_normr_part1_train;
   infodata(k).mid1.part2_normr_train_extrap = imid1_extrap_normr_part2_train;
   infodata(k).mid1.part3_normr_train_extrap = imid1_extrap_normr_part3_train;
   infodata(k).mid1.part4_normr_train_extrap = imid1_extrap_normr_part4_train;

   infodata(k).mid1.part1_test = imid1_part1_test;
   infodata(k).mid1.part2_test = imid1_part2_test;
   infodata(k).mid1.part3_test = imid1_part3_test;
   infodata(k).mid1.part4_test = imid1_part4_test;

   infodata(k).mid1.part1_test_extrap = imid1_extrap_part1_test;
   infodata(k).mid1.part2_test_extrap = imid1_extrap_part2_test;
   infodata(k).mid1.part3_test_extrap = imid1_extrap_part3_test;
   infodata(k).mid1.part4_test_extrap = imid1_extrap_part4_test;

   infodata(k).mid1.part1_normr_test_extrap = imid1_extrap_normr_part1_test;
   infodata(k).mid1.part2_normr_test_extrap = imid1_extrap_normr_part2_test;
   infodata(k).mid1.part3_normr_test_extrap = imid1_extrap_normr_part3_test;
   infodata(k).mid1.part4_normr_test_extrap = imid1_extrap_normr_part4_test;



   % Second MID Parts 1-4:

   infodata(k).mid2.part1_train = imid2_part1_train;
   infodata(k).mid2.part2_train = imid2_part2_train;
   infodata(k).mid2.part3_train = imid2_part3_train;
   infodata(k).mid2.part4_train = imid2_part4_train;

   infodata(k).mid2.part1_train_extrap = imid2_extrap_part1_train;
   infodata(k).mid2.part2_train_extrap = imid2_extrap_part2_train;
   infodata(k).mid2.part3_train_extrap = imid2_extrap_part3_train;
   infodata(k).mid2.part4_train_extrap = imid2_extrap_part4_train;

   infodata(k).mid2.part1_normr_train_extrap = imid2_extrap_normr_part1_train;
   infodata(k).mid2.part2_normr_train_extrap = imid2_extrap_normr_part2_train;
   infodata(k).mid2.part3_normr_train_extrap = imid2_extrap_normr_part3_train;
   infodata(k).mid2.part4_normr_train_extrap = imid2_extrap_normr_part4_train;

   infodata(k).mid2.part1_test = imid2_part1_test;
   infodata(k).mid2.part2_test = imid2_part2_test;
   infodata(k).mid2.part3_test = imid2_part3_test;
   infodata(k).mid2.part4_test = imid2_part4_test;

   infodata(k).mid2.part1_test_extrap = imid2_extrap_part1_test;
   infodata(k).mid2.part2_test_extrap = imid2_extrap_part2_test;
   infodata(k).mid2.part3_test_extrap = imid2_extrap_part3_test;
   infodata(k).mid2.part4_test_extrap = imid2_extrap_part4_test;

   infodata(k).mid2.part1_normr_test_extrap = imid2_extrap_normr_part1_test;
   infodata(k).mid2.part2_normr_test_extrap = imid2_extrap_normr_part2_test;
   infodata(k).mid2.part3_normr_test_extrap = imid2_extrap_normr_part3_test;
   infodata(k).mid2.part4_normr_test_extrap = imid2_extrap_normr_part4_test;



   % First and Second MID Parts 1-4:

   infodata(k).mid12.part1_train = imid12_part1_train;
   infodata(k).mid12.part2_train = imid12_part2_train;
   infodata(k).mid12.part3_train = imid12_part3_train;
   infodata(k).mid12.part4_train = imid12_part4_train;

   infodata(k).mid12.part1_train_extrap = imid12_extrap_part1_train;
   infodata(k).mid12.part2_train_extrap = imid12_extrap_part2_train;
   infodata(k).mid12.part3_train_extrap = imid12_extrap_part3_train;
   infodata(k).mid12.part4_train_extrap = imid12_extrap_part4_train;

   infodata(k).mid12.part1_normr_train_extrap = imid12_extrap_normr_part1_train;
   infodata(k).mid12.part2_normr_train_extrap = imid12_extrap_normr_part2_train;
   infodata(k).mid12.part3_normr_train_extrap = imid12_extrap_normr_part3_train;
   infodata(k).mid12.part4_normr_train_extrap = imid12_extrap_normr_part4_train;

   infodata(k).mid12.part1_test = imid12_part1_test;
   infodata(k).mid12.part2_test = imid12_part2_test;
   infodata(k).mid12.part3_test = imid12_part3_test;
   infodata(k).mid12.part4_test = imid12_part4_test;

   infodata(k).mid12.part1_test_extrap = imid12_extrap_part1_test;
   infodata(k).mid12.part2_test_extrap = imid12_extrap_part2_test;
   infodata(k).mid12.part3_test_extrap = imid12_extrap_part3_test;
   infodata(k).mid12.part4_test_extrap = imid12_extrap_part4_test;

   infodata(k).mid12.part1_normr_test_extrap = imid12_extrap_normr_part1_test;
   infodata(k).mid12.part2_normr_test_extrap = imid12_extrap_normr_part2_test;
   infodata(k).mid12.part3_normr_test_extrap = imid12_extrap_normr_part3_test;
   infodata(k).mid12.part4_normr_test_extrap = imid12_extrap_normr_part4_test;

% pause(1)


end % (for k)




function plot_mid_proj_info(ista, imid1, imid2, imid12, datalength, numbins)

   % Calculate some mean values (take the mean across different binsizes):
   mnsta = mean(ista(2:end,:),1);
   mnmid1 = mean(imid1(2:end,:),1);
   mnmid2 = mean(imid2(2:end,:),1);
   mnmid12 = mean(imid12(2:end,:),1);

   [mnsta(:) mnmid1(:) mnmid2(:) mnmid12(:)]

pause

   % Rows of ista correspond to a particular binsize. The columns
   % of ista correspond to a datalength value, and correspond to the
   % entries in the vector datalength.

   figure;

   c = {'b' 'g' 'r' 'c' 'm' 'k' 'b' 'g' 'r' 'c' 'm' 'k' 'b' 'g' 'r' 'c' 'm' 'k'};
   c = c(1:length(numbins));
   d = {'-' '-' '-' '-' '-' '-' ':' ':' ':' ':' ':' ':' '-.' '-.' '-.' '-.' '-.' '-.'};

   for i = 1:length(numbins)
      binlabel{i} = sprintf('%.0f', numbins(i));
   end % (for i)

   datafrac = 50:10:100;

   for i = 1:length(datafrac)
      xticklabel{i} = sprintf('1/%.0f ', datafrac(i));
   end % (for i)



   % Plot the STA information data:
   subplot(3,2,1);
   for i = 1:size(ista,1)
      hold on;
      plot(1./datalength, ista(i,:), [c{i} d{i}]);
   end % (for i)
   plot(1./datalength, mnsta,  'ks-');
   title('I(sta)');
   xlim([0.0095 1/49]);
   set(gca,'tickdir', 'out');
   set(gca,'xtick', fliplr(1./datafrac));
   set(gca,'xticklabel', fliplr(xticklabel));
   legend(binlabel);



   % Plot the MID1 information data:
   subplot(3,2,3);
   for i = 1:size(imid1,1)
      hold on;
      plot(1./datalength, imid1(i,:), [c{i} d{i}]);
   end % (for i)
   plot(1./datalength, mnmid1,  'rs-');
   title('I(1)');
   xlim([0.0095 1/49]);
   set(gca,'tickdir', 'out');
   set(gca,'xtick', fliplr(1./datafrac));
   set(gca,'xticklabel', fliplr(xticklabel));


   % Plot the MID2 information data:
   subplot(3,2,2);
   for i = 1:size(imid2,1)
      hold on;
      plot(1./datalength, imid2(i,:), [c{i} d{i}]);
   end % (for i)
   plot(1./datalength, mnmid2,  'bs-');
   title('I(2)');
   xlim([0.0095 1/49]);
   set(gca,'tickdir', 'out');
   set(gca,'xtick', fliplr(1./datafrac));
   set(gca,'xticklabel', fliplr(xticklabel));


   % Plot the MID12 information data:
   subplot(3,2,4);
   for i = 1:size(imid12,1)
      hold on;
      plot(1./datalength, imid12(i,:), [c{i} d{i}]);
   end % (for i)
   plot(1./datalength, mnmid12,  'gs-');
   xlabel('1 / %Data');
   title('I(1,2)');
   xlim([0.0095 1/49]);
   set(gca,'tickdir', 'out');
   set(gca,'xtick', fliplr(1./datafrac));
   set(gca,'xticklabel', fliplr(xticklabel));


   % Plot the mean information values:
   subplot(3,2,5);
   hold on;
   plot(1./datalength, mnsta,  'ks-');
   plot(1./datalength, mnmid1,  'rs-');
   plot(1./datalength, mnmid2,  'bs-');
   plot(1./datalength, mnmid12,  'gs-');
   plot(1./datalength, mnmid1+mnmid2,'cs-');
   xlabel('1/%Data');
   legend('STA', 'MID1', 'MID2', 'MID12','MID1+MID2');
   xlim([0.0095 1/49]);
   set(gca,'tickdir', 'out');
   set(gca,'xtick', fliplr(1./datafrac));
   set(gca,'xticklabel', fliplr(xticklabel));

%    subplot(3,2,6);
%    hold on;
%    plot(datalength, mnsta,  'ks-');
%    plot(datalength, mnmid1,  'rs-');
%    plot(datalength, mnmid2,  'bs-');
%    plot(datalength, mnmid12,  'gs-');
%    xlabel('%Data');

return;





function [istainf, imid1inf, imid2inf, imid12inf, normr_sta, normr_mid1, normr_mid2, normr_mid12] = ...
      get_info_extrapolation(cell_num, exp, site, chan, model, ista, imid1, imid2, imid12, datalength, numbins)


% cell_num
% exp
% site
% chan
% model
% pause

istainf = zeros( size(ista,1), 1);
imid1inf = zeros( size(imid1,1), 1);
imid2inf = zeros( size(imid2,1), 1);
imid12inf = zeros( size(imid12,1), 1);

normr_sta = zeros( size(ista,1), 1);
normr_mid1 = zeros( size(imid1,1), 1);
normr_mid2 = zeros( size(imid2,1), 1);
normr_mid12 = zeros( size(imid12,1), 1);

invdatafrac = fliplr( 1 ./ datalength );

for i = 1:length(numbins)

   nbins = numbins(i);

   infosta = ista(i,:);
   infomid1 = imid1(i,:);
   infomid2 = imid2(i,:);
   infomid12 = imid12(i,:);

   invsta   = fliplr(infosta);
   invmid1  = fliplr(infomid1);
   invmid2  = fliplr(infomid2);
   invmid12 = fliplr(infomid12);

   xfit = linspace(0,max(invdatafrac));

   index = 1:5;%length(invdatafrac);

   [betasta,s]  = polyfit(invdatafrac(index), invsta(index),1);
   ysta = polyval(betasta,xfit); % for plotting
%    extrap_fit_sta = polyval(betasta, invdatafrac); % for r^2 calc
   nr_sta = s.normr;

% [length(ysta) length(invsta) length(extrap_fit_sta) normres_sta]

   [betamid1,s] = polyfit(invdatafrac(index), invmid1(index),1);
   ymid1 = polyval(betamid1,xfit);
   nr_mid1 = s.normr;


   [betamid2,s] = polyfit(invdatafrac(index), invmid2(index),1);
   ymid2 = polyval(betamid2,xfit);
   nr_mid2 = s.normr;


   [betamid12,s]= polyfit(invdatafrac(index), invmid12(index),1);
   ymid12 = polyval(betamid12,xfit);
   nr_mid12 = s.normr;


   istainf(i) = betasta(2);
   normr_sta(i) = nr_sta;

   imid1inf(i) = betamid1(2);
   normr_mid1(i) = nr_mid1;

   imid2inf(i) = betamid2(2);
   normr_mid2(i) = nr_mid2;

   imid12inf(i) = betamid12(2);
   normr_mid12(i) = nr_mid12;


   minmin = min([infosta infomid1 infomid2 infomid12]);
   maxmax = max([infosta infomid1 infomid2 infomid12]);
   range = maxmax - minmin;

   % Plot the mean information values:
   clf;
   hold on;
   plot(invdatafrac, invsta,  'ks');
   plot(invdatafrac, invmid1,  'rs');
   plot(invdatafrac, invmid2,  'bs');
   plot(invdatafrac, invmid12,  'gs');

   plot(xfit, ysta, 'k-');
   plot(xfit, ymid1, 'r-');
   plot(xfit, ymid2, 'b-');
   plot(xfit, ymid12, 'g-');

   datafrac = 50:10:100;

   for j = 1:length(datafrac)
      xticklabel{j} = sprintf('1/%.0f ', datafrac(j));
   end % (for i)

   xlabel('1/%Data');
   title(sprintf('#%.0f: %s - site%.0f - chan %.0f - model %.0f, Numbins = %.0f', cell_num, exp, site, chan, model, numbins(i)));
   legend('STA', 'MID1', 'MID2', 'MID12');
   xlim([0.0095 1/49]);
   ylim([0 maxmax+0.05*range]);
%    ylim([minmin-0.05*range maxmax+0.05*range]);
   set(gca,'tickdir', 'out');
   set(gca,'xtick', fliplr(1./datafrac));
   set(gca,'xticklabel', fliplr(xticklabel));

pause(0.5);

end

return;






















%    % Now we have raw projection values. Next we want to normalize the
%    % values with respect to the standard deviation of p(proj|no spike)
%    mn = mean(projrand);
%    sd = std(projrand);
%    projspk_scaled = (projspk - mn) ./ sd;
%    projrand_scaled = (projrand - mn) ./ sd;
% 
% 
%    % Finally, we calculate probability distributions for different
%    % bin sizes
%    [binsize, info] = get_1d_filter_info(projspk_scaled, projrand_scaled)
% 
% 
%    subplot(3,2,1);
%    hist(projspk,25);
% 
%    subplot(3,2,3);
%    hist(projrand,25);
% 
%    subplot(3,2,5);
%    [nspk, xspk] = hist(projspk_scaled, 25);
%    [nr, xr] = hist(projrand_scaled, 25);
%    plot(xspk, nspk./sum(nspk), 'r-', xr, nr./sum(nr), 'k-');
%    legend('spk', 'no spk');
% 
%    subplot(3,2,2);
%    plot(binsize, info, 'ko-');
%    title('versus bin size');
% 
%    subplot(3,2,4);
%    plot(1 ./ binsize, info, 'ko-');
%    title('versus 1 / bin size');
% 
% 
%    % Now we have raw projection values. Next we want to normalize the
%    % values with respect to the standard deviation of p(proj|no spike)
%    mn = mean(projrand);
%    sd = std(projrand);
%    projspk_scaled = (projspk - mn) ./ sd;
%    projrand_scaled = (projrand - mn) ./ sd;
% 
% 
%    % Finally, we calculate probability distributions for different
%    % bin sizes
%    [binsize, info] = get_2d_filter_info(projspk_scaled, projrand_scaled)



