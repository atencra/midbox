function info = calculate_info(proj)
% calculate_info - compute mutual information for filters 
%    from projection values 
%
% info = calculate_info(proj)
% ----------------------------------------------------------
%
% calculate_info computes mutual information values for
% filters as a function of dataset size and bin size. The
% dataset sizes that are used are 50, 60, 70, 80, 90, 92.5,
% 95, 97.5, and 100 % of the data. Binsize is a function of
% the number of bins used to make a histogram for the 
% probability distribution of projection values. The number
% of bins is 11 to 25 in steps of 2.
%
% proj is a struct array holding all the projection data. It 
% is created by using the call:
%
% [proj] = calculate_projection(mid, spk, trigger, specfile);
%
% caa 12/18/06


% proj = struct array with the following fields:
% 
%     exp
%     site
%     chan
%     model
%     depth
%     position
%     stim
%     atten
%     spl
%     sm
%     tm
%     mdb
%     trigger
%     spetrand
%     randspklocation
%     spklocation
%     xstarand
%     xstaspk

%     x1x2rand_mean
%     x1x2spk_mean
%     x1xrandspk_mean
%     x1x2rand_part1
%     x1x2spk_part1
%     x1xrandspk_part1
%     x1x2rand_part2
%     x1x2spk_part2
%     x1xrandspk_part2
%     x1x2rand_part3
%     x1x2spk_part3
%     x1xrandspk_part3
%     x1x2rand_part4
%     x1x2spk_part4
%     x1xrandspk_part4



datalength = [50 60 70 80 90 92.5 95 97.5 100];
binboundary = 6;
numbins = (11:2:25)+1;
numreps = 1; % was used earlier for st error measure; not used now

for i = 1:length(proj)

   trigger = proj(i).trigger;
   randspklocation = proj(i).randspklocation;
   spklocation = proj(i).spklocation;

   % Assign STA data:
   % --------------------------------------------
   xstarand_mean = proj(i).xstarand_mean; % Nx1
   xstaspk_mean = proj(i).xstaspk_mean; % Nx1

   xstarand_part1 = proj(i).xstarand_part1;
   xstaspk_part1 = proj(i).xstaspk_part1; 

   xstarand_part2 = proj(i).xstarand_part2; 
   xstaspk_part2 = proj(i).xstaspk_part2; 

   xstarand_part3 = proj(i).xstarand_part3; 
   xstaspk_part3 = proj(i).xstaspk_part3; 

   xstarand_part4 = proj(i).xstarand_part4; 
   xstaspk_part4 = proj(i).xstaspk_part4; 


   % Assign MID data:
   % --------------------------------------------
   x1x2rand_mean = proj(i).x1x2rand_mean; % Nx2
   x1x2spk_mean = proj(i).x1x2spk_mean; % Nx2
   x1xrandspk_mean = proj(i).x1xrandspk_mean; % Nx2

   x1x2rand_part1 = proj(i).x1x2rand_part1;
   x1x2spk_part1 = proj(i).x1x2spk_part1;
   x1xrandspk_part1 = proj(i).x1xrandspk_part1;

   x1x2rand_part2 = proj(i).x1x2rand_part2;
   x1x2spk_part2 = proj(i).x1x2spk_part2;
   x1xrandspk_part2 = proj(i).x1xrandspk_part2;

   x1x2rand_part3 = proj(i).x1x2rand_part3;
   x1x2spk_part3 = proj(i).x1x2spk_part3;
   x1xrandspk_part3 = proj(i).x1xrandspk_part3;

   x1x2rand_part4 = proj(i).x1x2rand_part4;
   x1x2spk_part4 = proj(i).x1x2spk_part4;
   x1xrandspk_part4 = proj(i).x1xrandspk_part4;


   % First make the calculations over the complete data set for each filter:

   for j = 1:length(datalength)

      % Get the appropriate fraction of the spike projection data.
      triggerfrac = floor( datalength(j)/100 * length(trigger) );
      indspkfrac = find( spklocation <= triggerfrac );

      xstaspk_mean_frac  = xstaspk_mean(indspkfrac);
      xstaspk_part1_frac = xstaspk_part1(indspkfrac);
      xstaspk_part2_frac = xstaspk_part2(indspkfrac);
      xstaspk_part3_frac = xstaspk_part3(indspkfrac);
      xstaspk_part4_frac = xstaspk_part4(indspkfrac);

      x1x2spk_mean_frac    = x1x2spk_mean(indspkfrac,:);
      x1xrandspk_mean_frac = x1xrandspk_mean(indspkfrac,:);

      x1x2spk_part1_frac    = x1x2spk_part1(indspkfrac,:);
      x1xrandspk_part1_frac = x1xrandspk_part1(indspkfrac,:);

      x1x2spk_part2_frac    = x1x2spk_part2(indspkfrac,:);
      x1xrandspk_part2_frac = x1xrandspk_part2(indspkfrac,:);

      x1x2spk_part3_frac    = x1x2spk_part3(indspkfrac,:);
      x1xrandspk_part3_frac = x1xrandspk_part3(indspkfrac,:);

      x1x2spk_part4_frac    = x1x2spk_part4(indspkfrac,:);
      x1xrandspk_part4_frac = x1xrandspk_part4(indspkfrac,:);


      % Now get the data fractions for the cross-validation, or test,
      % information calculations:

      trigger_quarter = floor( 0.25 * length(trigger) ); % test data is only 1/4 of the total data length
      trigger_quarter_frac = floor( datalength(j)/100 * trigger_quarter ); % take a fraction of it

      index_part1 = find(spklocation <= trigger_quarter_frac);
      index_part2 = find(0.25*length(trigger) <= spklocation & spklocation <= trigger_quarter_frac+0.25*length(trigger) );
      index_part3 = find(0.5*length(trigger) <= spklocation & spklocation <= trigger_quarter_frac+0.5*length(trigger) );
      index_part4 = find(0.75*length(trigger) <= spklocation & spklocation <= trigger_quarter_frac+0.75*length(trigger) );

      xstaspk_part1_frac_test = xstaspk_part1(index_part1);
      xstaspk_part2_frac_test = xstaspk_part2(index_part2);
      xstaspk_part3_frac_test = xstaspk_part3(index_part3);
      xstaspk_part4_frac_test = xstaspk_part4(index_part4);

      x1x2spk_part1_frac_test = x1x2spk_part1(index_part1,:);
      x1x2spk_part2_frac_test = x1x2spk_part2(index_part2,:);
      x1x2spk_part3_frac_test = x1x2spk_part3(index_part3,:);
      x1x2spk_part4_frac_test = x1x2spk_part4(index_part4,:);


      % We also want to calculate the mutual information values that were
      % used to compute the filters. For part 1 this means using data from
      % the 2nd, 3rd, and 4th quarters of the dataset. For part 2 it is the
      % 1st, 3rd, and 4th, quarters. For part 3 is the 1st, 2nd, and 4th
      % quarters. For part 4 it is the 1st, 2nd, and 3rd quarters.

      trigger_index_1st = 1:floor(0.25*length(trigger));
      trigger_index_2nd = (floor(0.25*length(trigger))+1):floor(0.5*length(trigger));
      trigger_index_3rd = (floor(0.5*length(trigger))+1):floor(0.75*length(trigger));
      trigger_index_4th = (floor(0.75*length(trigger))+1):length(trigger);

      trigger_part_length = floor( 0.75 * length(trigger) ); % datalength for each part's calculations
      trigger_part_length_frac = floor( datalength(j)/100 * trigger_part_length );

      trigger_part1 = [trigger_index_2nd trigger_index_3rd trigger_index_4th];
      trigger_part1_frac = trigger_part1(1:trigger_part_length_frac);
      index_part1 = ismember(spklocation, trigger_part1_frac);

      xstaspk_part1_frac_train = xstaspk_part1(index_part1);
      x1x2spk_part1_frac_train = x1x2spk_part1(index_part1,:);

      trigger_part2 = [trigger_index_1st trigger_index_3rd trigger_index_4th];
      trigger_part2_frac = trigger_part2(1:trigger_part_length_frac);
      index_part2 = ismember(spklocation, trigger_part2_frac);

      xstaspk_part2_frac_train = xstaspk_part2(index_part2);
      x1x2spk_part2_frac_train = x1x2spk_part2(index_part2,:);

      trigger_part3 = [trigger_index_1st trigger_index_2nd trigger_index_4th];
      trigger_part3_frac = trigger_part3(1:trigger_part_length_frac);
      index_part3 = ismember(spklocation, trigger_part3_frac);

      xstaspk_part3_frac_train = xstaspk_part3(index_part3);
      x1x2spk_part3_frac_train = x1x2spk_part3(index_part3,:);

      trigger_part4 = [trigger_index_1st trigger_index_2nd trigger_index_3rd];
      trigger_part4_frac = trigger_part4(1:trigger_part_length_frac);
      index_part4 = ismember(spklocation, trigger_part4_frac);

      xstaspk_part4_frac_train = xstaspk_part4(index_part4);
      x1x2spk_part4_frac_train = x1x2spk_part4(index_part4,:);

      % We won't take fractions of the data to get the prior distribution. 
      % This will allow us to get a clean estimate of the prior.

      for k = 1:length(numbins)

         fprintf('Data==%.0f\tNbins==%.0f\n',datalength(j),numbins(k)-1);

         infosta_mean_temp = [];

         infosta_part1_temp = [];
         infosta_part1_test_temp = [];
         infosta_part1_train_temp = [];

         infosta_part2_temp = [];
         infosta_part2_test_temp = [];
         infosta_part2_train_temp = [];

         infosta_part3_temp = [];
         infosta_part3_test_temp = [];
         infosta_part3_train_temp = [];

         infosta_part4_temp = [];
         infosta_part4_test_temp = [];
         infosta_part4_train_temp = [];


         infox1x2_mean_temp = [];
         infox1xrand_mean_temp = [];

         infox1x2_part1_temp = [];
         infox1x2_part1_test_temp = [];
         infox1x2_part1_train_temp = [];
         infox1xrand_part1_temp = [];

         infox1x2_part2_temp = [];
         infox1x2_part2_test_temp = [];
         infox1x2_part2_train_temp = [];
         infox1xrand_part2_temp = [];

         infox1x2_part3_temp = [];
         infox1x2_part3_test_temp = [];
         infox1x2_part3_train_temp = [];
         infox1xrand_part3_temp = [];

         infox1x2_part4_temp = [];
         infox1x2_part4_test_temp = [];
         infox1x2_part4_train_temp = [];
         infox1xrand_part4_temp = [];


         nbins = numbins(k);

         % For the sta:
         % ------------------------------------------------
         [ista] = get_1d_filter_info(xstaspk_mean_frac, xstarand_mean, nbins);
         infosta_mean_temp = [infosta_mean_temp; ista];
         clear('ista');


         % STA: Part 1
         % ------------------------------------------------
         [ista] = get_1d_filter_info(xstaspk_part1_frac, xstarand_part1, nbins);
         infosta_part1_temp = [infosta_part1_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part1_frac_test, xstarand_part1, nbins);
         infosta_part1_test_temp = [infosta_part1_test_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part1_frac_train, xstarand_part1, nbins);
         infosta_part1_train_temp = [infosta_part1_train_temp; ista];
         clear('ista');


         % STA: Part 2
         % ------------------------------------------------
         [ista] = get_1d_filter_info(xstaspk_part2_frac, xstarand_part2, nbins);
         infosta_part2_temp = [infosta_part2_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part2_frac_test, xstarand_part2, nbins);
         infosta_part2_test_temp = [infosta_part2_test_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part2_frac_train, xstarand_part2, nbins);
         infosta_part2_train_temp = [infosta_part2_train_temp; ista];
         clear('ista');


         % STA: Part 3
         % ------------------------------------------------
         [ista] = get_1d_filter_info(xstaspk_part3_frac, xstarand_part3, nbins);
         infosta_part3_temp = [infosta_part3_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part3_frac_test, xstarand_part3, nbins);
         infosta_part3_test_temp = [infosta_part3_test_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part3_frac_train, xstarand_part3, nbins);
         infosta_part3_train_temp = [infosta_part3_train_temp; ista];
         clear('ista');


         % STA: Part 4
         % ------------------------------------------------
         [ista] = get_1d_filter_info(xstaspk_part4_frac, xstarand_part4, nbins);
         infosta_part4_temp = [infosta_part4_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part4_frac_test, xstarand_part4, nbins);
         infosta_part4_test_temp = [infosta_part4_test_temp; ista];
         clear('ista');

         [ista] = get_1d_filter_info(xstaspk_part4_frac_train, xstarand_part4, nbins);
         infosta_part4_train_temp = [infosta_part4_train_temp; ista];
         clear('ista');



         % For the mean mids:
         % ------------------------------------------------
         [info12, info1, info2] = get_2d_filter_info(x1x2spk_mean_frac, x1x2rand_mean, nbins);
         infox1x2_mean_temp = [infox1x2_mean_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1xrandspk_mean_frac, x1x2rand_mean, nbins);
         infox1xrand_mean_temp = [infox1xrand_mean_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');


         % For part 1:
         % ------------------------------------------------
         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part1_frac, x1x2rand_part1, nbins);
         infox1x2_part1_temp = [infox1x2_part1_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part1_frac_test, x1x2rand_part1, nbins);
         infox1x2_part1_test_temp = [infox1x2_part1_test_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part1_frac_train, x1x2rand_part1, nbins);
         infox1x2_part1_train_temp = [infox1x2_part1_train_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1xrandspk_part1_frac, x1x2rand_part1, nbins);
         infox1xrand_part1_temp = [infox1xrand_part1_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');


         % For part 2:
         % ------------------------------------------------
         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part2_frac, x1x2rand_part2, nbins);
         infox1x2_part2_temp = [infox1x2_part2_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part2_frac_test, x1x2rand_part2, nbins);
         infox1x2_part2_test_temp = [infox1x2_part2_test_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part2_frac_train, x1x2rand_part2, nbins);
         infox1x2_part2_train_temp = [infox1x2_part2_train_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1xrandspk_part2_frac, x1x2rand_part2, nbins);
         infox1xrand_part2_temp = [infox1xrand_part2_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');


         % For part 3:
         % ------------------------------------------------
         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part3_frac, x1x2rand_part3, nbins);
         infox1x2_part3_temp = [infox1x2_part3_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part3_frac_test, x1x2rand_part3, nbins);
         infox1x2_part3_test_temp = [infox1x2_part3_test_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part3_frac_train, x1x2rand_part3, nbins);
         infox1x2_part3_train_temp = [infox1x2_part3_train_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1xrandspk_part3_frac, x1x2rand_part3, nbins);
         infox1xrand_part3_temp = [infox1xrand_part3_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');


         % For part 4:
         % ------------------------------------------------
         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part4_frac, x1x2rand_part4, nbins);
         infox1x2_part4_temp = [infox1x2_part4_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part4_frac_test, x1x2rand_part4, nbins);
         infox1x2_part4_test_temp = [infox1x2_part4_test_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1x2spk_part4_frac_train, x1x2rand_part4, nbins);
         infox1x2_part4_train_temp = [infox1x2_part4_train_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');

         [info12, info1, info2] = get_2d_filter_info(x1xrandspk_part4_frac, x1x2rand_part4, nbins);
         infox1xrand_part4_temp = [infox1xrand_part4_temp; info12 info1 info2];
         clear('info12', 'info1', 'info2');


         % Assign the data to cell array data structures:
         % ------------------------------------------------
         infosta_mean{j}{k}        = infosta_mean_temp; % save the data into a cell array

         infosta_part1{j}{k}       = infosta_part1_temp;
         infosta_part1_test{j}{k}  = infosta_part1_test_temp;
         infosta_part1_train{j}{k} = infosta_part1_train_temp;

         infosta_part2{j}{k}       = infosta_part2_temp;
         infosta_part2_test{j}{k}  = infosta_part2_test_temp;
         infosta_part2_train{j}{k} = infosta_part2_train_temp;

         infosta_part3{j}{k}       = infosta_part3_temp;
         infosta_part3_test{j}{k}  = infosta_part3_test_temp;
         infosta_part3_train{j}{k} = infosta_part3_train_temp;

         infosta_part4{j}{k}       = infosta_part4_temp;
         infosta_part4_test{j}{k}  = infosta_part4_test_temp;
         infosta_part4_train{j}{k} = infosta_part4_train_temp;


         infomid1mid2_mean{j}{k} = infox1x2_mean_temp; % save the data into a cell array
         infomid1mid2rand_mean{j}{k} = infox1xrand_mean_temp; % save the data into a cell array

         infomid1mid2_part1{j}{k} = infox1x2_part1_temp;
         infomid1mid2_part1_test{j}{k} = infox1x2_part1_test_temp;
         infomid1mid2_part1_train{j}{k} = infox1x2_part1_train_temp;
         infomid1mid2rand_part1{j}{k} = infox1xrand_part1_temp;

         infomid1mid2_part2{j}{k} = infox1x2_part2_temp;
         infomid1mid2_part2_test{j}{k} = infox1x2_part2_test_temp;
         infomid1mid2_part2_train{j}{k} = infox1x2_part2_train_temp;
         infomid1mid2rand_part2{j}{k} = infox1xrand_part2_temp;

         infomid1mid2_part3{j}{k} = infox1x2_part3_temp;
         infomid1mid2_part3_test{j}{k} = infox1x2_part3_test_temp;
         infomid1mid2_part3_train{j}{k} = infox1x2_part3_train_temp;
         infomid1mid2rand_part3{j}{k} = infox1xrand_part3_temp;

         infomid1mid2_part4{j}{k} = infox1x2_part4_temp;
         infomid1mid2_part4_test{j}{k} = infox1x2_part4_test_temp;
         infomid1mid2_part4_train{j}{k} = infox1x2_part4_train_temp;
         infomid1mid2rand_part4{j}{k} = infox1xrand_part4_temp;

      end % (for k)

   end % (for j)


   % Save the data in the struct array info
   % ------------------------------------------------
   info(i).exp = proj(i).exp;
   info(i).site = proj(i).site;
   info(i).chan = proj(i).chan;
   info(i).model = proj(i).model;
   info(i).depth = proj(i).depth;
   info(i).position = proj(i).position;
   info(i).stim = proj(i).stim;
   info(i).atten = proj(i).atten;
   info(i).spl = proj(i).spl;
   info(i).sm = proj(i).sm;
   info(i).tm = proj(i).tm;
   info(i).mdb = proj(i).mdb;
   info(i).datalength = datalength;
   info(i).numbins = numbins-1;
   info(i).numreps = numreps;

   % Assign the STA information values:
   info(i).sta_mean        = infosta_mean;

   info(i).sta_part1       = infosta_part1;
   info(i).sta_part1_test  = infosta_part1_test;
   info(i).sta_part1_train = infosta_part1_train;

   info(i).sta_part2       = infosta_part2;
   info(i).sta_part2_test  = infosta_part2_test;
   info(i).sta_part2_train = infosta_part2_train;

   info(i).sta_part3       = infosta_part3;
   info(i).sta_part3_test  = infosta_part3_test;
   info(i).sta_part3_train = infosta_part3_train;

   info(i).sta_part4       = infosta_part4;
   info(i).sta_part4_test  = infosta_part4_test;
   info(i).sta_part4_train = infosta_part4_train;


   % Assign the MID information values:
   info(i).mid1mid2_mean        = infomid1mid2_mean;
   info(i).mid1mid2rand_mean    = infomid1mid2rand_mean;

   info(i).mid1mid2_part1       = infomid1mid2_part1;
   info(i).mid1mid2_part1_test  = infomid1mid2_part1_test;
   info(i).mid1mid2_part1_train = infomid1mid2_part1_train;
   info(i).mid1mid2rand_part1   = infomid1mid2rand_part1;

   info(i).mid1mid2_part2       = infomid1mid2_part2;
   info(i).mid1mid2_part2_test  = infomid1mid2_part2_test;
   info(i).mid1mid2_part2_train = infomid1mid2_part2_train;
   info(i).mid1mid2rand_part2   = infomid1mid2rand_part2;

   info(i).mid1mid2_part3       = infomid1mid2_part3;
   info(i).mid1mid2_part3_test  = infomid1mid2_part3_test;
   info(i).mid1mid2_part3_train = infomid1mid2_part3_train;
   info(i).mid1mid2rand_part3   = infomid1mid2rand_part3;

   info(i).mid1mid2_part4       = infomid1mid2_part4;
   info(i).mid1mid2_part4_test  = infomid1mid2_part4_test;
   info(i).mid1mid2_part4_train = infomid1mid2_part4_train;
   info(i).mid1mid2rand_part4   = infomid1mid2rand_part4;


   % Clear out the temporary data variables
   clear('infosta_mean');

   clear('infosta_part1');
   clear('infosta_part1_test');
   clear('infosta_part1_train');

   clear('infosta_part2');
   clear('infosta_part2_test');
   clear('infosta_part2_train');

   clear('infosta_part3');
   clear('infosta_part3_test');
   clear('infosta_part3_train');

   clear('infosta_part4');
   clear('infosta_part4_test');
   clear('infosta_part4_train');

   clear('infomid1mid2_mean');
   clear('infomid1mid2rand_mean');

   clear('infomid1mid2_part1');
   clear('infomid1mid2_part1_test');
   clear('infomid1mid2_part1_train');
   clear('infomid1mid2rand_part1');

   clear('infomid1mid2_part2');
   clear('infomid1mid2_part2_test');
   clear('infomid1mid2_part2_train');
   clear('infomid1mid2rand_part2');

   clear('infomid1mid2_part3');
   clear('infomid1mid2_part3_test');
   clear('infomid1mid2_part3_train');
   clear('infomid1mid2rand_part3');

   clear('infomid1mid2_part4');
   clear('infomid1mid2_part4_test');
   clear('infomid1mid2_part4_train');
   clear('infomid1mid2rand_part4');

   % Notate the progress
   fprintf('%s - site%.0f - chan%.0f - model%.0f  -  %.0f of %.0f completed.\n', ...
      proj(i).exp, proj(i).site, proj(i).chan, proj(i).model, i, length(proj));
   pause(0.5);

end % (for i)



function [info] = get_1d_filter_info(projspk, projrand, numbins)

   %bins = linspace(-binboundary, binboundary, numbins); % equally spaced bins

   % We have raw projection values. Next we want to normalize the
   % values with respect to the standard deviation of p(proj|no spike)
   mn = mean(projrand);
   sd = std(projrand);
   projspk_scaled = (projspk - mn) ./ sd;
   projrand_scaled = (projrand - mn) ./ sd;

%    minbin = min([min(projspk_scaled) min(projrand_scaled)]);
%    maxbin = max([max(projspk_scaled) max(projrand_scaled)]);
%    maxmax = max([abs(minbin) maxbin]);
%    bins = linspace(-maxmax, maxmax, numbins); % equally spaced bins


   minbin = min( projspk_scaled );
   maxbin = max( projspk_scaled );
   binedges = linspace(minbin, maxbin, numbins); % equally spaced bins


   % Calculate probability distributions
   [nspk] = histc(projspk_scaled, binedges);
   nspk = nspk(1:end-1);
   pxspk = nspk ./ sum(nspk);

   [nrand] = histc(projrand_scaled, binedges);
   nrand = nrand(1:end-1);
   px = nrand ./ sum(nrand);


   if ( sum(nspk) )
      index = pxspk>0 & px>0;
      info = sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) );
   else
      info = 0;
   end

return;



function [info12, info1, info2, x1bins, x2bins] = get_2d_filter_info(projspk, projrand, numbins)

%    bincenters = linspace(-binboundary, binboundary, numbins); % equally spaced bins
%    binedges = center2edge(bincenters);

   % We have raw projection values. Next we want to normalize the
   % all values by subtracting out the mean and dividing by the 
   % standard deviation of the prior distribution, p(proj|no spike), i.e.
   % projrand


% x0 = min(min(x1_mtx));
% x1 = max(max(x1_mtx));
% maxmax = max([abs(x0) x1]);
% xedges = linspace(-maxmax, maxmax, Nbins);
% 
% y0 = min(min(x2_mtx));
% y1 = max(max(x2_mtx));
% maxmax = max([abs(y0) y1]);
% yedges = linspace(-maxmax, maxmax, Nbins);


   mn1 = mean(projrand(:,1));
   mn2 = mean(projrand(:,2));

   sd1 = std(projrand(:,1));
   sd2 = std(projrand(:,2));

   projspk_scaled(:,1) = ( projspk(:,1) - mn1 ) ./ sd1;
   projspk_scaled(:,2) = ( projspk(:,2) - mn2 ) ./ sd2;

   projrand_scaled(:,1) = ( projrand(:,1) - mn1 ) ./ sd1;
   projrand_scaled(:,2) = ( projrand(:,2) - mn2 ) ./ sd2;


%    x1minbin = min([min(projspk_scaled(:,1)) min(projrand_scaled(:,1))]);
%    x1maxbin = max([max(projspk_scaled(:,1)) max(projrand_scaled(:,1))]);
%    x1maxmax = max([abs(x1minbin) x1maxbin]);
%   x1bins = linspace(-x1maxmax, x1maxmax, numbins); % equally spaced bin centers
%   x1binedges = center2edge(x1bins); % make bin edges for hist2d function

   x1min = min( projspk_scaled(:,1) );
   x1max = max( projspk_scaled(:,1) );
   x1bins = linspace(x1min, x1max, numbins); % equally spaced bin edges
   x1binedges = x1bins; %center2edge(x1bins); % make bin edges for hist2d function


%    x2minbin = min([min(projspk_scaled(:,2)) min(projrand_scaled(:,2))]);
%    x2maxbin = max([max(projspk_scaled(:,2)) max(projrand_scaled(:,2))]);
%    x2maxmax = max([abs(x2minbin) x2maxbin]);
%    x2bins = linspace(-x2maxmax, x2maxmax, numbins); % equally spaced bin centers
%    x2binedges = center2edge(x2bins); % make bin edges for hist2d function

   x2min = min( projspk_scaled(:,2) );
   x2max = max( projspk_scaled(:,2) );
   x2bins = linspace(x2min, x2max, numbins); % equally spaced bin edges
   x2binedges = x2bins; %center2edge(x2bins); % make bin edges for hist2d function


   
   % projspk_scaled and projrand_scaled contain data for mid1 and mid2.
   % The mid1 values are in the first column, and the mid2 values in 
   % the second column.

   % Now we need to bin the projection values and obtain the probability
   % distributions:

   nspk = hist2d(projspk_scaled, x1binedges, x2binedges);
   px1x2spk = nspk ./ sum(sum(nspk)); % normalize


   nrand = hist2d(projrand_scaled, x1binedges, x2binedges);
   px1x2 = nrand ./ sum(sum(nrand)); % normalize


   if ( sum(nspk) )

      % Now calculate the 2D information using column vectors:
      pxspk = px1x2spk(:);
      px = px1x2(:);
      index = pxspk>0 & px>0; % log2(0) is tricky! so don't do it!
      info12 = sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) );

      % Now calculate the mid1 1D information:
      pxspk = sum(px1x2spk,2); % sum across columns
      px = sum(px1x2,2);
      index = pxspk>0 & px>0;
      info1 = sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) );

      % Calculate the mid2 1D information:
      pxspk = sum(px1x2spk,1); % sum across rows
      px = sum(px1x2,1);

      index = pxspk>0 & px>0;
      info2 = sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) );

   else

      info12 = 0;
      info1 = 0;
      info2 = 0;

   end

return;










