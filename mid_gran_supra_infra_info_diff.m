function mid_gran_supra_infra_info_diff(infocell)
%mid_gran_supra_infra_info_diff - Do MID information values vary with
%depth?
%
% For each penetration, the difference between information values in gran,
% supragranular, and infragranular layers is plotted. This means that 
% granular-supragranular, granular-infragranular, and supragranular-infragranular
% comparisons are made. From these comparisons, the statistical difference
% is computed using a signed-rank test.
%
% infocell : cell array, with each element containing the information 
%				values for a single penetration.
%
% exportfig(gcf,'fig678_layer4_layer23_layer5_layer6_btmf.eps', ...
%'color', 'bw', 'width', 6.3, 'height', 3.94, 'fontmode', 'fixed', 'fontsize', 8);
%
%   caa 12/1/08


% We need to process the following pairs:
% 
% gran_supra
% gran_infra
% supra_infra
% 
% That's all.


% var's to hold mtf best frequency data

contrib_gran_supra_pairs = [];
contrib_gran_infra_pairs = [];
contrib_supra_infra_pairs = [];

synergy_gran_supra_pairs = [];
synergy_gran_infra_pairs = [];
synergy_supra_infra_pairs = [];


% Define the layers as:

gran = [700 1100];
supra = [0 600];
infra = [1200 1900];


for i = 1:length(infocell)

   info_params = infocell{i};

	position = [];
	info_sta = [];
	info1 = [];
	info2 = [];
	info_both = [];
	contrib = [];
	synergy = [];

	for j = 1:length(info_params)

	   position(j) = info_params(j).position;

		info_sta(j) = mean( info_params(j).sta.information );
		info1(j) = mean( info_params(j).mid1.information );
		info2(j) = mean( info_params(j).mid2.information );
		info_both(j) = mean( info_params(j).mid12.information );

	end % (for j)


	% Error check to get rid of bad INF values:

	index = ~isinf(info_sta) & ~isinf(info1) & ...
			  ~isinf(info2) & ~isinf(info_both);

	info_sta = info_sta(index);
	info1 = info1(index);
	info2 = info2(index);
	info_both = info_both(index);

	position = position(index);

	contrib = 100 * info1 ./ info_both;
	synergy = 100 * info_both ./ (info1 + info2);



	% Parameter we really care about are asi and sepindex


   % ==============================================================
   %       all possible pairs of bmf, width, and cutoff freq
   % ==============================================================

   % process pairs gran-supra, gran-infra

   indgran = find( position >= gran(1) & position <= gran(2) );

   if ( ~isempty(indgran) )

      contrib_gran = contrib(indgran);
      synergy_gran = synergy(indgran);

      % ----- Supragranular Layer Data -----
      indsupra = find( position <= supra(2) );

      if ( ~isempty(indsupra) )

	      contrib_supra = contrib(indsupra);
			contrib_supra = contrib_supra(:);

		   synergy_supra = synergy(indsupra);
			synergy_supra = synergy_supra(:);


         for k = 1:length(contrib_gran)

            contrib_gran_supra_pairs = [contrib_gran_supra_pairs; ...
						contrib_gran(k)*ones(size(contrib_supra)) contrib_supra];

            synergy_gran_supra_pairs = [synergy_gran_supra_pairs; ...
						synergy_gran(k)*ones(size(synergy_supra)) synergy_supra];

         end % (for k = 1:length(contrib_gran) )

      end % (if)


      % ----- Infragranular Layer Data -----
      indinfra = find( position >= infra(1) & position <= infra(2) );

      if ( ~isempty(indsupra) )

	      contrib_infra = contrib(indinfra);
			contrib_infra = contrib_infra(:);

		   synergy_infra = synergy(indinfra);
			synergy_infra = synergy_infra(:);

         for k = 1:length(contrib_gran)

            contrib_gran_infra_pairs = [contrib_gran_infra_pairs; ...
						contrib_gran(k)*ones(size(contrib_infra)) contrib_infra];

            synergy_gran_infra_pairs = [synergy_gran_infra_pairs; ...
						synergy_gran(k)*ones(size(synergy_infra)) synergy_infra];

         end % (for k = 1:length(contrib_gran) )

      end % (if)

   end % (if)


   % process supra-infra pairs
   indsupra = find( position <= supra(2) );

   if ( ~isempty(indsupra) )

      contrib_supra = contrib(indsupra);
      synergy_supra = synergy(indsupra);


      % ----- Infragranular Layer Data -----
      indinfra = find(position >= infra(1) & position <= infra(2) );

      if ( ~isempty(indinfra) )

         contrib_infra = contrib(indinfra); 
			contrib_infra = contrib_infra(:);

         synergy_infra = synergy(indinfra); 
			synergy_infra = synergy_infra(:);

         for k = 1:length(contrib_supra)

            contrib_supra_infra_pairs = [contrib_supra_infra_pairs; ...
						contrib_supra(k)*ones(size(contrib_infra)) contrib_infra];

            synergy_supra_infra_pairs = [synergy_supra_infra_pairs; ...
						synergy_supra(k)*ones(size(synergy_infra)) synergy_infra];

         end % (for k = 1:length(asi_gran) )

      end % (if)

   end % (if)

end % (for i)



comb{1} = 'gran_supra';
comb{2} = 'gran_infra';
comb{3} = 'supra_infra';

titlecomb{1} = 'Gran-Supra';
titlecomb{2} = 'Gran-Infra';
titlecomb{3} = 'Supra-Infra';


%--------------------------------------------------------------------------------
%          Plot Contrib and Synergy Differences for all pair-wise combinations
%--------------------------------------------------------------------------------

close all;

f1 = figure;

for i = 1:length(comb)

   % Plot Contribution differences
   %-----------------------------------------
   eval(['pairs = contrib_' comb{i} '_pairs;']);
   newpairs = zeros(size(pairs));
   newpairs(:,1) = pairs(:,1);
   newpairs(:,2) = pairs(:,2);
   contrib_diff = newpairs(:,1) - newpairs(:,2);
   diff_median = median(contrib_diff);
   diff_mean = mean(contrib_diff);

	subplot(2,3,i);

   edges = linspace(-50,50,17);
   xtick = -50:10:50; %-2:0.5:2;
   xticklabel = xtick;
   xlim = [-50 50];

   num = histc(contrib_diff, edges);
   hb = bar(edges, num, 'histc');
   set(hb, 'facecolor', [0.8 0.8 0.8]);
   set(gca,'xlim', xlim);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   ylim = [0 1.1*max(num)];
   set(gca,'ylim', ylim);
   set(gca, 'tickdir', 'out');
	set(gca,'ticklength', [0.02 0.02]);
   hold on;
   plot([0 0], ylim, 'k-', 'linewidth', 1);
   [pvalue, ha, stats] = signrank( newpairs(:,1), newpairs(:,2) );
   title(sprintf([titlecomb{i} ',md=%.3f,mn=%.3f,n=%.0f,p=%.4f'], ...
		diff_median, diff_mean, size(newpairs,1), pvalue));
   if ( i == 1 ), ylabel(sprintf('MID1 Contribution Differences\n\n# Pairs')); end
	if ( i == 2 ), xlabel('MID1 Contribution Differences'); end;


   % Plot Synergy differences
   %-----------------------------------------
   eval(['pairs = synergy_' comb{i} '_pairs;']);
   newpairs = zeros(size(pairs));
   newpairs(:,1) = pairs(:,1);
   newpairs(:,2) = pairs(:,2);
   synergy_diff = newpairs(:,1) - newpairs(:,2);
   diff_median = median(synergy_diff);
   diff_mean = mean(synergy_diff);


	subplot(2,3,i+3);

   edges = linspace(-100,100,17);
   xtick = -100:20:100; %-2:0.5:2;
   xticklabel = xtick;
   xlim = [-100 100];

   num = histc(synergy_diff, edges);
   hb = bar(edges, num, 'histc');
   set(hb, 'facecolor', [0.8 0.8 0.8]);
   set(gca,'xlim', xlim);
   set(gca,'xtick', xtick, 'xticklabel', xticklabel);
   ylim = [0 1.1*max(num)];
   set(gca,'ylim', ylim);
   set(gca, 'tickdir', 'out');
	set(gca,'ticklength', [0.02 0.02]);
   hold on;
   plot([0 0], ylim, 'k-', 'linewidth', 1);
   [pvalue, ha, stats] = signrank( newpairs(:,1), newpairs(:,2) );
   title(sprintf([titlecomb{i} ',md=%.3f,mn=%.3f,n=%.0f,p=%.4f'], ...
		diff_median, diff_mean, size(newpairs,1), pvalue));
   if ( i == 1 ), ylabel(sprintf('Synergy Differences\n\n# Pairs')); end;
	if ( i == 2 ), xlabel('Synergy Differences'); end;

end % for i = 1:length(comb)
orient landscape;
print_mfilename(mfilename);

% Having to resize the figure every time is annoying, so I set it to a good
% size here
figure(f1);
set(gcf,'position', [100 200 1100 700]);

return






