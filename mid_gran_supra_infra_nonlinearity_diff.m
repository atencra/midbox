function mid_gran_supra_infra_nonlinearity_diff(fiocell, midcell)
%mid_gran_supra_infra_nonlinearity_diff - plots temporal,
% spectral modulation transfer function widths at different positions 
% in the cortical microcircuit.
%
% The plot compares the average mtf width of layer 4 cells to the widths in
% layer 2/3, layer5, and layer 6.
%
% rtf_params_database : contains all the rtf_params struct arrays from all the
% data.
%
% exportfig(gcf,'fig678_layer4_layer23_layer5_layer6_btmf.eps', 'color', 'bw', 'width', 6.3, 'height', 3.94, 'fontmode', 'fixed', 'fontsize', 8);
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

asi_gran_supra_pairs = [];
asi_gran_infra_pairs = [];
asi_supra_infra_pairs = [];

sep_gran_supra_pairs = [];
sep_gran_infra_pairs = [];
sep_supra_infra_pairs = [];


% Define the layers as:

gran = [700 1100];
supra = [0 600];
infra = [1200 1900];


for i = 1:length(fiocell)

   fio_params = fiocell{i};
   mid_params = midcell{i};

	position = [];
	asi = [];
	sep = [];

	for j = 1:length(fio_params)

	   position(j) = mid_params(j).position;

	   % STA f(x) params
	   x = fio_params(j).sta.x;
		fx = fio_params(j).sta.fx;
		indexright = find(x>0);
		indexleft = find(x<0);
		right = sum( fx(indexright) - min(fx) );
		left = sum( fx(indexleft) - min(fx) );
		sta_fx_asi(j) = (right - left) / (right + left);


	   % MID1 f(x) params
		x = fio_params(j).v1.x;
		fx = fio_params(j).v1.fx;
		indexright = find(x>0);
		indexleft = find(x<0);
		right = sum( fx(indexright) - min(fx) );
		left = sum( fx(indexleft) - min(fx) );
		asi(j) = (right - left) / (right + left);


	   % MID2 f(x) params
	   x = fio_params(j).v2.x;
		fx = fio_params(j).v2.fx;
		indexright = find(x>0);
		indexleft = find(x<0);
		right = sum( fx(indexright) - min(fx) );
		left = sum( fx(indexleft) - min(fx) );
		v2_fx_asi(j) = (right - left) / (right + left);
 

		% MID12 f(x) params
	   fx1 = mid_params(j).rpdx1x2px_pxt_2.ior1_mean;
		fx2 = mid_params(j).rpdx1x2px_pxt_2.ior2_mean;
		fx1x2 = mid_params(j).rpdx1x2px_pxt_2.ior12_mean;
		[u,s,v] = svd(fx1x2);
		singvals = sum(s);
		eigvals = singvals .^ 2;
		sep(j) = 1 - eigvals(1) / (sum(eigvals)+eps);

	end % (for j)


	% Parameter we really care about are asi and sepindex


   % ==============================================================
   %       all possible pairs of bmf, width, and cutoff freq
   % ==============================================================

   % process pairs gran-supra, gran-infra

   indgran = find( position >= gran(1) & position <= gran(2) );

   if ( ~isempty(indgran) )

      asi_gran = asi(indgran);
      sep_gran = sep(indgran);

      % ----- Supragranular Layer Data -----
      indsupra = find( position <= supra(2) );

      if ( ~isempty(indsupra) )

	      asi_supra = asi(indsupra);
			asi_supra = asi_supra(:);

		   sep_supra = sep(indsupra);
			sep_supra = sep_supra(:);


         for k = 1:length(asi_gran)

            asi_gran_supra_pairs = [asi_gran_supra_pairs; ...
									asi_gran(k)*ones(size(asi_supra)) asi_supra];

            sep_gran_supra_pairs = [sep_gran_supra_pairs; ...
									sep_gran(k)*ones(size(sep_supra)) sep_supra];

         end % (for k = 1:length(btmf_4) )

      end % (if)


      % ----- Infragranular Layer Data -----
      indinfra = find( position >= infra(1) & position <= infra(2) );

      if ( ~isempty(indsupra) )

	      asi_infra = asi(indinfra);
			asi_infra = asi_infra(:);

		   sep_infra = sep(indinfra);
			sep_infra = sep_infra(:);

         for k = 1:length(asi_gran)

            asi_gran_infra_pairs = [asi_gran_infra_pairs; ...
									asi_gran(k)*ones(size(asi_infra)) asi_infra];

            sep_gran_infra_pairs = [sep_gran_infra_pairs; ...
									sep_gran(k)*ones(size(sep_infra)) sep_infra];

         end % (for k = 1:length(asi_gran) )

      end % (if)

   end % (if)


   % process supra-infra pairs
   indsupra = find( position <= supra(2) );

   if ( ~isempty(indsupra) )

      asi_supra = asi(indsupra);
      sep_supra = sep(indsupra);


      % ----- Infragranular Layer Data -----
      indinfra = find(position >= infra(1) & position <= infra(2) );

      if ( ~isempty(indinfra) )

         asi_infra = asi(indinfra); 
			asi_infra = asi_infra(:);

         sep_infra = sep(indinfra); 
			sep_infra = sep_infra(:);

         for k = 1:length(asi_supra)

            asi_supra_infra_pairs = [asi_supra_infra_pairs; ...
									asi_supra(k)*ones(size(asi_infra)) asi_infra];

            sep_supra_infra_pairs = [sep_supra_infra_pairs; ...
									sep_supra(k)*ones(size(sep_infra)) sep_infra];

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
%          Plot best TMF and tMTF width for all pair-wise combinations
%--------------------------------------------------------------------------------

close all;

f1 = figure;

for i = 1:length(comb)

   % Plot ASI differences
   %-----------------------------------------
   eval(['pairs = asi_' comb{i} '_pairs;']);
   newpairs = zeros(size(pairs));
   newpairs(:,1) = pairs(:,1);
   newpairs(:,2) = pairs(:,2);
   asi_diff = newpairs(:,1) - newpairs(:,2);
   diff_median = median(asi_diff);
   diff_mean = mean(asi_diff);

%    fprintf('ASI DIFF Median = %.4f\n', diff_median);
%    fprintf('ASI DIFF Mean = %.4f\n', diff_mean);

	subplot(2,3,i);

   edges = linspace(-2,2,17);
   xtick = -2:2; %-2:0.5:2;
   xticklabel = xtick;
   xlim = [-2 2];

   num = histc(asi_diff, edges);
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
   if ( i == 1 ), ylabel(sprintf('ASI Differences\n\n# Pairs')); end
	if ( i == 2 ), xlabel('MID1 ASI Differences'); end;


   % Plot SEP differences
   %-----------------------------------------
   eval(['pairs = sep_' comb{i} '_pairs;']);
   newpairs = zeros(size(pairs));
   newpairs(:,1) = pairs(:,1);
   newpairs(:,2) = pairs(:,2);
   sep_diff = newpairs(:,1) - newpairs(:,2);
   diff_median = median(sep_diff);
   diff_mean = mean(sep_diff);

%    fprintf('SEP DIFF Median = %.4f\n', diff_median);
%    fprintf('SEP DIFF Mean = %.4f\n', diff_mean);

	subplot(2,3,i+3);

   edges = linspace(-1,1,17);
   xtick = -1:0.5:1; %-2:0.5:2;
   xticklabel = xtick;
   xlim = [-1 1];

   num = histc(sep_diff, edges);
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
   if ( i == 1 ), ylabel(sprintf('Inseparability Differences\n\n# Pairs')); end;
	if ( i == 2 ), xlabel('2D Nonlinearity Inseparability Index Differences'); end;

end % for i = 1:length(comb)
orient landscape;
print_mfilename(mfilename);

% Having to resize the figure every time is annoying, so I set it to a good
% size here
figure(f1);
set(gcf,'position', [100 200 1100 700]);



return






