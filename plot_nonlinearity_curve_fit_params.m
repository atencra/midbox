function plot_nonlinearity_curve_fit_params(params)
%
%
% You get "params" from the function
%
% [params, nmse] = plot_nonlinearity_params(fiopos, midpos)
%
% 


% fx = fx ./ .005; % 5 ms time bins - convert to firing rate
% a0 = [1 1 1 1];
% [a, resnorm] = lsqcurvefit(@sachs_abbas_func, a0, x, fx);
% params = [params; a(:)'];
% xfit = linspace(min(x), max(x), 100);
% fxfit = sachs_abbas_func(a, xfit);

%    subplot(3,3,i);
%    hold on;
%    plot(x,fx, 'ko', 'markerfacecolor', 'k');
%    plot(xfit,fxfit, 'k-');
%    xlim([-6 6]);



close all;

% Let's plot some example nonlinearity functions

% Vary the growth rate but keep everything else constant

xfit = linspace(-6, 10, 10000);

figure;

m = {'-' '--' '-' '--' '-' '--'};
c = [208	209 230; 166 189 219; 116 169 207; 54 144 192; 5 112 176; 3 78 123] ./ 256;
c = flipud(c);
% plot(unqdiff, uniquemean, [m{i}], 'linewidth', 2, 'color', c(i,:));


% Vary the growth rate
subplot(1,2,1);
hold on;
mr = 12;


a1 = [0 1 2 0];
fx1 = sachs_abbas_func(a1, xfit);
plot(xfit(fx1<mr), fx1(fx1<mr), [m{1}], 'linewidth', 2, 'color', c(1,:));
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);

a1 = [0 2 2 0];
fx1 = sachs_abbas_func(a1, xfit);
plot(xfit(fx1<mr), fx1(fx1<mr), [m{2}], 'linewidth', 2, 'color', c(2,:));
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);

a2 = [0 3 2 0];
fx2 = sachs_abbas_func(a2, xfit);
plot(xfit(fx2<mr), fx2(fx2<mr), [m{3}], 'linewidth', 2, 'color', c(3,:));

a3 = [0 5 2 0];
fx = sachs_abbas_func(a3, xfit);
plot(xfit(fx<mr), fx(fx<mr), [m{4}], 'linewidth', 2, 'color', c(4,:));
title('Growth Rate');

a3 = [0 7 2 0];
fx3 = sachs_abbas_func(a3, xfit);
plot(xfit(fx3<mr), fx3(fx3<mr), [m{5}], 'linewidth', 2, 'color', c(5,:));
title('Growth Rate');

ylim([0 13]);
xlim([-5 7]);
legend('1', '2', '3', '5', '7', 'location', 'northwest');


% % Vary the width
% subplot(1,3,2);
% hold on;
% 
% a1 = [0 2 1 0];
% fx1 = sachs_abbas_func(a1, xfit);
% plot(xfit, fx1, 'k-', 'linewidth', 2);
% set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
% 
% a2 = [0 2 3 0];
% fx2 = sachs_abbas_func(a2, xfit);
% plot(xfit, fx2, 'c-', 'linewidth', 2);
% 
% a3 = [0 2 5 0];
% fx3 = sachs_abbas_func(a3, xfit);
% plot(xfit, fx3, 'b-', 'linewidth', 2);
% temp = [fx1(:); fx2(:); fx3(:)];
% ylim([0 1.05*max(temp)]);
% xlim(1.05*[min(xfit) max(xfit)]);
% title('Width');
% 
% legend('1', '3', '5', 'location', 'northwest');



% Vary the transition point
subplot(1,2,2);
hold on;


a1 = [0 2 2 -2];
fx = sachs_abbas_func(a1, xfit);
plot(xfit(fx<mr), fx(fx<mr), [m{1}], 'linewidth', 2, 'color', c(1,:));

a1 = [0 2 2 0];
fx = sachs_abbas_func(a1, xfit);
plot(xfit(fx<mr), fx(fx<mr), [m{2}], 'linewidth', 2, 'color', c(2,:));

a2 = [0 2 2 2];
fx = sachs_abbas_func(a2, xfit);
plot(xfit(fx<mr), fx(fx<mr), [m{3}], 'linewidth', 2, 'color', c(3,:));


a3 = [0 2 2 4];
fx = sachs_abbas_func(a3, xfit);
plot(xfit(fx<mr), fx(fx<mr), [m{4}], 'linewidth', 2, 'color', c(4,:));
ylim([0 13]);
xlim([-5 8]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition (SD)');

legend('-2', '0', '2', '4', 'location', 'northwest');

return;



figure;


a1 = 0;
a2 = [1 3 5 7];
a3 = 2;
a4 = [-2 0 2 4];

for i = 1:length(a2)

   for j = 1:length(a4)


      subplot(length(a2), length(a4), (i-1)*length(a2)+ j );
      a = [a1 a2(i) a3 a4(j)];
      fx = sachs_abbas_func(a, xfit);
      plot(xfit, fx, 'k-', 'linewidth', 2);
      set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
      xlim([-6.25 6.25]);
      ylim([0 70]); %1.05 * max(fx)]);

   end % (for i)

end % (for j)

set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);






figure;

subplot(2,3,1);
hist(params(:,2), 25);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Growth Rate');

subplot(2,3,2);
hist(params(:,3), 25);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Width');

subplot(2,3,3);
hist(params(:,4), 25);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition (SD)');

subplot(2,3,4);
plot(params(:,2), params(:,3), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Width Vs. Growth Rate');
xlabel('Growth Rate');
ylabel('Width');

subplot(2,3,5);
plot(params(:,2), params(:,4), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition Vs. Growth Rate');
xlabel('Growth Rate');
ylabel('Transition (SD)');

subplot(2,3,6);
plot(params(:,3), params(:,4), 'ko', 'markerfacecolor', 'k');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
title('Transition Vs. Width');
xlabel('Width');
ylabel('Transition (SD)');

return;
