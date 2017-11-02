function plot_dmr_params_from_mid_projections(x0spktrain_params, x0train_params)


close all;

% training set 1

figure;

for i = 1:length(x0spktrain_params)

x0params1 = x0spktrain_params{i};
tmf = x0params1(:,1);
smf = x0params1(:,2);


subplot(4,2, (i-1)*2+1 );
hist(tmf, 50);
title( sprintf('TMF: Training set %.0f', i) );

subplot(4,2, i*2 );
hist(smf, 50);
title( sprintf('SMF: Training set %.0f', i) );

end



figure;

for i = 1:length(x0spktrain_params)

x0params1 = x0spktrain_params{i};
tmf = x0params1(:,1);
smf = x0params1(:,2);

subplot(2,2, i );
plot(tmf, smf, '.');
title( sprintf('SMF vs TMF: Training set %.0f', i) );

end






figure;

for i = 1:length(x0train_params)

x0params1 = x0train_params{i};
tmf = x0params1(:,1);
smf = x0params1(:,2);


subplot(4,2, (i-1)*2+1 );
hist(tmf, 50);
title( sprintf('Prior TMF: Training set %.0f', i) );

subplot(4,2, i*2 );
hist(smf, 50);
title( sprintf('Prior SMF: Training set %.0f', i) );

end

return;

