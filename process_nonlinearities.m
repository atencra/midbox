function fio_params = process_nonlinearities(mid)
% process_nonlinearities - analyze nonlinearities using curve fits
%    and summary statistics
%
% caa 4/19/06


for i = 1:length(mid)

   sta_fio_params = sta_nonlinearity_params(mid(i).rpx1pxpxt_sta);
   v1_fio_params = v1v2_nonlinearity_params(mid(i).rpdx1x2px_pxt_2, 1);
   v2_fio_params = v1v2_nonlinearity_params(mid(i).rpdx1x2px_pxt_2, 2);

   % Now compute some simple correlation coefficients
   fx_sta = mid(i).rpx1pxpxt_sta.ior_mean;
   fx_v1 = mid(i).rpdx1x2px_pxt_2.ior1_mean;
   fx_v2 = mid(i).rpdx1x2px_pxt_2.ior2_mean;

   [r_sta_v1, p_sta_v1] = corrcoef(fx_sta, fx_v1);
   [r_sta_v2, p_sta_v2] = corrcoef(fx_sta, fx_v2);
   [r_v1_v2, p_v1_v2] = corrcoef(fx_v1, fx_v2);


   % Assign the output data
   fio_params(i).location = mid(i).location;
   fio_params(i).unit = mid(i).unit;

   fio_params(i).sta = sta_fio_params;
   fio_params(i).v1 = v1_fio_params;
   fio_params(i).v2 = v2_fio_params;

   fio_params(i).r_sta_v1 = r_sta_v1(1,2);
   fio_params(i).p_sta_v1 = p_sta_v1(1,2);

   fio_params(i).r_sta_v2 = r_sta_v2(1,2);
   fio_params(i).p_sta_v2 = p_sta_v2(1,2);

   fio_params(i).r_v1_v2 = r_v1_v2(1,2);
   fio_params(i).p_v1_v2 = p_v1_v2(1,2);

end % (for)

return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   Function Declarations   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fio_params = sta_nonlinearity_params(fio_sta_data)

xdata = fio_sta_data.x;
fx = fio_sta_data.ior_mean;
fx = fx + 0.001;
fxmax = max(fx);
px = fio_sta_data.px_mean;

logistic = inline('x(1) + fxmax * ( 1 + x(2) * exp( -x(3) .* xdata ) ).^(-1)', 'x', 'xdata', 'fxmax');
logisticparams = lsqcurvefit(logistic, [mean(fx) 1 1], xdata, fx, [], [], [], fxmax);

sachsabbas = inline('x(1) +  x(2)*x(3).*log10( 1 + 10.^( (xdata - x(4)) / x(3) ) )', 'x', 'xdata', 'fxmax');
saparams = lsqcurvefit(sachsabbas, [mean(fx) 1 1 1], xdata, fx, [], [], [], fxmax);

nakarushton = inline('x(1) +  fxmax * xdata.^x(2) ./ ( xdata.^x(2) + x(3).^x(2) ) ', 'x', 'xdata', 'fxmax');
nrparams = lsqcurvefit(nakarushton, [mean(fx) 1 1], xdata, fx, [], [], [], fxmax);

dxi = 0.01;
xi = xdata(1):dxi:xdata(end);
logisticfit = logistic(logisticparams, xi, fxmax);
safit = sachsabbas(saparams, xi, fxmax);
nrfit = nakarushton(nrparams, xi, fxmax);

fxfit = interp1(xdata, fx, xi, 'linear');
pxfit = interp1(xdata, px, xi, 'linear');

eflogf = dxi * trapz( fxfit .* log2(fxfit) .* pxfit ); % E( f log2 f )
ef = dxi * trapz( fxfit .* pxfit ); % E( f )
ef2 = dxi * trapz( fxfit.^2 .* pxfit); % E( f^2 )

% Now get the derivative of f
xd = unique( sort( [xi xdata] ) );
yd = interp1(xdata, fx, xd, 'linear');

for i = 1:length(xdata)
   index = find( xd == xdata(i) );
   if ( i == 1 ) % second point minus first point
      fxderiv(i) = ( yd(index+1) - yd(index) ) / ( xd(index+1) - xd(index) );
   elseif ( i > 1 & i < length(xdata) ) % i+1 minus i-1 point
      fxderiv(i) = ( yd(index+1) - yd(index-1) ) / ( xd(index+1) - xd(index-1) );
   else % last point minus second to last point
      fxderiv(i) = ( yd(index) - yd(index-1) ) / ( xd(index) - xd(index-1) );
   end
end % (for)

fxderivfit = interp1(xdata, fxderiv, xi, 'linear');
efdf = dxi * trapz( fxderivfit ./ fxfit .* pxfit); % E( f' / f )

%[ef ef2 eflogf efdf]
clf;
hold on;
plot(xdata, fx, 'go');
plot(xi, fxfit, 'r-');
plot(xdata, fxderiv, 'bo');
plot(xi, fxderivfit, 'c-');
plot([xdata(1) xdata(end)], [0 0], 'k-');
%pause;

fio_params.x = xdata;
fio_params.fx = fx;
fio_params.dfx = fxderiv;
fio_params.ef = ef;
fio_params.ef2 = ef2;
fio_params.eflogf = eflogf;
fio_params.efdf = efdf;

return;




function fio_params = v1v2_nonlinearity_params(fio_data, filter)

if ( filter == 1 )
   xdata = fio_data.x1;
   fx = fio_data.ior1_mean;
   fxmax = max(fx);
   px = fio_data.px1_mean;
elseif ( filter == 2 )
   xdata = fio_data.x2;
   fx = fio_data.ior2_mean;
   fxmax = max(fx);
   px = fio_data.px2_mean;
else
   error('Wrong filter.');
end

fx = fx + 0.001;

dxi = 0.01;
xi = xdata(1):dxi:xdata(end);
fxfit = interp1(xdata, fx, xi, 'linear');
pxfit = interp1(xdata, px, xi, 'linear');

eflogf = dxi * trapz( fxfit .* log2(fxfit) .* pxfit ); % E( f log2 f )
ef = dxi * trapz( fxfit .* pxfit ); % E( f )
ef2 = dxi * trapz( fxfit.^2 .* pxfit); % E( f^2 )


% Now get the derivative of f(x)
xd = unique( sort( [xi xdata] ) );
yd = interp1(xdata, fx, xd, 'linear');

for i = 1:length(xdata)
   index = find( xd == xdata(i) );
   if ( i == 1 ) % second point minus first point
      fxderiv(i) = ( yd(index+1) - yd(index) ) / ( xd(index+1) - xd(index) );
   elseif ( i > 1 & i < length(xdata) ) % i+1 minus i-1 point
      fxderiv(i) = ( yd(index+1) - yd(index-1) ) / ( xd(index+1) - xd(index-1) );
   else % last point minus second to last point
      fxderiv(i) = ( yd(index) - yd(index-1) ) / ( xd(index) - xd(index-1) );
   end
end % (for)

fxderivfit = interp1(xdata, fxderiv, xi, 'linear');
efdf = dxi * trapz( fxderivfit ./ fxfit .* pxfit); % E( f' / f )

%[ef ef2 eflogf efdf]
clf;
hold on;
plot(xdata, fx, 'go');
plot(xi, fxfit, 'r-');
plot(xdata, fxderiv, 'bo');
plot(xi, fxderivfit, 'c-');
plot([xdata(1) xdata(end)], [0 0], 'k-');
legend('f(x)', 'fit', 'df(x)', 'fit', 2);
%pause;

fio_params.x = xdata;
fio_params.fx = fx;
fio_params.dfx = fxderiv;
fio_params.ef = ef;
fio_params.ef2 = ef2;
fio_params.eflogf = eflogf;
fio_params.efdf = efdf;

return;




