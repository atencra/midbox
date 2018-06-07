function fio = mid_filtstr_to_nonlinearity(filtstr, stimulus)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fio = filtstr;

for i = 1:length(filtstr)

    
    fprintf('\nProcessing unit %d of %d\n', i, length(filtstr))
   
    locator = filtstr(i).locator; 
    x0 = filtstr(i).x0;

    %tempstim = stim_mat(x0:x0+25-1,:);

    stimulus_reduced = stimulus(x0:x0+25-1,:);
    
    sta = filtstr(i).v_sta;
    mid1 = filtstr(i).v1;
    mid2 = filtstr(i).v2;
   
    try 
        %[xprior, xposterior] = ne_sta_stimulus_projection(sta, locator(i,:), tempstim);
        %[px,pspk,pxspk,xbinedges] = calc_px_pspk_pxspk(xprior,xposterior, 15);


        [xprior, xposterior] = mid_filter_locator_stimulus_to_projection(sta, locator, stimulus_reduced);
        [xbins, pspk, px, pxspk, pspkx] = mid_projection_to_nonlinearity(xprior, xposterior);
        fio(i).sta_xbins = xbins;
        fio(i).sta_pspk = pspk;
        fio(i).sta_px = px;
        fio(i).sta_pxspk = pxspk;
        fio(i).sta_pspkx = pspkx;
        clear('xprior', 'xposterior', 'xbins', 'px', 'pxspk', 'pspkx');

        
        %[xprior, xposterior] = ne_sta_stimulus_projection(mid1, locator(i,:), tempstim);
        %[px,pspk,pxspk,xbinedges] = calc_px_pspk_pxspk(xprior,xposterior, 15);
        %filtstr(i).mid1_xbins = edge2center(xbinedges);
        %filtstr(i).mid1_fio = pspk .* pxspk ./ px;
        
        [xprior, xposterior] = mid_filter_locator_stimulus_to_projection(mid1, locator, stimulus_reduced);
        [xbins, pspk, px, pxspk, pspkx] = mid_projection_to_nonlinearity(xprior, xposterior);
        fio(i).mid1_xbins = xbins;
        fio(i).mid1_pspk = pspk;
        fio(i).mid1_px = px;
        fio(i).mid1_pxspk = pxspk;
        fio(i).mid1_pspkx = pspkx;
        clear('xprior', 'xposterior', 'xbins', 'px', 'pxspk', 'pspkx');

        
        
        %[xprior, xposterior] = ne_sta_stimulus_projection(mid2, locator(i,:), tempstim);
        %[px,pspk,pxspk,xbinedges] = calc_px_pspk_pxspk(xprior,xposterior, 15);
        %filtstr(i).mid2_xbins = edge2center(xbinedges);
        %filtstr(i).mid2_fio = pspk .* pxspk ./ px;

        [xprior, xposterior] = mid_filter_locator_stimulus_to_projection(mid2, locator, stimulus_reduced);
        [xbins, pspk, px, pxspk, pspkx] = mid_projection_to_nonlinearity(xprior, xposterior);
        fio(i).mid2_xbins = xbins;
        fio(i).mid2_pspk = pspk;
        fio(i).mid2_px = px;
        fio(i).mid2_pxspk = pxspk;
        fio(i).mid2_pspkx = pspkx;
        clear('xprior', 'xposterior', 'xbins', 'px', 'pxspk', 'pspkx');

    catch

        fio(i).sta_xbins = [];
        fio(i).sta_pspk = [];
        fio(i).sta_px = [];
        fio(i).sta_pxspk = [];
        fio(i).sta_pspkx = [];
        
        fio(i).mid1_xbins = [];
        fio(i).mid1_pspk = [];
        fio(i).mid1_px = [];
        fio(i).mid1_pxspk = [];
        fio(i).mid1_pspkx = [];

        fio(i).mid2_xbins = [];
        fio(i).mid2_pspk = [];
        fio(i).mid2_px = [];
        fio(i).mid2_pxspk = [];
        fio(i).mid2_pspkx = [];

    end

     
end

fprintf('\n')

return;




