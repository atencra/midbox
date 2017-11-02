% Script to make information calculations for the MID data.
%
% caa 1/4/07


% sprfile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min_6_carriers_per_octave.spr'
% 
% load 2002-8-26-site16-2414um-20db-dmr1-fs17857-spkcmb.mat
% load nonlinearity_params_516.mat
% keep mid spk trigger sprfile

% Need to replace the following with something similar, but for feature
% selectivity. Maybe:

calculate_mid2_fsi(mid, spk, trigger, sprfile);

% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);

% save projection_information_data_516 proj info



% 
% 
% load 2002-8-27-site17-2398um-20db-dmr1-fs17857-spkcmb.mat                                 
% load nonlinearity_params_517.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_517 proj info infodata
% 
% 
% load 2002-8-27-site19-2400um-20db-dmr1-fs17857-spkcmb.mat
% load nonlinearity_params_519.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_519 proj info infodata
% 
% 
% load 2003-3-5-site32-2357um-30db-dmr2-fs24038-spkcmb.mat
% load nonlinearity_params_532.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_532 proj info infodata
% 
% 
% load 2003-3-5-site34-2371um-30db-dmr1-fs24038-spkcmb.mat
% load nonlinearity_params_534.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_534 proj info infodata
% 
% 
% load 2003-3-5-site37-2337um-30db-dmr1-fs24038-spkcmb.mat
% load nonlinearity_params_537.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_537 proj info infodata
% 
% 
% load 2003-4-8-site15-2332um-20db-dmr1-fs27173-spkcmb.mat
% load nonlinearity_params_515.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_515 proj info infodata



% sprfile = 'dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min_6_carriers_per_octave.spr'
% 
% 
% load 2004-1-14-site2-2352um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_602.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_602 proj info infodata
% 
% 
% load 2004-1-14-site3-2388um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_603.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_603 proj info infodata
% 
% 
% load 2003-11-12-site4-2355um-30db-dmr1-fs20833-spkcmb.mat
% load nonlinearity_params_604.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_604 proj info infodata
% 
% 
% load 2004-1-14-site8-2395um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_608.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_608 proj info infodata
% 
% 
% 
% load 2004-1-14-site9-2389um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_609.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_609 proj info infodata
% 
% 
% load 2004-1-14-site11-2391um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_611.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_611 proj info infodata
% 
% 
% load 2004-1-14-site14-2392um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_614.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_614 proj info infodata
% 
% 
% load 2004-1-14-site16-2390um-30db-dmr1-fs18115-spkcmb.mat
% load nonlinearity_params_616.mat
% keep mid spk trigger sprfile
% proj = calculate_projection(mid, spk, trigger, sprfile);
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_616 proj info infodata




% load projection_information_data_515.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_515 proj info infodata
% 
% 
% load projection_information_data_516.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_516 proj info infodata
% 
% 
% load projection_information_data_517.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_517 proj info infodata
% 
% 
% load projection_information_data_519.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_519 proj info infodata
% 
% 
% load projection_information_data_532.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_532 proj info infodata
% 
% 
% load projection_information_data_534.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_534 proj info infodata
% 
% 
% load projection_information_data_537.mat
% keep proj
% info = calculate_info(proj);
% infodata = plot_proj_info(info);
% save projection_information_data_537 proj info infodata



load projection_information_data_602.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_602 proj info infodata


load projection_information_data_603.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_603 proj info infodata


load projection_information_data_604.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_604 proj info infodata


load projection_information_data_608.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_608 proj info infodata


load projection_information_data_609.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_609 proj info infodata


load projection_information_data_611.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_611 proj info infodata


load projection_information_data_614.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_614 proj info infodata


load projection_information_data_616.mat
keep proj
info = calculate_info(proj);
infodata = plot_proj_info(info);
save projection_information_data_616 proj info infodata


exit




