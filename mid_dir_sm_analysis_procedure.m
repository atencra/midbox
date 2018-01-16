function mid_dir_sm_analysis_procedure(varargin)

library('tmpmidbox');

options = struct('batch', 0);

options = input_options(options, varargin);

stimfolder = 'I:\SM_MIDs\Stimuli';
sprfile = 'dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8.spr';

if ( options.batch )

    d = dir('.');
    isub = [d(:).isdir]; %# returns logical vector
    nameFolds = {d(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    outerfolder = pwd;

    for i = 1:length(nameFolds)

        cd(nameFolds{i});

        mid_dir_sm_file_names_to_file_struct('savepath', '.');

        mid_dir_sm_file_struct_to_filters;

        mid_dir_plot_filter_fio;

        mid_dir_plot_filter_fio2d;

        mid_dir_filter_to_fio_info('stimfolder', stimfolder, 'batch', 1);

        cd(outerfolder);

    end

else

    mid_dir_file_names_to_file_struct('savepath', '.');

    mid_dir_file_struct_to_filters;

    mid_dir_plot_filter_fio;

    mid_dir_plot_filter_fio2d;

    mid_dir_filter_to_fio_info('stimfolder', stimfolder);

end


mid_dir_plot_train_test_filters_nonlinearities.m

mid_sm_plot_projinfo(projinfo);

return;




void  meanresp(double *sta,double *stimuli, unsigned long Nn,unsigned long fsize,int locator[],unsigned long Nspikes,  unsigned long Ntrials){
  unsigned long i,k;
  double temp;
  for(i=1;i<=Nn;i++)   sta[i]=0;
  double rbar;
  rbar=(double)Nspikes/(double)Ntrials;
  for(k=1;k<=Ntrials;k++){
    temp= ((double)locator[k]-rbar)/(double)Ntrials;
    for(i=1;i<=Nn;i++){
      sta[i]+=stimuli[(k-1)*fsize+i]*temp;
    }
  }
}










