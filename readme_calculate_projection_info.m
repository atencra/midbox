readme_calculate_projection_info

Description of steps used to calculate mutual information from
the STA and MIDs that Tatyana Sharpee computed.

The set of commands that need to be issued is:


   proj = calculate_projection(mid, spk, trigger, sprfile);
   info = calculate_info(proj);
   infodata = plot_proj_info(info);
   plot_infodata(infodata);



These commands are explained below:



(1) Calculate projections.
--------------------------------------------------------------------------

Run the following command

      proj = calculate_projection(mid, spk, trigger, sprfile)

This function is stored in

      C:\MATLAB65\work\mid\calc_info_code

mid and spk are struct arrays, while trigger is a vector.

mid holds the filters, while spk and trigger contain the spike trains
and stimulus trigger times. The struct array mid is saved in files of the
form nonlinearity_params_*.mat

spk and trigger are saved in files such as 

      2003-4-8-site15-2332um-20db-dmr1-fs27173-spkcmb.mat  

All of these files are stored in directories, such as 
site515, site532, etc. in C:\MATLAB65\work\mid\

proj is a struct array holding all the projection data.

The file that holds the stimulus envelope, sprfile, is either  

      sprfile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min_DFt22_DFf7.spr';

   or

      sprfile = 'dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min_DFt10_DFf8.spr'

The main difference between the files is the sampling rate, and, thus,
the frequency extent.


(2) Calculate information.
--------------------------------------------------------------------------

Run the command

      info = calculate_info(proj)

to compute the mutual information for the STA and MIDs from the
projection values 


This function calculates mutual information for all the projection
values, and also over only portions of the data. Portions are used
since Tatyana calculated four estimates of the filters: one for
each 3/4 of the data. Thus the first estimate is for quarters 2-4;
the second for quarters 1, 3, 4; the third for quarters 1, 2, 4;
and the last estimate for quarters 1-3.

The mutual information is calculated over these 3/4 data set lengths.
Also, the mi is calculated, for each filter, over the quarter of the
data that was not used to calculate the filter. So for filter 1, for
which quarters 2-4 were used to estimate, the information from one
quarter would be that from quarter 1.

In the language of learning theory, the first information calculations
are over the training data, and the last calculations are over the 
test data.

Further, the information is calculated over different lengths of the data.
Thus, each data set that is used to calculate information is broken into
50, 60, 70, 80, 90, 92.5, 95, 97.5, and 100 of the data. The information
is calculated over these portions so that extrapolation procedures may
be used in the following analysis ...


(3) Get the extrapolated mutual information valus.
--------------------------------------------------------------------------

Run the command

      infodata = plot_proj_info(info)

This function takes the previous information calculations and extrapolates
the values to estimate the information based on infinite data set size.

The results of this procedure are stored in the struct array infodata.




(4) Plot the information.
--------------------------------------------------------------------------

Run the command

      plot_infodata(infodata)

To get a plot of the extrapolated information values.

Then function ends by plotting the information without extrapolation
to infinite data set size versus that obtained using the extrapolation 
routines.









