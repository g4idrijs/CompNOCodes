/*********************************************************************
 * NAME     : motion.c
 * ABSTRACT : Managing motion compensation in files.
 * CREATED  : 07 Sep 2000, Svetoslav Nikolov
 *********************************************************************/
 
#include "../h/focus.h"
#include "../h/sys_params.h" 
#include "../h/geometry.h"
#include "../h/error.h"

#include <math.h>


/**********************************************************************
 * FUNCTION : delay_line_linear
 *
 * ABSTRACT : Delay the samples of a whole line, using linear
 *            interpolation between the samples.
 *
 * ARGUMENTS: sys - Pointer to a structure with the parameters of the 
 *                  system
 *            times - Array with times after which the next delay is valid
 *            delays - Array with the delays
 *            no_delays - Length of "times" and "delays"
 *            src_no_samples-  The length of the input signal
 *            src_start_time - The starting time of the input signal
 *            dest_start_time - Starting time of the output signal
 *            dest_no_samples - Number of samples in the output signal
 *
 * RETURNS :  Pointer to an array with the output signal. The array is 
 *            allocated inside the function
 **********************************************************************/

double* delay_line_linear(TSysParams *sys, double *times, double *delays,
                          ui32 no_delays, double *src, ui32 src_no_samples,
                          double src_start_time, double dest_start_time,
                          ui32 dest_no_samples)
{
  ui32 oi;               /*  Output absolute sample  */
  ui32 ii;               /*  Input sample*/
  ui32 o_abs_s;          /*  Absolute time output sample */
  ui32 i_start_s;        /* The sample number corresponding to src_start_time*/
  ui32 delay_no = 0;
  
  double * out;
  double delay;          /* Current delay              */
  double sample_delay;   /* The delay in samples       */
  int int_delay;         /* Integer delay              */
  double a1, a2;         /* Interpolation coefficients */
  double cur_time;

  PFUNC

  if (sys == NULL){
    errprintf("%s", "sys = NULL \n");
    goto dll_fail_1;
  }
  
  if (times == NULL){
     errprintf("%s","times = NULL \n");
     goto dll_fail_1;
  }
  
  if (delays == NULL){
     errprintf("%s","delays = NULL \n");
     goto dll_fail_1;
  }
  
  if (src == NULL){
     errprintf("%s"," src = NULL \n");
     goto dll_fail_1;
  }
  

  
  out = (double*)calloc(dest_no_samples, sizeof(double));
  if (out == NULL){
     errprintf("%s","Cannot allocate memory for the output line \n");
     goto dll_fail_1;
  }
  
  
  delay = delays[0];
  
  sample_delay = delay * sys->fs;

  int_delay = (int) floor(sample_delay);
    
  a1 = sample_delay - int_delay;
  a2 = 1.0 - a1;

  o_abs_s = (int)floor(dest_start_time*sys->fs);
  i_start_s = (int)floor(src_start_time*sys->fs);  
  dprintf("src_no_samples :%d \n", src_no_samples);
  for (oi = 0; oi < dest_no_samples; oi ++, o_abs_s++){
     cur_time = o_abs_s/sys->fs;
     
     /*
      *   Change the delays if necessary
      */
     if (delay_no < no_delays - 1){
        if (cur_time > times[delay_no + 1]){
  	        dprintf("%s\n","Changing the delay");
           delay_no ++;
           delay = delays[delay_no];
           sample_delay = delay * sys->fs;
           int_delay = (int) floor(sample_delay);
           a1 = sample_delay - int_delay;
           a2 = 1.0 - a1;
        }
     }
     ii = o_abs_s - int_delay - i_start_s;
     
     if (ii < src_no_samples && (ii-1) < src_no_samples)
        out[oi] = a2*src[ii] + a1*src[ii-1];
  }
  return out;
  
/* 
 *  Error handling part
 */  
dll_fail_1:
  return NULL;
}

/*********************************************************************
 * FUNCTION : delay_line_filter
 * ABSTRACT : Delay the samples of a whole line, using a filter bank
 *********************************************************************/
 
double* delay_line_filter(TSysParams *sys, TFilterBank* fb, 
                          double *times, double *delays,
                          ui32 no_delays, double *src, ui32 src_no_samples,
                          double src_start_time, double dest_start_time,
                          ui32 dest_no_samples)
{
  
  ui32 oi;
  ui32 ii;
 
  ui32 o_abs_s;          /*  Absolute time output sample */
  ui32 i_start_s;        /* The sample number corresponding to src_start_time*/
  ui32 delay_no = 0;
  ui32 jj;
  double * out;
  double delay;          /* Current delay              */
  double filter_delay;   /* The constant delay, introduced by the first coefficient of the filter*/
  double sample_delay;   /* The delay in samples       */
  int int_delay;         /* Integer delay              */
  double cur_time;
  int bank_no;
  
  PFUNC
  

  
  if (sys == NULL){
    errprintf("%s", "sys = NULL \n");
    goto dlf_fail_1;
  }
  
  if (fb == NULL){
     errprintf("%s", "The pointer to the filter bank is NULL \n");
     goto dlf_fail_1;
  }
  
  if (times == NULL){
     errprintf("%s","times = NULL \n");
     goto dlf_fail_1;
  }
  
  if (delays == NULL){
     errprintf("%s","delays = NULL \n");
     goto dlf_fail_1;
  }
  
  if (src == NULL){
     errprintf("%s","src = NULL \n");
     goto dlf_fail_1;
  }
  
  
  if (fb->Nf == 0 || fb->Ntaps == 0){
     errprintf("%s","There is no filter bank set \n");
     goto dlf_fail_1;
  }


  out = (double*)calloc(dest_no_samples, sizeof(double));
  if (out == NULL){
     errprintf("%s","Cannot allocate memory for the output line \n");
     goto dlf_fail_1;
  }
  
  filter_delay = (fb->Ntaps/2.0 + 1.0/(2.0*fb->Nf) - 1) / sys->fs;
  
  src_start_time -= filter_delay;
  /*  dest_start_time += filter_delay; */
  delay = delays[0];
  
  sample_delay = delay * sys->fs;

  int_delay = (int) floor(sample_delay);
      
   
  o_abs_s = (int)floor(dest_start_time*sys->fs);
  i_start_s = (int)floor(src_start_time*sys->fs);  
  dprintf("src_no_samples :%d \n", src_no_samples);
  bank_no = (int)((sample_delay - int_delay)*fb->Nf);
  /* bank_no = fb->Nf - bank_no -1;*/
  
  dprintf("sample_delay %f \n", sample_delay);
  dprintf("int_delay %d \n", int_delay);
  dprintf("Filter bank #%d \n", bank_no);
  for (oi = 0; oi < dest_no_samples; oi ++, o_abs_s++){
     cur_time = o_abs_s/sys->fs;
     
     /*
      *   Change the delays if necessary
      */
     if (delay_no < no_delays - 1){
        if (cur_time > times[delay_no + 1]){
           delay_no ++;
           delay = delays[delay_no];
           sample_delay = delay * sys->fs;
           int_delay = (int) floor(sample_delay);
        }
     }
    
     /*ii = o_abs_s - int_delay - i_start_s+1;*/
     ii = o_abs_s - int_delay - i_start_s + 1;
     if (ii < src_no_samples && (ii-fb->Ntaps) < src_no_samples){
        out[oi] = 0;
        for (jj = 0; jj < fb->Ntaps; jj++){
           out[oi] += fb->bank[bank_no][fb->Ntaps-jj+1]*src[ii - jj];
        }
     }
  }
  return out;

  
dlf_fail_1:
   return  NULL;
}





/**********************************************************************
 * FUNCTION  : delay_and_add
 * ABSTRACT  : 
 **********************************************************************/
 
double *delay_and_add(TSysParams *sys, TFilterBank* fb, 
                      double *times1, double *delays1,
                      ui32 no_delays1, double *src1,
                      double *times2, double *delays2,
                      ui32 no_delays2, double *src2,
                      double start_time, ui32 no_samples)
{
   
   double *out;
   



   PFUNC
   if (sys == NULL){
      errprintf("%s", "'sys' == NULL \n");
      return NULL;
   }  
   
   if (fb == NULL){
      errprintf("%s", "'filter bank' == NULL \n");
      return NULL;
   }
   
   if (times1 == NULL){
      errprintf("%s", "times1 == NULL \n");
      return NULL;
   }
   
   if (delays1 == NULL){
      errprintf("%s", "'delays1' == NULL \n");
      return NULL;
   }
   
   if (src1 == NULL){
      errprintf("%s", "'src1' == NULL \n");
      return NULL;
   }
   
   if (times2 == NULL){
      errprintf("%s","'times2' == NULL \n");
      return NULL;
   }
   
   if (delays2 == NULL){
      errprintf("%s","'delays2' == NULL \n");
      return NULL;
   }
   
   if (src2 == NULL){
       errprintf("%s","'src2' == NULL \n");
       return NULL;
   }
   
   
   out = (double*)calloc(no_samples, sizeof(double));
   if (out == NULL){
      errprintf("%s", "Can not allocate memory for output line \n");
      return NULL;
   }
   return out;
   
}
