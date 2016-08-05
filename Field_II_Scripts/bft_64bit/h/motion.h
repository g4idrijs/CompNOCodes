#ifndef __motion_h
 #define __motion_h
/**********************************************************************
 * NAME     : motion.h
 * ABSTRACT : Header file for the motion compensation functions.
 **********************************************************************/
 
#ifdef __cpluscplus 
  extern "C"{
#endif
double* delay_line_linear(TSysParams *sys, double *times, double *delays,
                          ui32 no_delays, double *src, ui32 src_no_samples,
                          double src_start_time, double dest_start_time,
                          ui32 dest_no_samples);
                          
double* delay_line_filter(TSysParams *sys, TFilterBank* fb, 
                          double *times, double *delays,
                          ui32 no_delays, double *src, ui32 src_no_samples,
                          double src_start_time, double dest_start_time,
                          ui32 dest_no_samples);

#ifdef __cplusplus
  };
#endif
 
#endif
