#ifndef __focus_h
  #define __focus_h
/*********************************************************************
 * NAME     : focus.h
 * ABSTRACT : Managing focus time lines
 * CREATED  : Jan-Feb 2000, Svetoslav Nikolov
 * MODIFIED : 04 Sep. 2000, Svetoslav Nikolov - added support for filter
 *                          banks.
 *********************************************************************/

#include "transducer.h" 
#include "sys_params.h"


#define MAX_SAMPLE_NO   65365

/*
 *  Definition of focus delay
 */
 
typedef struct delay{
   double time;      /* Time after which the delay is valid          */
   si32 *d;          /* The delay. One delay per channel             */
   double *a;        /* Weighting coefficient for linear interpolation */
}TDelay;



/*
 *    Definition of filter bank
 */
typedef struct filter_bank{
  double *  coefs;         /* Coefficients */
  double ** bank;          
  ui32 Nf;
  ui32 Ntaps;
}TFilterBank;



/*
 *   Definition of apodization
 */
typedef struct apodization{
  double time;      /* Time after which the associated apodization is valid */
  double* a;        /* Apodization coefficient. One per channel             */
} TApodization;




/*
 *   Focus time-line - collection of delays. 
 *   If the line is to be dynamically focused, then the delays are neglected,
 *   and the focus is re-calculated for every sample. The beam is defined 
 *   by the directional angles xy, and xz;
 */
 
typedef struct focus_time_line{
   ui32 no_times;          /* How many delays are associated with this line */
   ui32 dynamic;           /* Whether this is dynamic apodization           */
   ui32 pixel;             /* Whether the focusing is pixel based           */
   
   TPoint3D center;        /* Center of the focus.                          */
   TPoint3D *pixels;       /* Pixels in which to focus the data             */
   double dir_xz;          /* Direction in XZ                               */
   double dir_yz;          /* Direction in YZ                               */
   TTransducer* xdc;       /* Used in the dynamic focusing                  */
   TDelay *delay;          /* Array of delays. One entry per focal zone     */
}TFocusTimeLine;


/*
 *  Collection of focus time-lines makes one whole image
 */

typedef struct{
   ui32 no_focus_time_lines;
   TFocusTimeLine* ftl;
   ui32 use_filter_bank;   /* Whether to use filter bank for delays calculation */
   TFilterBank filter_bank;  /* Filter bank, used to calculate the delays  */
}TFocusLineCollection;


/*
 *  Apodization time-line - Collection of apodization values.
 */


typedef struct apodization_time_line{
  ui32 no_times;
  TApodization* a;
}TApoTimeLine;


/*
 *  Collection of apodization time lines - makes a whole image
 */
typedef struct{
   ui32 no_apo_time_lines;
   TApoTimeLine * atl;
}TApoLineCollection;



/**********************************************************************
 *                                                                    *
 *       Here follow the definitions of the exported functions        *
 *                                                                    *
 **********************************************************************/


#ifdef __cplusplus
  extern"C"{
#endif

void del_focus_time_line(TFocusTimeLine* p);

TFocusLineCollection* new_focus_line_collection();

void del_focus_line_collection(TFocusLineCollection* f);
void del_apo_line_collection(TApoLineCollection* alc); 

void set_no_lines(TApoLineCollection* alc, TFocusLineCollection *flc,
                                                      ui32 no_lines);
                                                      
void set_center_focus(TFocusLineCollection *flc, TPoint3D *p, ui32 line_no);

void set_dynamic_focus(TFocusLineCollection *flc, TTransducer* xdc,
                             ui32 line_no, double dir_xz, double dir_yz);

void set_focus_times(TFocusLineCollection *flc, TSysParams* sys, 
                     TTransducer* xdc, double* times, double *delays, 
                     ui32 no_times,  ui32 line_no);

void set_focus(TFocusLineCollection *flc, TSysParams* sys, 
                     TTransducer* xdc, double* times, TPoint3D *points, 
                     ui32 no_times,  ui32 line_no);

void set_focus_2way(TFocusLineCollection *flc, TSysParams* sys, 
                     TTransducer* xdc, double* times, TPoint3D *points, 
                     ui32 no_times,  ui32 line_no);

void set_focus_pixel(TFocusLineCollection *flc, TSysParams* sys, 
                      TTransducer* xdc, TPoint3D *points, 
                      ui32 no_points,  ui32 line_no);


void set_apodization(TApoLineCollection *alc, TSysParams* sys, 
                     TTransducer* xdc, double* times, double *apo, 
                     ui32 no_times,  ui32 line_no);



void set_filter_bank( TFocusLineCollection *flc, ui32 Nf, ui32 Ntaps, 
                                                         double *coef);
                     

#ifdef __cplusplus
  };
#endif


#endif
