#ifndef __beamform_h
#define __beamform_h
/**********************************************************************
 * NAME     : beamform.h
 * ABSTRACT : Functions for beamforing.
 **********************************************************************/
#include "types.h"
#include "focus.h"
#include "geometry.h"
#include <stdio.h>
#include <malloc.h>

#ifdef __cplusplus
  extern"C"{
#endif

typedef struct {
  TFocusLineCollection *flc;
  TApoLineCollection* alc;
  TSysParams* sys;
  double time;
  double **rf_data;
  ui32 no_samples;
  TPoint3D* elem;
  double** line;
  ui32 i;
} BFT_ThreadData;


double* beamform_apo_line_dynamic(TFocusTimeLine *ftl, TApoTimeLine* atl,
        TSysParams* sys, double time,  double **rf_data, ui32 no_samples);

double** beamform_image(TFocusLineCollection *flc, TApoLineCollection* alc,
   TSysParams* sys, double time, double **rf_data, ui32 no_samples, ui32 element_no, TPoint3D* xmt);

double* beamform_apo_line_times(TFocusTimeLine *ftl, TApoTimeLine* atl,
                            TSysParams* sys, double time, 
                            double **rf_data, ui32 no_samples);

double** apodize_fix(TApoTimeLine *atl, double **rf_data,
                                     ui32 no_samples, ui32 no_channels);

double* beamform_line_times(TFocusTimeLine *ftl, TSysParams* sys,
                        double time, double **rf_data, ui32 no_samples);
                                     

double* sum_lines_time(TFocusTimeLine *ftl, TApoTimeLine* atl, TSysParams* sys, 
                      double* rf_line1, ui32 element1, 
                      double* rf_line2, ui32 element2,
                      double time,      ui32 no_samples);

double **sum_images(TFocusLineCollection *flc, TApoLineCollection *alc,
                    TSysParams* sys,
                    double **rf1, ui32 element1,
                    double **rf2, ui32 element2, 
                    double time, ui32 no_samples);

void add_images(TFocusLineCollection *flc, TApoLineCollection *alc,
                    TSysParams* sys, double **hi_res,
                    double **lo_res, ui32 element, 
                    double time, ui32 no_samples);

void sub_images(TFocusLineCollection *flc, TApoLineCollection *alc,
                    TSysParams* sys, double **hi_res,
                    double **lo_res, ui32 element, 
                    double time, ui32 no_samples);


#ifdef __cplusplus
  };
#endif

#endif
