#ifndef __mex_beamform_h
#define __mex_beamform_h

/*********************************************************************
 * Header file for the beamformer
 *
 *********************************************************************/
#define V4_COMPAT
#include "mex.h"
#include "beamform.h"
#include "focus.h"
#include "transducer.h"

#include <math.h>

#define BFT_INIT             0
#define BFT_END              1
#define BFT_PARAM            2
#define BFT_NO_LINES         3
#define BFT_TRANSDUCER       4
#define BFT_XDC_FREE         5
#define BFT_CENTER_FOCUS     6
#define BFT_FOCUS            7
#define BFT_FOCUS_TIMES      8
#define BFT_APODIZATION      9
#define BFT_DYNAMIC_FOCUS    10
#define BFT_BEAMFORM         11
#define BFT_SUM_IMAGES       12
#define BFT_ADD_IMAGES       13
#define BFT_SUM_APODIZATION  14
#define BFT_SUB_IMAGES       15
#define BFT_FOCUS_2WAY       16
#define BFT_FOCUS_PIXEL      17
#define BFT_FILTER           18
#define BFT_DELAY            19
#define BFT_DELAY_FILTER     20
#define BFT_XDC_SET          21

#endif
