#ifndef __sys_params_h
  #define __sys_params_h
/*********************************************************************
 * NAME     : sys_params.h
 * ABSTRACT : Definitions of system parameters
 *
 *********************************************************************/

#include "types.h"

typedef struct sys_params{
   double  c;              /*  Speed of sound      */
   double  fs;             /*  Sampling frequency  */
}TSysParams;

#endif
