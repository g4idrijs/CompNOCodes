#ifndef __geometry_h
/*********************************************************************
 * NAME     : geometry.h
 * ABSTRACT : Geometrical functions and primitives definitions
 *
 *********************************************************************/

#include "types.h"
#include <math.h>

#ifdef __cplusplus
  extern"C"{
#endif

double distance( TPoint3D *p1, TPoint3D *p2);

#ifdef __cplusplus
  };
#endif

#endif
