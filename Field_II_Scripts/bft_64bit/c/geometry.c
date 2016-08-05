/******************************************************************
 * NAME     : geometry.c
 * ABSTRACT : Geometrical functions and primitives definitions
 ******************************************************************/

#include "../h/geometry.h"

double distance( TPoint3D *p1, TPoint3D *p2)
{ double dx,dy,dz;
  dx = (p1->x - p2->x)*(p1->x - p2->x);
  dy = (p1->y - p2->y)*(p1->y - p2->y);
  dz = (p1->z - p2->z)*(p1->z - p2->z);
  return sqrt(dx + dy + dz);
}
