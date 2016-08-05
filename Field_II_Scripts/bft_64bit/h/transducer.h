#ifndef __transducer_h
  #define __transducer_h


#include "types.h"

typedef struct transducer{
   ui16 no_elements;              /* Number of elements             */
   TPoint3D* c;                   /*  Center of the transducer      */
   struct transducer *next;
}TTransducer;


#ifdef __cplusplus
  extern"C"{
#endif

TTransducer* bft_transducer(ui32 no_elements, TPoint3D *p);
void bft_transducer_set(TTransducer* x, ui32 no_elements, TPoint3D *p);
void bft_free_xdc(TTransducer* x);
void bft_free_all_xdc();
si32 is_xdc_valid(TTransducer* x);
void assert_xdc(TTransducer *x);


#ifdef __cplusplus
  };
#endif


#endif
