/*********************************************************************
 * NAME      : transducer.c
 * ABSTRACT  : Manimpulating transducers.
 *
 *********************************************************************/
 
#include "../h/transducer.h"
#include <stdlib.h>

static TTransducer *xdc = NULL;   /* Pointer to the transducer 
                                   * definitions 
                                   */


/*********************************************************************
 *  bft_transducer  : Add a new transducer definition to the chain
 *     of transducer definitions.
 *********************************************************************/
TTransducer* bft_transducer(ui32 no_elements, TPoint3D *p)
{
   TTransducer *x = NULL;
   ui32 i;
   
   x = xdc;
   xdc = malloc(sizeof(TTransducer));
   assert(xdc!=NULL);
   xdc->next = x;
   
   xdc->no_elements = no_elements;
   xdc->c = (TPoint3D*)malloc(no_elements * sizeof(TPoint3D));
   
   for (i = 0; i < no_elements; i++) {
       xdc->c[i].x = p[i].x;
       xdc->c[i].y = p[i].y;
       xdc->c[i].z = p[i].z;
   }
   
   return xdc;     
}


/*********************************************************************
 *  bft_transducer_set  : Set new coordinates for the transducer elements
 *
 *********************************************************************/
void bft_transducer_set(TTransducer* x, ui32 no_elements, TPoint3D *p)
{
   
   ui32 i;
	if (is_xdc_valid(x)){   
		if (x->no_elements == no_elements){
		   for (i = 0; i < no_elements; i++) {
      		x->c[i].x = p[i].x;
       		x->c[i].y = p[i].y;
	       	x->c[i].z = p[i].z;
   		}
		}
   }
}


/*********************************************************************
 *  bft_free_xdc  - Free a transducer definition. 
 *********************************************************************/
void bft_free_xdc(TTransducer* x)
{
   TTransducer *c, *p=NULL;
   
   if (xdc == NULL) return;
#ifdef DEBUG
  printf("Freeing at address %x \n", (ui32)x);
#endif   
   if (x == xdc) {
      xdc = xdc->next;
      free(x->c);
      free(x);
   } else {
      c = xdc;
      while ((c != NULL) && (c!=x)){ p = c; c = c->next;}
      if (c==x){
         c = x->next;
         p->next = c;
         free(x->c);
         free(x);
      }else{
         printf("bft_free_xdc: Cannot find object to free \n");
      }
   }
  
}


/*********************************************************************
 *  bft_free_all_xdc - Free all transducer definitions
 *********************************************************************/
 
void bft_free_all_xdc()
{
  while(xdc!=NULL) bft_free_xdc(xdc);
}



/**********************************************************************
 * FUNCTION  : is_xdc_valid
 * ABSTRACT  : Check if a pointer points to a transducer definition in 
 *             in the transducers chain
 **********************************************************************/
si32 is_xdc_valid(TTransducer* x)
{
  TTransducer *c;

  if (xdc!= NULL && x!= NULL){
     c = xdc;
     while ((c != NULL) && (c!=x)){c = c->next;}
     if (c == xdc) return TRUE;
   }
   return FALSE;
}


/*********************************************************************
 * FUNCTION : assert_xdc
 * ABSTRACT : The function checks for valid transducer definition, and
 *            if there isn't one, makes the program abort.
 *********************************************************************/
void assert_xdc(TTransducer *x)
{
  if (!is_xdc_valid(x)) abort();
}

