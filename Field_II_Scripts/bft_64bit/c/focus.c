/**********************************************************************
 * NAME     : focus.c
 * ABSTRACT : Managing focus time - lines.
 * 
 * CREATED  : 2000 Svetoslav Nikolov
 *
 * Ver 1.4  : Sep 2000 - Added delays using filer banks.
 * Ver 1.3  : May 2000 
 **********************************************************************/
 
#include "../h/focus.h"
#include "../h/sys_params.h"
#include "../h/geometry.h"
#include "../h/error.h"

#include <math.h>



/*********************************************************************
 * FUNCTION  : del_apodization
 * ABSTRACT  : Delete one apodization
 *********************************************************************/
void del_apodization(TApodization *a)
{
  PFUNC
  if (a->a!=NULL) free(a->a);
  a->a = NULL;
  a->time = 0.0;
}



/*********************************************************************
 * FUNCTION  : del_apo_time_line
 *********************************************************************/
void del_apo_time_line(TApoTimeLine* al)
{
  PFUNC
  if (al->a!=NULL) {
     for(;al->no_times>0;al->no_times--)
        del_apodization(al->a + al->no_times - 1);
     free(al->a);
  }
  
  al->no_times = 0;
  al->a = NULL;
}


/*********************************************************************
 * FUNCTION  : del_apo_line_collection
 *********************************************************************/
void del_apo_line_collection(TApoLineCollection* alc) 
{
   PFUNC
   if (alc->atl!=NULL){
     for(; alc->no_apo_time_lines > 0; alc->no_apo_time_lines--)
        del_apo_time_line(alc->atl + alc->no_apo_time_lines - 1);
     free(alc->atl);
   }
   alc->atl = NULL;
   alc->no_apo_time_lines = 0;
}

 
/*********************************************************************
 * FUNCTION  : del_delay
 * ABSTRACT  : Delete one delay
 *********************************************************************/
void del_delay(TDelay* d) 
{
   free(d->d); d->d = NULL;
   free(d->a); d->a = NULL;
   
}

/*********************************************************************
 * FUNCTION  : del_filter_bank
 * ABSTRACT  : delete the contents of a filter bank
 *********************************************************************/
void del_filter_bank(TFilterBank* fb)
{
  PFUNC
 /* if (fb->bank != NULL) free(fb->bank);
  if (fb->coefs != NULL) free(fb->coefs);
 */ 
  fb->Ntaps = 0;
  fb->Nf = 0;
  fb->bank = NULL;
  fb->coefs = NULL;
}



/*********************************************************************
 * FUNCTION  : del_focus_time_line
 * ABSTRACT  : delete a focus time line. 
 * ARGUMETNT : Pointer to the time line to be deleted.
 *********************************************************************/
void del_focus_time_line(TFocusTimeLine* p)
{ 
   PFUNC
   if(p->delay!=NULL)
   {
     for(; p->no_times > 0; p->no_times--)
        del_delay(p->delay + p->no_times-1);
     free(p->delay);
     p->delay = NULL;
   }
   if(p->pixels!=NULL) free(p->pixels);
}


/**********************************************************************
 *  FUNCTION  : del_focus_line_collection
 *  ABSTRACT  : 
 **********************************************************************/
void del_focus_line_collection(TFocusLineCollection* f)
{
  /*
   *   Release the memory assigned to the focus time lines
   */
  PFUNC
  if (f->ftl != NULL)
  {
    for( ;f->no_focus_time_lines>0; f->no_focus_time_lines --)
       del_focus_time_line(f->ftl + f->no_focus_time_lines - 1);
    free(f->ftl);
  }
  
  /*
   *   Release the memory assigned to the Filter bank
   */
  
  if (f->filter_bank.Nf > 0) del_filter_bank(&f->filter_bank);
  f->ftl = NULL;               /*  These 2 lines are just in case  */
  f->no_focus_time_lines = 0;  /* something somewhere went wrong   */
}


/**********************************************************************
 * FUNCTION  : set_no_lines
 **********************************************************************/
void set_no_lines(TApoLineCollection* alc, TFocusLineCollection *flc,
                                                          ui32 no_lines)
{
  PFUNC
  if(flc->no_focus_time_lines > 0) del_focus_line_collection(flc);
  flc->ftl = (TFocusTimeLine*) calloc(no_lines, sizeof(TFocusTimeLine));
  flc->no_focus_time_lines = no_lines;
  flc->use_filter_bank = 0;
  
  if(alc->no_apo_time_lines > 0) del_apo_line_collection(alc);
  alc->atl = (TApoTimeLine*) calloc(no_lines, sizeof(TApoTimeLine));
  alc->no_apo_time_lines = no_lines;
}


/*********************************************************************
 * FUNCTION  : set_center_focus
 *********************************************************************/ 
void set_center_focus(TFocusLineCollection *flc, TPoint3D *p, ui32 line_no)
{
   PFUNC
   if (line_no < flc->no_focus_time_lines){
      flc->ftl[line_no].center.x = p->x;
      flc->ftl[line_no].center.y = p->y;
      flc->ftl[line_no].center.z = p->z;
   }else{
      errprintf("%s","Error: line_no is out of range \n");
   }
}


/*********************************************************************
 * FUNCTION  : set_dynamic_focus
 *********************************************************************/ 
void set_dynamic_focus(TFocusLineCollection *flc, TTransducer* xdc,
                             ui32 line_no, double dir_xz, double dir_yz)
{ 
   PFUNC
   assert_xdc(xdc);
   if (line_no < flc->no_focus_time_lines){
      flc->ftl[line_no].dir_xz = dir_xz;
      flc->ftl[line_no].dir_yz = dir_yz;
      flc->ftl[line_no].dynamic = TRUE;
      flc->ftl[line_no].pixel = FALSE;
      flc->ftl[line_no].xdc = xdc;
   }else{
      errprintf("%s,","\"line_no\" is out of range \n");
   }
}


/*********************************************************************
 * FUNCTION : set_focus_times(flc,sys,xdc,times,delays,no_times,line_no)
 * ABSTRACT : Set the delays for focusing one line
 * ARGUMENTS: flc - Pointer to TFocusLineCollection
 *            sys - Pointer to the system parameters
 *            xdc - Pointer to transducer definition
 *            times - Pointer to an array defining the time after which
 *                    the associated delay is valid
 *            delays - array with the delays. One row per time.
 *            line_no - Number of line for which we set the delays
 *********************************************************************/
void set_focus_times(TFocusLineCollection *flc, TSysParams* sys, 
                     TTransducer* xdc, double* times, double *delays, 
                     ui32 no_times,  ui32 line_no)
{
   ui32 i, j;
   double sample_delay;
   PFUNC
   assert_xdc(xdc);
   if (line_no < flc->no_focus_time_lines){
      if(flc->ftl[line_no].no_times > 0) 
         del_focus_time_line(flc->ftl + line_no);
      flc->ftl[line_no].delay = (TDelay*)calloc(no_times + 1,sizeof(TDelay));
      
		
		flc->ftl[line_no].no_times = no_times;
      flc->ftl[line_no].dynamic = FALSE;
      flc->ftl[line_no].pixel = FALSE;
      flc->ftl[line_no].xdc = xdc;
      for (i = 0; i < no_times; i++ )
      {
         flc->ftl[line_no].delay[i].d = 
                       (si32*) malloc(xdc->no_elements*sizeof(si32));
         assert(flc->ftl[line_no].delay[i].d);
         flc->ftl[line_no].delay[i].a = 
                       (double*) malloc(xdc->no_elements*sizeof(double));
         assert(flc->ftl[line_no].delay[i].a);
         flc->ftl[line_no].delay[i].time = *times * sys->fs;
         for(j = 0; j < xdc->no_elements; j ++)
         {  
            sample_delay = *delays ++;
            sample_delay *= sys->fs;
            flc->ftl[line_no].delay[i].d[j] = (si32)floor(sample_delay);
            flc->ftl[line_no].delay[i].a[j] = ceil(sample_delay) - sample_delay;
         }
         times ++;
      }
      
      flc->ftl[line_no].delay[no_times].d = (si32*) calloc(xdc->no_elements,sizeof(si32));
      assert(flc->ftl[line_no].delay[no_times].d);
      flc->ftl[line_no].delay[no_times].a = (double*) calloc(xdc->no_elements,sizeof(double));
      assert(flc->ftl[line_no].delay[no_times].a);
      flc->ftl[line_no].delay[no_times].time = MAX_SAMPLE_NO;
   }else{
      errprintf("%s","\"line_no\" is out of range \n");
   }
}                     


/**********************************************************************
 * FUNCTION  : set_focus(flc,sys,xdc,times,points,no_times,line_no)
 * ABSTRACT  : Set the focus points.
 * ARGUMENTS : flc - Focus Line Collection
 *             sys - System parameters
 *             xdc - Transducer definition
 *             times - Time after which the associated point will be 
 *                     the next focus point.
 *             points - Focal points
 *             no_times - Number of focal zones.
 *             line_no - Number of line which we are setting.
 **********************************************************************/
void set_focus(TFocusLineCollection *flc, TSysParams* sys, 
                     TTransducer* xdc, double* times, TPoint3D *points, 
                     ui32 no_times,  ui32 line_no)
{
   ui32 i, j;
   double sample_delay;
   TPoint3D *center;

   PFUNC
   assert_xdc(xdc);
   if (line_no < flc->no_focus_time_lines){
      if(flc->ftl[line_no].no_times > 0) 
              del_focus_time_line(flc->ftl + line_no);
      flc->ftl[line_no].delay = (TDelay*)calloc(no_times+1,sizeof(TDelay));
      flc->ftl[line_no].no_times = no_times;
      flc->ftl[line_no].dynamic = FALSE;
      flc->ftl[line_no].pixel = FALSE;
      flc->ftl[line_no].xdc = xdc;
      
      center = &flc->ftl[line_no].center;
      
      for (i = 0; i < no_times; i++ )
      {
         flc->ftl[line_no].delay[i].d = (si32*) malloc(xdc->no_elements*sizeof(si32));
         assert(flc->ftl[line_no].delay[i].d);
         flc->ftl[line_no].delay[i].a = (double*) malloc(xdc->no_elements*sizeof(double));
         assert(flc->ftl[line_no].delay[i].a);

         flc->ftl[line_no].delay[i].time = *times * sys->fs;
         for(j = 0; j < xdc->no_elements; j ++)
         {  
            sample_delay = distance(center, points)*sys->fs;
            sample_delay -= distance(xdc->c+j, points)*sys->fs;
            sample_delay = sample_delay / sys->c;
             
            flc->ftl[line_no].delay[i].d[j] = (si32)floor(sample_delay);
            flc->ftl[line_no].delay[i].a[j] = sample_delay - floor(sample_delay);
         }
         
         points ++;
         times ++;
      }
      flc->ftl[line_no].delay[no_times].d = (si32*) calloc(xdc->no_elements,sizeof(si32));
      assert(flc->ftl[line_no].delay[no_times].d);
      flc->ftl[line_no].delay[no_times].a = (double*) calloc(xdc->no_elements,sizeof(double));
      assert(flc->ftl[line_no].delay[no_times].a);
      flc->ftl[line_no].delay[no_times].time = MAX_SAMPLE_NO;
   }else{
      errprintf("%s", "\"line_no\" is out of range \n");
   }
   
}

/**********************************************************************
 * FUNCTION  : set_focus_2way(flc,sys,xdc,times,points,no_times,line_no)
 * ABSTRACT  : Set the focus points. The delays are set for 2way movement
 * ARGUMENTS : flc - Focus Line Collection
 *             sys - System parameters
 *             xdc - Transducer definition
 *             times - Time after which the associated point will be 
 *                     the next focus point.
 *             points - Focal points
 *             no_times - Number of focal zones.
 *             line_no - Number of line which we are setting.
 **********************************************************************/
void set_focus_2way(TFocusLineCollection *flc, TSysParams* sys, 
                      TTransducer* xdc, double* times, TPoint3D *points, 
                      ui32 no_times,  ui32 line_no)
{
   ui32 i, j;
   double sample_delay;
   TPoint3D *center;
   
   PFUNC
   assert_xdc(xdc);
   if (line_no < flc->no_focus_time_lines){
      if(flc->ftl[line_no].no_times > 0) 
              del_focus_time_line(flc->ftl + line_no);
      flc->ftl[line_no].delay = (TDelay*)calloc(no_times+1,sizeof(TDelay));
      flc->ftl[line_no].no_times = no_times;
      flc->ftl[line_no].dynamic = FALSE;
      flc->ftl[line_no].pixel = FALSE;
      flc->ftl[line_no].xdc = xdc;
      
      center = &flc->ftl[line_no].center;
      
      for (i = 0; i < no_times; i++ )
      {
         flc->ftl[line_no].delay[i].d = (si32*) malloc(xdc->no_elements*sizeof(si32));
         assert(flc->ftl[line_no].delay[i].d);
         flc->ftl[line_no].delay[i].a = (double*) malloc(xdc->no_elements*sizeof(double));
         assert(flc->ftl[line_no].delay[i].a);

         flc->ftl[line_no].delay[i].time = *times * sys->fs;
         for(j = 0; j < xdc->no_elements; j ++)
         {  
            sample_delay = distance(center, points)*sys->fs;
            sample_delay -= distance(xdc->c+j, points)*sys->fs;
            sample_delay = 2*sample_delay / sys->c;
             
            flc->ftl[line_no].delay[i].d[j] = (si32)floor(sample_delay);
            flc->ftl[line_no].delay[i].a[j] = sample_delay - floor(sample_delay);
         }
         
         points ++;
         times ++;
      }
      flc->ftl[line_no].delay[no_times].d = (si32*) calloc(xdc->no_elements,sizeof(si32));
      assert(flc->ftl[line_no].delay[no_times].d);
      flc->ftl[line_no].delay[no_times].a = (double*) calloc(xdc->no_elements,sizeof(double));
      assert(flc->ftl[line_no].delay[no_times].a);
      flc->ftl[line_no].delay[no_times].time = MAX_SAMPLE_NO;
   }else{
      errprintf("%s", "\"line_no\" is out of range \n");
   }
   
}


/**********************************************************************
 * FUNCTION  : set_focus_pixel(flc,sys,xdc, points, no_times,line_no)
 * ABSTRACT  : Set the focus points. The focal points correspond to 
 *             irregular grid on the screen. Therefore no starting time
 *             is set. The offset of the signals is calculated 
 * ARGUMENTS : flc - Focus Line Collection
 *             sys - System parameters
 *             xdc - Transducer definition
 *             points - Focal points
 *             no_times - Number of focal zones.
 *             line_no - Number of line which we are setting.
 **********************************************************************/
void set_focus_pixel(TFocusLineCollection *flc, TSysParams* sys, 
                      TTransducer* xdc, TPoint3D *points, 
                      ui32 no_times,  ui32 line_no)
{
   ui32 i;
   
   PFUNC
   assert_xdc(xdc);
   if (line_no < flc->no_focus_time_lines){
      if(flc->ftl[line_no].no_times > 0) 
              del_focus_time_line(flc->ftl + line_no);
              
      flc->ftl[line_no].no_times = no_times;
      flc->ftl[line_no].dynamic = FALSE;
      flc->ftl[line_no].pixel = TRUE;
      flc->ftl[line_no].xdc = xdc;
      flc->ftl[line_no].pixels = (TPoint3D *)malloc(no_times * sizeof(TPoint3D));
      if (flc->ftl[line_no].pixels == NULL){
         errprintf("%s", "Cannot allocate memory \n");
         assert(flc->ftl[line_no].pixels);
      }
      for (i = 0; i < no_times; i++ )
         flc->ftl[line_no].pixels[i] = *points ++;

   }else{
      errprintf("%s", "\"line_no\" is out of range \n");
   }
}


/*********************************************************************
 * FUNCTION : set_apodization(alc,sys,xdc,times,apo ,no_times,line_no)
 * ABSTRACT : Set the delays for focusing one line
 * ARGUMENTS: alc - Pointer to TApoLineCollection
 *            sys - Pointer to the system parameters
 *            xdc - Pointer to transducer definition
 *            times - Pointer to an array defining the time after which
 *                    the associated apodization is valid
 *            apo  - array with the apodizations. One row per time.
 *            line_no - Number of line for which we set the delays
 *********************************************************************/
void set_apodization(TApoLineCollection *alc, TSysParams* sys, 
                     TTransducer* xdc, double* times, double *apo, 
                     ui32 no_times,  ui32 line_no)
{
  ui32 i, j;
  int different;
  int zero_times;
  int to_delete;
		
  PFUNC
    assert_xdc(xdc);
	
  if (line_no < alc->no_apo_time_lines){
    zero_times = alc->atl[line_no].no_times == 0;
    different = alc->atl[line_no].no_times != no_times;
    to_delete = !zero_times && different;
    if(!zero_times && different) del_apo_time_line(alc->atl + line_no);
		
    if (zero_times || to_delete){
      alc->atl[line_no].a = (TApodization*)calloc(no_times + 1,sizeof(TApodization));
      alc->atl[line_no].no_times = no_times;
      for (i = 0; i < no_times; i++ )
      	{
	  alc->atl[line_no].a[i].a = 
	    (double*) malloc(xdc->no_elements*sizeof(double));
	  assert(alc->atl[line_no].a[i].a);
         
	  alc->atl[line_no].a[i].time = *times * sys->fs;
         
	  for(j = 0; j < xdc->no_elements; j ++)
	    alc->atl[line_no].a[i].a[j] = *apo++;
	  times ++;
	}
    }else{
      for (i = 0; i < no_times; i++ )
      	{
	  alc->atl[line_no].a[i].time = *times * sys->fs;
	  for(j = 0; j < xdc->no_elements; j ++)
	    alc->atl[line_no].a[i].a[j] = *apo++;
	  times ++;
	}
    }
    alc->atl[line_no].a[no_times].a= (double*) calloc(xdc->no_elements,sizeof(double));
    assert(alc->atl[line_no].a[no_times].a);
    alc->atl[line_no].a[no_times].time = MAX_SAMPLE_NO;
  }else{
    errprintf("%s", "\"line_no\" is out of range \n");
  }
}                     



/**********************************************************************
 * FUNCTION : set filter bank
 * ABSTRACT : Setting the filter bank. Filter banks are used for 
 *            phase shift in the beamformation process. The process can
 *            be viewed as resampling the signal at a higher frequency
 *            by zero-padding, then shifting with a number of samples
 *            and then resampling back to the original frequency. 
 *            
 *            Because most of the samples in the up-sampling are == 0
 *            the process can be viewed as applying a simple filter 
 *            on the original signal. For each of the different delays
 *            a different filter is applied, thus a "bank of filters"
 *            is needed.
 *
 * ARGUMENTS: Nf    - Resolution of the fine delay. 
 *            Ntaps - Number of filter coefficients used per per delay.
 *                    The total number off coefficients at the highest 
 *                    frequency must be Nf*Ntaps - 1
 *            coefs - The actual coefficients in the filter.
 **********************************************************************/
void set_filter_bank( TFocusLineCollection *flc, ui32 Nf, ui32 Ntaps,
                                                         double *coef)
{
  int i,j;
  
  PFUNC  

  if (flc->use_filter_bank){
    del_filter_bank(&flc->filter_bank);
    flc->use_filter_bank = 0;
  }

  flc->filter_bank.Nf = Nf;
  flc->filter_bank.Ntaps = Ntaps;
  
  flc->filter_bank.coefs = calloc(Nf*Ntaps,sizeof(double));
  if (flc->filter_bank.coefs == NULL){
     errprintf("%s","Cannot allocate memory for the filter coefficients\n");
     assertp(flc->filter_bank.coefs);
  }
  
  flc->filter_bank.bank = calloc(Nf, sizeof(double*));
  if (flc->filter_bank.bank == NULL){
     free(flc->filter_bank.coefs); flc->filter_bank.coefs = NULL;
     errprintf("%s","Cannot allocate memory for the filter banks. \n");
     assertp(flc->filter_bank.bank);
  }
  
  for (i = Nf-1; i>-1; i--){
     flc->filter_bank.bank[i] = flc->filter_bank.coefs + i*Ntaps;
  }
  
  
  for (i = 0; i < Nf; i++ )
     for (j = 0; j < Ntaps; j++)
       flc->filter_bank.bank[i][j] = coef[i + j*Nf];
  

  flc->use_filter_bank = 1;
}

