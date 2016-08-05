/*********************************************************************
 * NAME     : mex_beamform
 * ABSTRACT : Matlab interface to the beamforming toolbox
 * CREATED  : Feb 2000, Svetoslav Nikolov
 *********************************************************************/
 
#include "../h/mex_beamform.h" 
#include "../h/error.h"
#include "../h/motion.h"
#include <signal.h>
#include <string.h>

#ifdef SPECIAL_CASE
	#include <unistd.h>
#endif

static int mexNoEntries = 0;
static TSysParams sys;
static TFocusLineCollection *flc;
static TApoLineCollection *alc;
static TApoLineCollection *salc;   /* Sum apo-line collection*/

static int initialized = FALSE;


/******************************************************************
 * FUNCTION : mexBFTExit - This is the function taking care of the 
 *            clean-up process
 ******************************************************************/
void mexBFTExit(void)
{
   if (mexNoEntries == 0) return;

   mexNoEntries = 0;
   if (initialized == FALSE) return;

#ifdef  MALLOC_CHECK_
  printf("MALLOC_CHECK_ is %d \n", MALLOC_CHECK_);
#endif
   if (flc != NULL){
#ifdef DEBUG   
      printf("Freeing Focusing settings \n");
#endif      
      del_focus_line_collection(flc);
      free(flc); flc = NULL;
   }
   
   if (alc != NULL){
#ifdef DEBUG   
      printf("Freeing all apodization settings \n");
#endif      
      del_apo_line_collection(alc);
      free(alc); alc = NULL;
   }

   if (salc != NULL){
#ifdef DEBUG   
      printf("Freeing all summation apodization settings \n");
#endif      
      del_apo_line_collection(salc);
      free(alc); salc = NULL;
   }
   
#ifdef DEBUG   
   printf("Freeing all transducers \n");
#endif   
   bft_free_all_xdc();
   initialized = FALSE;
#ifdef SPECIAL_CASE
   nice(0);
#endif 
   printf("\t**************************************************\n");
   printf("\t       Exiting the BeamForming Toolbox\n");
   printf("\t**************************************************\n");
}


/*********************************************************************
 * FUNCTION : atAbort
 * ABSTRACT : The function is called when a signal "SIGABRT" is caught
 *            This signal is sent by the function abort(), which on its
 *            order is called by assert(). This is one of the simplest
 *            mechanisms for error trapping in case of memory allocation
 *            problems.
 *              atAbort() makes Matlab release the memory allocated by
 *            the toolbox.
 *********************************************************************/
 
void atAbort(int sig_no)
{ 
   printf("\t***********************************************************\n");
   printf("\t*                                                         *\n");
   printf("\t* BeamForming Toolbox  message :                          *\n");
   printf("\t*                                                         *\n");
   printf("\t* Signal 'abort' caught                                   *\n");
   printf("\t* One of the probable causes is bad memory pointer        *\n");
   printf("\t* Check if you have supplied wrong transducer definition  *\n");
   printf("\t*                                                         *\n");
   printf("\t***********************************************************\n");
   mexErrMsgTxt("\n");
}

/*********************************************************************
 * FUNCTION : bft_init
 * ABSTRACT : Initialization of the beamforming toolbox
 *
 *********************************************************************/
void bft_init(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{ 
  printf("\t**************************************************************\n");
  printf("\t*                                                            *\n");
  printf("\t*               Beamforming  Toolbox                         *\n");
  printf("\t*                       by                                   *\n");
  printf("\t*                Svetoslav Nikolov                           *\n");
  printf("\t*                                                            *\n");
  printf("\t*              Version 1.4,  May 12, 2003                    *\n");
  printf("\t*                                                            *\n");
  printf("\t**************************************************************\n");
  
  /*
   *   Initialize the physical constants
   */
  sys.fs = 40e6;
  sys.c = 1540.0;
  
  /*
   *  Allocate the memory, necessary for the beamforming and apodization
   *  data. Allocate memory for at least one line
   */
  flc = (TFocusLineCollection *) calloc(1, sizeof(TFocusLineCollection));
  assert(flc!= NULL);
  alc = (TApoLineCollection*) calloc(1, sizeof(TApoLineCollection));
  assert(alc != NULL);
  set_no_lines(alc, flc, 1);

  salc = (TApoLineCollection*) calloc(1, sizeof(TApoLineCollection));
  assert(salc != NULL);
  set_no_lines(salc, flc, 1);
  flc->use_filter_bank = 0;
  initialized = TRUE;  
}

/*********************************************************************
 * FUNCTION : bft_end
 * ABSTRACT : Initialization of the beamforming toolbox
 *
 *********************************************************************/
void bft_end(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   mexBFTExit();
   initialized = 0;
}

/*********************************************************************
 * FUNCTION : bft_param
 * ABSTRACT : Set one system parameter.
 *********************************************************************/
void bft_param(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   ui32 len;
   char param_name[80];
   
   if(!initialized)
      mexErrMsgTxt("\nToolbox is not initialized.\n");
   
   if(nrhs!=3)
      mexErrMsgTxt("\nExpecting 'name_of_param' and 'value_of_param'\n");

   if (!mxIsChar(prhs[1]))
      mexErrMsgTxt("\n'name_of_param' must be string \n");

   if (!mxIsDouble(prhs[2])) 
      mexErrMsgTxt("\n'value_of_param' must be of type 'double'\n");

   if ((len = mxGetN(prhs[1]))>79)
       mexErrMsgTxt("\nThe string is too long.\n");
   
   if(mxGetString(prhs[1], param_name, len+1))
        mexErrMsgTxt("\nBad string argument\n");
   
   if(mxGetN(prhs[2])>1 || mxGetM(prhs[2])>1)
       mexErrMsgTxt("\nThe value of the parameter must be scalar\n");
   
   if (!strcmp(param_name,"c")){
      sys.c = mxGetScalar(prhs[2]);
   }else if(!strcmp(param_name,"fs")){
      sys.fs = mxGetScalar(prhs[2]);
   }else{
      printf("\nUnknown parameter name '%s'\n ",param_name);
      mexErrMsgTxt("");
   }
   
}


/*********************************************************************
 * FUNCTION : bft_no_lines
 * ABSTRACT : Allocate memory for given amount of focusing lines.
 *********************************************************************/
void bft_no_lines(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   ui32 no_lines;
   if(!initialized)
      mexErrMsgTxt("\nToolbox is not initialized.\n");
   
   if(nrhs!=2)
      mexErrMsgTxt("\nExpecting 'no_of_lines' and 'value_of_param'\n");

   if (mxGetN(prhs[1])>1 || mxGetM(prhs[1])>1)
      mexErrMsgTxt("\n'no_of_lines' must be scalar \n");
   no_lines = (ui32)(floor)(mxGetScalar(prhs[1]) + 0.5);
   set_no_lines(alc, flc, no_lines);
   set_no_lines(salc, flc, no_lines);
}


/*********************************************************************
 *  bft_transducer - Define a transducer.
 *********************************************************************/
void call_bft_transducer(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  TPoint3D *centers;
  ui32 no_elements;
  double* ret_val;
  int dim[2];
  
  if (!initialized)
     mexErrMsgTxt("\nToolbox is not initialized\n");
  
  if (nlhs!=1)   
     mexErrMsgTxt("\nMust return a handle \n");
     
  if (nrhs!=2)
     mexErrMsgTxt("\nExpecting an array with the centers of the elements\n");
  
  if (mxGetM(prhs[1])!=3)
     mexErrMsgTxt("The points must be defined with their [x y z] coordinates\n");
     
  no_elements = mxGetN(prhs[1]);
  centers = (TPoint3D*)mxGetPr(prhs[1]);

  xdc = bft_transducer(no_elements, centers);
  
  /*
   *     Return the values back to the main program
   */  
   
  dim[0] = 1;  dim[1] = 1;
  plhs[0] = mxCreateNumericArray(2, dim, mxDOUBLE_CLASS, mxREAL);
  assert(plhs[0]!=NULL);
  ret_val = mxGetPr(plhs[0]);
  *ret_val = (uint64)xdc;
}



/******************************************************************
 * FUNCTION  : bft_xdc_free
 * ABSTRACT  : Free the memory of the transducer
 *
 ******************************************************************/

void call_bft_free_xdc( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{ 
  TTransducer *xdc;
  uint64 address;
  
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=2)
     mexErrMsgTxt("\nExpecting a pointer to the transducer\n");
     
  if (mxGetN(prhs[1]) > 1 || mxGetN(prhs[1])> 1)     
     mexErrMsgTxt("\nExpecting scalar for a tranducer handle \n");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  bft_free_xdc(xdc);
}


/******************************************************************
 * FUNCTION  : bft_center_focus
 * ABSTRACT  : Set the reference point for the for the focus 
 *             calculations
 ******************************************************************/
void bft_center_focus(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   TPoint3D *p;
   ui32 line_no;
   
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=3)
     mexErrMsgTxt("\nExpecting arrays of points and line number\n");
  
  if (mxGetM(prhs[1])!=3)
     mexErrMsgTxt("\nThere must be 3 coordinates (x,y,z) per point \n");
 
  if (mxGetM(prhs[2])>1 || mxGetN(prhs[2])>1)
     mexErrMsgTxt("\nExpecting scalar for number of line \n");
  
  p = (TPoint3D*)mxGetPr(prhs[1]);
  line_no = (ui32)floor(mxGetScalar(prhs[2])) - 1;
  set_center_focus(flc, p, line_no);
}


/*******************************************************************
 * FUNCTION : bft_focus
 * ABSTRACT : 
 *******************************************************************/ 
void bft_focus(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  TPoint3D *p;
  double *times;
  uint64 address;
  ui32 no_times;
  ui32 line_no;
  ui32 m,n;
  
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=5){
      printf("\nExpecting a pointer to aperture, 'times'");
      printf(" and 'point' vectors, and 'line_no'\n");
      mexErrMsgTxt("");
   } 
  
  
  /*  Get 'xdc'  */
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1)
     mexErrMsgTxt("The pointer to the transducer must be only one");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  
  
  /* Get 'times' and 'no_times' */
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m > 1 && n > 1)
     mexErrMsgTxt("'times' must be vector, not a mxArray");
  no_times = m*n;
  times = mxGetPr(prhs[2]);
    
  
  /*Get 'points' */
  m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
  if (m!=3 || n!=no_times){
     printf("The length of 'times' is %d \n", no_times);
     printf("The size of the 'point' array must be 3 by %d \n", no_times);
     mexErrMsgTxt("");
  }
  
  p = (TPoint3D*) mxGetPr(prhs[3]);
  
  /* Get 'line_no' */
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1){
     mexErrMsgTxt("'line_no' must be a single number\n" );
  }
  line_no = (ui32)((ui32)floor(mxGetScalar(prhs[4])) & 0xffff) - 1;
  set_focus(flc, &sys, xdc, times, p, no_times, line_no);
}


/*******************************************************************
 * FUNCTION : bft_focus_pixel
 * ABSTRACT : Set the focusing based on pixels
 *******************************************************************/ 
void bft_focus_pixel(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  TPoint3D *p;
  uint64 address;
  ui32 no_times;
  ui32 line_no;
  ui32 m;
  
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=4){
      printf("\nExpecting a pointer to aperture, ");
      printf(" 'point' vectors, and 'line_no'\n");
      mexErrMsgTxt("");
   } 
  
  
  /*  Get 'xdc'  */
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1)
     mexErrMsgTxt("The pointer to the transducer must be only one");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  
 
  /*Get 'points' */
  
  m = mxGetM(prhs[2]); no_times = mxGetN(prhs[2]);
  if (m!=3){
     printf("The size of the 'point' array must be 3 by number of points \n");
     mexErrMsgTxt("");
  }
  
  p = (TPoint3D*) mxGetPr(prhs[2]);
  
  /* Get 'line_no' */
  if (mxGetM(prhs[3])>1 || mxGetN(prhs[3])>1){
     mexErrMsgTxt("'line_no' must be a single number\n" );
  }
  line_no = (ui32)((ui32)floor(mxGetScalar(prhs[3])) & 0xffff) - 1;
  set_focus_pixel(flc, &sys, xdc, p, no_times, line_no);
}



/*******************************************************************
 * FUNCTION : bft_focus_2way
 * ABSTRACT : 
 *******************************************************************/ 
void bft_focus_2way(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  TPoint3D *p;
  double *times;
  uint64 address;
  ui32 no_times;
  ui32 line_no;
  ui32 m,n;
  
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=5){
      printf("\nExpecting a pointer to aperture, 'times'");
      printf(" and 'point' vectors, and 'line_no'\n");
      mexErrMsgTxt("");
   } 
  
  
  /*  Get 'xdc'  */
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1)
     mexErrMsgTxt("The pointer to the transducer must be only one");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  
  
  /* Get 'times' and 'no_times' */
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m > 1 && n > 1)
     mexErrMsgTxt("'times' must be vector, not a mxArray");
  no_times = m*n;
  times = mxGetPr(prhs[2]);
    
  
  /*Get 'points' */
  m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
  if (m!=3 || n!=no_times){
     printf("The length of 'times' is %d \n", no_times);
     printf("The size of the 'point' array must be 3 by %d \n", no_times);
     mexErrMsgTxt("");
  }
  
  p = (TPoint3D*) mxGetPr(prhs[3]);
  
  /* Get 'line_no' */
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1){
     mexErrMsgTxt("'line_no' must be a single number\n" );
  }
  line_no = (ui32)((ui32)floor(mxGetScalar(prhs[4])) & 0xffff) - 1;
  set_focus_2way(flc, &sys, xdc, times, p, no_times, line_no);
}

/********************************************************************
 * FUNCTION  : bft_focus_times
 * ABSTRCT   : Set a focus time - line
 ********************************************************************/
void bft_focus_times(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  double *delays;
  double *times;
  uint64 address;
  ui32 no_times;
  ui32 line_no;
  ui32 m,n;
   
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=5){
      printf("\nExpecting a pointer to aperture, 'times'");
      printf(" and 'point' vectors, and 'line_no'\n");
      mexErrMsgTxt("");
   } 
  
  
  /*  Get 'xdc'  */
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1)
     mexErrMsgTxt("The pointer to the transducer must be only one");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  
  
  /* Get 'times' and 'no_times' */
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m > 1 && n > 1)
     mexErrMsgTxt("'times' must be vector, not a mxArray");
  no_times = m*n;
  times = mxGetPr(prhs[2]);
    
  
  /*Get 'points' */
  m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
  if (m!=xdc->no_elements || n!=no_times){
     printf("The length of 'times' is %d \n", no_times);
     printf("The number of elements in the transducer is %d \n",xdc->no_elements);
     printf("The size of the 'delays' mxArray must be %d by %d \n",xdc->no_elements, no_times);
     mexErrMsgTxt("");
  }
  delays = mxGetPr(prhs[3]);
  
  /* Get 'line_no' */
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1){
     mexErrMsgTxt("'line_no' must be a single number\n" );
  }
  line_no = (ui32)((ui32)floor(mxGetScalar(prhs[4])) & 0xffff) - 1;
  set_focus_times(flc, &sys, xdc, times, delays, no_times, line_no);
   
   
}
 

/*******************************************************************
 * FUNCTION : bft_apodization
 * ABSTRACT : Set the apodization time line 
 *******************************************************************/ 
void bft_apodization(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  double *apodization;
  double *times;
  uint64 address;
  ui32 no_times;
  ui32 line_no;
  ui32 m,n;
   
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=5){
      printf("\nExpecting a pointer to aperture, 'times'");
      printf(" 'apodization' mxArray and 'line_no'\n");
      mexErrMsgTxt("");
   } 
  
  
  /*  Get 'xdc'  */
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1)
     mexErrMsgTxt("The pointer to the transducer must be only one");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  
  
  /* Get 'times' and 'no_times' */
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m > 1 && n > 1)
     mexErrMsgTxt("'times' must be vector, not a mxArray");
  no_times = m*n;
  times = mxGetPr(prhs[2]);
    
  
  /*Get 'points' */
  m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
  if (m!=xdc->no_elements || n!=no_times){
     printf("The length of 'times' is %d \n", no_times);
     printf("The number of elements in the transducer is %d \n",xdc->no_elements);
     printf("The size of the 'apodization' mxArray must be %d by %d \n",xdc->no_elements, no_times);
     mexErrMsgTxt("");
  }
  apodization = mxGetPr(prhs[3]);
  
  /* Get 'line_no' */
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1){
     mexErrMsgTxt("'line_no' must be a single number\n" );
  }
  line_no = (ui32)((ui32)floor(mxGetScalar(prhs[4])) & 0xffff) - 1;

  set_apodization(alc, &sys, xdc, times,apodization, no_times, line_no);

}




/*******************************************************************
 * FUNCTION : bft_sum_apodization
 * ABSTRACT : Set the apodization time line 
 *******************************************************************/ 
void bft_sum_apodization(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  double *apodization;
  double *times;
  uint64 address;
  ui32 no_times;
  ui32 line_no;
  ui32 m,n;
   
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=5){
      printf("\nExpecting a pointer to aperture, 'times'");
      printf(" 'apodization' mxArray and 'line_no'\n");
      mexErrMsgTxt("");
   } 
  
  
  /*  Get 'xdc'  */
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1)
     mexErrMsgTxt("The pointer to the transducer must be only one");
  
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*)address;
  
  
  /* Get 'times' and 'no_times' */
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m > 1 && n > 1)
     mexErrMsgTxt("'times' must be vector, not a mxArray");
  no_times = m*n;
  times = mxGetPr(prhs[2]);
    
  
  /*Get 'points' */
  m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
  if (m!=xdc->no_elements || n!=no_times){
     printf("The length of 'times' is %d \n", no_times);
     printf("The number of elements in the transducer is %d \n",xdc->no_elements);
     printf("The size of the 'apodization' mxArray must be %d by %d \n",xdc->no_elements, no_times);
     mexErrMsgTxt("");
  }
  apodization = mxGetPr(prhs[3]);
  
  /* Get 'line_no' */
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1){
     mexErrMsgTxt("'line_no' must be a single number\n" );
  }
  line_no = (ui32)((ui32)floor(mxGetScalar(prhs[4])) & 0xffff) - 1;

  set_apodization(salc, &sys, xdc, times, apodization, no_times, line_no);
}

/*******************************************************************
 * FUNCTION : bft_dynamic_focus
 * ABSTRACT : This function should implement dynamic focusing
 *******************************************************************/
void bft_dynamic_focus(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  TTransducer *xdc;
  uint64 address;
  double dir_xz;
  double dir_yz;
  ui32 line_no;
  
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=5){
      printf("\nExpecting a pointer to aperture, 'times'");
      printf(" and 'dir_zx', and 'dir_zy'\n");
      mexErrMsgTxt("");
  } 
  
  if (mxGetM(prhs[1])>1 || mxGetN(prhs[1])>1){
     mexErrMsgTxt("The pointer must be a single value \n");
  }
  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*) address;
  
  if (mxGetM(prhs[2])>1 || mxGetN(prhs[2])>1){
     mexErrMsgTxt("'dir_xz' must be a single value \n");
  }
  dir_xz = mxGetScalar(prhs[2]);
  

  if (mxGetM(prhs[3])>1 || mxGetN(prhs[3])>1){
     mexErrMsgTxt("'dir_yz' must be a single value \n");
  }
  dir_yz = mxGetScalar(prhs[3]);
  
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1){
     mexErrMsgTxt("'line_no' must be a single value");
  }
  line_no = (ui32)floor(mxGetScalar(prhs[4]))-1; 
  
  set_dynamic_focus(flc, xdc, line_no,  dir_xz, dir_yz);

}

/*******************************************************************
 * FUNCTION : bft_beamform
 * ABSTRACT : Beamform the image.
 *******************************************************************/
void bft_beamform(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   double Time;        /* Starting time of the first sample              */
   ui32 no_samples;    /* Number of samples per RF line                  */
   ui32 no_elements;   /* Number of elements that have recorded this line*/
   ui32 element_no=-1;  /* No of element, which is used in transmit       */
   double *ptr;        /* Pointer to the array passed by Matlab          */
   double **rf_data;   /* 2D array, passed to the beamforming  routine   */
   double **bf_data;   /* 2D array with the beamformed data              */  
	TPoint3D *xmt=NULL;
   ui32 i; 
   
    
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
   
  if (nrhs!=3 && nrhs!=4)
      mexErrMsgTxt("\nExpecting  'time', 'rf_data' and (optionally) 'element_no'\n");
  
  if (mxGetM(prhs[1])> 1 || mxGetN(prhs[1])>1)
      mexErrMsgTxt("\nExpecting a single value for 'time' \n");

  if (nrhs==4){
    if((mxGetM(prhs[3]) * mxGetN(prhs[3]))==1)
    	element_no = (ui32)floor(mxGetScalar(prhs[3])) - 1;
	 else if((mxGetM(prhs[3]) * mxGetN(prhs[3]))==3)
	 	xmt = (TPoint3D*)mxGetPr(prhs[3]);
	 else
	 	mexErrMsgTxt("The transmitting aperture must be given either as coordinates or as an index\n");
	 
  }
  Time = mxGetScalar(prhs[1]);
  
  no_elements = mxGetN(prhs[2]);
  no_samples = mxGetM(prhs[2]);
  ptr = mxGetPr(prhs[2]);
  
  rf_data = (double**)calloc(no_elements, sizeof(double*));
  if (rf_data == NULL)
     mexErrMsgTxt("Cannot allocate memory \n");
  
  for (i = 0; i < no_elements; i++)
     rf_data[i] = ptr + i*no_samples;
  
  bf_data = beamform_image(flc, alc, &sys, Time, rf_data, no_samples, element_no, xmt);
  
  if (bf_data == NULL)
     mexErrMsgTxt("Beamforming is unsuccessful \n");
  
  free(rf_data);
  if ((flc->no_focus_time_lines == 1) && (flc->ftl[0].pixel == TRUE)){
     plhs[0] = mxCreateDoubleMatrix(flc->ftl[0].no_times,1,mxREAL);
     no_samples = flc->ftl[0].no_times;
  }else
     plhs[0] = mxCreateDoubleMatrix(no_samples,flc->no_focus_time_lines,mxREAL);

  ptr = mxGetPr(plhs[0]);
  
  for (i = 0; i<flc->no_focus_time_lines; i++){
     memcpy(ptr, bf_data[i], no_samples*sizeof(double));
     ptr += no_samples;
     free(bf_data[i]);
  }
  free(bf_data);
   
}

/*******************************************************************
 * FUNCTION : bft_sum_images
 * ABSTRACT : Sum 2 low resolution images in one high resoltuion.
 *******************************************************************/ 
void bft_sum_images(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  ui32 element1;
  ui32 element2;
  ui32 no_samples;
  ui32 m,n,i;
  double **rf1;
  double **rf2;
  double ** hi_res;
  double *ptr1, *ptr2;
  double time;

  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
  
  if (nrhs != 6)
     mexErrMsgTxt("Expecting 'image1', 'ele1', 'image2', 'ele2',and 'time'\n");
   
  m = mxGetM(prhs[1]); n= mxGetN(prhs[1]);
  no_samples = m; 
  if ( n!= flc->no_focus_time_lines )
     mexErrMsgTxt("The number of columns of 'image1' must be equal to the number of lines.\n");
       
  
  if(mxGetM(prhs[2])>1 || mxGetN(prhs[2])>1)
     mexErrMsgTxt(" Expecting a scalar for 'ele1'");
  
  m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
  if (m!= no_samples) 
      mexErrMsgTxt("The two images must have the same number of samples per line\n");
 
  if ( n!= flc->no_focus_time_lines )
     mexErrMsgTxt("The number of columns of 'image2' must be equal to the number of lines.\n");
  
  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1)
     mexErrMsgTxt("Expecting scalar for 'ele2'");

  if (mxGetM(prhs[5])>1 || mxGetN(prhs[5])> 1)
     mexErrMsgTxt("Expecting a scalar for 'time'\n");
  
  element1 = (ui32)floor(mxGetScalar(prhs[2]));
  element2 = (ui32)floor(mxGetScalar(prhs[4]));
  time = mxGetScalar(prhs[5]);
  
  rf1 = (double**)malloc(flc->no_focus_time_lines*sizeof(double*));
  assert(rf1);
  
  rf2 = (double**)malloc(flc->no_focus_time_lines * sizeof(double*));
  if (rf2 == NULL) {free(rf1); abort();}
  
  ptr1 = mxGetPr(prhs[1]);
  ptr2 = mxGetPr(prhs[3]);
  
  for (i = 0; i < flc->no_focus_time_lines; i++){
     rf1[i] = ptr1 + no_samples*i;
     rf2[i] = ptr2 + no_samples*i;
  }
  
  
  hi_res = sum_images(flc, alc, &sys,rf1, element1, rf2, element2,time, no_samples );

  plhs[0] = mxCreateDoubleMatrix(no_samples,flc->no_focus_time_lines,mxREAL);
  ptr1 = mxGetPr(plhs[0]);
  
  for (i = 0; i<flc->no_focus_time_lines; i++){
     memcpy(ptr1, hi_res[i], no_samples*sizeof(double));
     ptr1 += no_samples;
     free(hi_res[i]);
  }
  free(hi_res);
  free(rf1);
  free(rf2);
}



/*******************************************************************
 * FUNCTION : bft_add_images
 * ABSTRACT : Sum 2 low resolution images in one high resoltuion.
 *******************************************************************/ 
void bft_add_images(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  ui32 element;
  ui32 no_samples;
  ui32 m,n,i;
  double **lo_res;
  double **hi_res;
  double *ptr1, *ptr2;
  double time;

  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
  
  if (nrhs != 5)
     mexErrMsgTxt("Expecting 'hi_res', 'lo_res', 'element', and 'time'\n");
   
  m = mxGetM(prhs[1]); n= mxGetN(prhs[1]);
  no_samples = m; 

  if ( n!= flc->no_focus_time_lines )
     mexErrMsgTxt("The number of columns of 'hi_res' must be equal to the number of lines.\n");
       
  
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m!= no_samples) 
      mexErrMsgTxt("The two images must have the same number of samples per line\n");
 
  if ( n!= flc->no_focus_time_lines )
     mexErrMsgTxt("The number of columns of 'lo_res' must be equal to the number of lines.\n");
  
  
  if (mxGetM(prhs[3])>1 || mxGetN(prhs[3])>1)
     mexErrMsgTxt("Expecting scalar for 'element'");

  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])> 1)
     mexErrMsgTxt("Expecting a scalar for 'time'\n");
  
  element = (ui32)floor(mxGetScalar(prhs[3]))-1;
  time = mxGetScalar(prhs[4]);
  
  hi_res = (double**)malloc(flc->no_focus_time_lines* sizeof(double*));
  assert(hi_res);
  
  lo_res = (double**)malloc(flc->no_focus_time_lines* sizeof(double*));
  if (lo_res == NULL) {free(hi_res); abort();}
 
  
  plhs[0] = mxCreateDoubleMatrix(no_samples,flc->no_focus_time_lines,mxREAL);
  
  memcpy(mxGetPr(plhs[0]), mxGetPr(prhs[1]), no_samples*flc->no_focus_time_lines*sizeof(double));
  
  ptr1 = mxGetPr(prhs[2]);
  ptr2 = mxGetPr(plhs[0]);
  
  for (i = 0; i < flc->no_focus_time_lines; i++){
     lo_res[i] = ptr1 + no_samples*i;
     hi_res[i] = ptr2 + no_samples*i;
  }
  
  add_images(flc, salc, &sys, hi_res, lo_res, element, time, no_samples);

  free(hi_res);
  free(lo_res);
}



/*******************************************************************
 * FUNCTION : bft_add_images
 * ABSTRACT : Sum 2 low resolution images in one high resoltuion.
 *******************************************************************/ 
void bft_sub_images(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  ui32 element;
  ui32 no_samples;
  ui32 m,n,i;
  double **lo_res;
  double **hi_res;
  double *ptr1, *ptr2;
  double time;

  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
  
  if (nrhs != 5)
     mexErrMsgTxt("Expecting 'hi_res', 'lo_res', 'element', and 'time'\n");
   
  m = mxGetM(prhs[1]); n= mxGetN(prhs[1]);
  no_samples = m; 

  if ( n!= flc->no_focus_time_lines )
     mexErrMsgTxt("The number of columns of 'hi_res' must be equal to the number of lines.\n");
       
  
  m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
  if (m!= no_samples) 
      mexErrMsgTxt("The two images must have the same number of samples per line\n");
 
  if ( n!= flc->no_focus_time_lines )
     mexErrMsgTxt("The number of columns of 'lo_res' must be equal to the number of lines.\n");
  
  
  if (mxGetM(prhs[3])>1 || mxGetN(prhs[3])>1)
     mexErrMsgTxt("Expecting scalar for 'element'");

  if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])> 1)
     mexErrMsgTxt("Expecting a scalar for 'time'\n");
  
  element = (ui32)floor(mxGetScalar(prhs[3]))-1;
  time = mxGetScalar(prhs[4]);
  
  hi_res = (double**)malloc(flc->no_focus_time_lines* sizeof(double*));
  assert(hi_res);
  
  lo_res = (double**)malloc(flc->no_focus_time_lines* sizeof(double*));
  if (lo_res == NULL) {free(hi_res); abort();}
 
  
  plhs[0] = mxCreateDoubleMatrix(no_samples,flc->no_focus_time_lines,mxREAL);
  
  memcpy(mxGetPr(plhs[0]), mxGetPr(prhs[1]), no_samples*flc->no_focus_time_lines*sizeof(double));
  
  ptr1 = mxGetPr(prhs[2]);
  ptr2 = mxGetPr(plhs[0]);
  
  for (i = 0; i < flc->no_focus_time_lines; i++){
     lo_res[i] = ptr1 + no_samples*i;
     hi_res[i] = ptr2 + no_samples*i;
  }
  
  add_images(flc, salc, &sys, hi_res, lo_res, element, time, no_samples);

  free(hi_res);
  free(lo_res);
}


/*****************************************************************
 * FUNCTION : bft_filter
 * ABSTRACT : Call the function "set_filter". This function is
 *            used to set the filter bank.
 *****************************************************************/
void bft_filter(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   double *coef;
   ui32 Ntaps;
   ui32 Nf;
   int m, n;
   
  if (!initialized)
      mexErrMsgTxt("\nToolbox is not initialized\n");
  
  if (nrhs != 4)
     mexErrMsgTxt("Expecting 'Nf', 'Ntaps', ' and 'h'\n");
   
  if (mxGetM(prhs[1]) > 1 || mxGetN(prhs[1]) > 1)
     mexErrMsgTxt("'Nf' must be scalar");
  
  if (mxGetM(prhs[2])>1 || mxGetN(prhs[2]) > 1)
     mexErrMsgTxt("'Ntaps' must be scalar");
  
  m=mxGetM(prhs[3]); n=mxGetN(prhs[3]);
  if (m > 1 && n > 1)
     mexErrMsgTxt("The impulse response 'h' must be a vector");
  
  Nf = (ui32)floor(mxGetScalar(prhs[1]) + 0.5);
  Ntaps  = (ui32)floor(mxGetScalar(prhs[2]) + 0.5);
  if ((m*n)!=Nf*Ntaps){
     eprintf("The size of  'h' is [%d x %d]\n", m, n);
     eprintf("Ntaps = %d\n Nf = %d \n", Ntaps, Nf);
     mexErrMsgTxt("It is necessary that \"m*n == Ntaps*Nf\"");
  }
  coef = mxGetPr(prhs[3]);
  set_filter_bank(flc, Nf, Ntaps, coef);
}

/*******************************************************************
 * FUNCTION  : bft_delay
 * ABSTRACT  : Delay a line.
 *******************************************************************/
void bft_delay(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
   double *times;
   double *delays;
   ui32 no_delays;
   double *src;
   ui32 src_no_samples;
   double src_start_time;
   double dest_start_time;
   ui32 dest_no_samples;
   double* dest;
   int m,n;
   
   double* ptr_to_result;
   
   PFUNC
   if (nrhs!=7){
      eprintf("%s","Expecting : 'in_line', 'delays', 'times', 'input_start_time' \n");
      eprintf("%s","            'output_start_time', 'no_output_samples'");
      mexErrMsgTxt("");
   }
   
   if (mxGetM(prhs[1])>1 && mxGetN(prhs[1])>1)
       mexErrMsgTxt("'in_line' must be a vector");
       
   if (mxGetM(prhs[2])>1 && mxGetN(prhs[2])>1)
       mexErrMsgTxt("'delays' must be a vector");
   
   if (mxGetM(prhs[3])>1 && mxGetN(prhs[3]))
       mexErrMsgTxt("'times' must be a vector");
   
   if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1)
       mexErrMsgTxt("'input_start_time' must be a scalar value");
   
   if (mxGetM(prhs[5])>1 || mxGetN(prhs[4])>1)
       mexErrMsgTxt("'output_start_time' must be a scalar value");
   
   if (mxGetM(prhs[6])>1 || mxGetN(prhs[6])>1)
       mexErrMsgTxt("'no_output_samples' must be a scalar");
   
   
   /*
    *  Argument #1 is a vector with the samples in the signal
    */
   m = mxGetM(prhs[1]); n = mxGetN(prhs[1]);
   src_no_samples = m*n;
   src = mxGetPr(prhs[1]);
   
   /*
    *  Argument #2 is a vector with the delays
    */
    
   m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
   no_delays = m*n;
   delays = mxGetPr(prhs[2]);
   
   /*
    *  Argument #3 is a vector with the times after which the 
    *  corresponding delay is active
    */
   m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
   if (no_delays != m*n){
      eprintf("The vector 'delays' has %d components \n", no_delays);
      eprintf("The vevtor 'times' has %d components \n", m*n);
      eprintf("These two vectors must be the same size \n");
      mexErrMsgTxt("");
   }
   

   times = mxGetPr(prhs[3]);
   

   /*
    *  Argument #4 is the starting time of the input signal
    */
   src_start_time = mxGetScalar(prhs[4]);
   

   /*
    *  Argument #5 is the starting time of output signal
    */
   
   dest_start_time = mxGetScalar(prhs[5]);
   
   

   /*
    *  Argument #6 is the number of samples in the output signal
    */
   
   dest_no_samples = (ui32)floor(mxGetScalar(prhs[6]));
   
   
   dest = delay_line_linear(&sys, times, delays, no_delays, src,
                             src_no_samples, src_start_time,
                             dest_start_time, dest_no_samples);
   if (dest!=NULL){
        plhs[0] = mxCreateDoubleMatrix(dest_no_samples,1,mxREAL);
        ptr_to_result = mxGetPr(plhs[0]);
        memcpy(ptr_to_result, dest, dest_no_samples * sizeof(double));
        free(dest);
    }else{
     mexErrMsgTxt("Could not apply delays");
   }                             
   
}

/*******************************************************************
 * FUNCTION : bft_delay_filter
 * ABSTRACT : Apply a delay over a whole line, usign a linear phase
 *            filter
 *******************************************************************/
void bft_delay_filter(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   double *times;
   double *delays;
   ui32 no_delays;
   double *src;
   ui32 src_no_samples;
   double src_start_time;
   double dest_start_time;
   ui32 dest_no_samples;
   double* dest;
   int m,n;
   
   double* ptr_to_result;
   
   PFUNC
   if (nrhs!=7){
      eprintf("%s","Expecting : 'in_line', 'delays', 'times', 'input_start_time' \n");
      eprintf("%s","            'output_start_time', 'no_output_samples'");
      mexErrMsgTxt("");
   }
   
   if (mxGetM(prhs[1])>1 && mxGetN(prhs[1])>1)
       mexErrMsgTxt("'in_line' must be a vector");
       
   if (mxGetM(prhs[2])>1 && mxGetN(prhs[2])>1)
       mexErrMsgTxt("'delays' must be a vector");
   
   if (mxGetM(prhs[3])>1 && mxGetN(prhs[3]))
       mexErrMsgTxt("'times' must be a vector");
   
   if (mxGetM(prhs[4])>1 || mxGetN(prhs[4])>1)
       mexErrMsgTxt("'input_start_time' must be a scalar value");
   
   if (mxGetM(prhs[5])>1 || mxGetN(prhs[4])>1)
       mexErrMsgTxt("'output_start_time' must be a scalar value");
   
   if (mxGetM(prhs[6])>1 || mxGetN(prhs[6])>1)
       mexErrMsgTxt("'no_output_samples' must be a scalar");
   
   
   /*
    *  Argument #1 is a vector with the samples in the signal
    */
   m = mxGetM(prhs[1]); n = mxGetN(prhs[1]);
   src_no_samples = m*n;
   src = mxGetPr(prhs[1]);
   
   /*
    *  Argument #2 is a vector with the delays
    */
    
   m = mxGetM(prhs[2]); n = mxGetN(prhs[2]);
   no_delays = m*n;
   delays = mxGetPr(prhs[2]);
   
   /*
    *  Argument #3 is a vector with the times after which the 
    *  corresponding delay is active
    */
   m = mxGetM(prhs[3]); n = mxGetN(prhs[3]);
   if (no_delays != m*n){
      eprintf("The vector 'delays' has %d components \n", no_delays);
      eprintf("The vevtor 'times' has %d components \n", m*n);
      eprintf("These two vectors must be the same size \n");
      mexErrMsgTxt("");
   }
   

   times = mxGetPr(prhs[3]);
   

   /*
    *  Argument #4 is the starting time of the input signal
    */
   src_start_time = mxGetScalar(prhs[4]);
   

   /*
    *  Argument #5 is the starting time of output signal
    */
   
   dest_start_time = mxGetScalar(prhs[5]); 
   
   

   /*
    *  Argument #6 is the number of samples in the output signal
    */
   
   dest_no_samples = (ui32)floor(mxGetScalar(prhs[6]));
   
   
   dest = delay_line_filter(&sys,  &flc->filter_bank, times, delays,
                                   no_delays, src,
                                   src_no_samples, src_start_time,
                                   dest_start_time, dest_no_samples);
   if (dest!=NULL){
        plhs[0] = mxCreateDoubleMatrix(dest_no_samples,1,mxREAL);
        ptr_to_result = mxGetPr(plhs[0]);
        memcpy(ptr_to_result, dest, dest_no_samples * sizeof(double));
        free(dest);
    }else{
     mexErrMsgTxt("Could not apply delays");
   }                             
  
}


/*******************************************************************
 BFT_XDC_SET - Set new coordinates for the transducer elements
 *******************************************************************/
void bft_xdc_set(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{

  TTransducer *xdc;
  TPoint3D *centers;
  ui32 no_elements;
  uint64 address;
  
   
  if (!initialized)
     mexErrMsgTxt("\nToolbox is not initialized\n");
  
     
  if (nrhs!=3)
     mexErrMsgTxt("\nExpecting an array with the centers of the elements\n");

  address = (uint64)mxGetScalar(prhs[1]);
  xdc = (TTransducer*) address;
  
  if (mxGetM(prhs[2])!=3)
     mexErrMsgTxt("The points must be defined with their [x y z] coordinates\n");
     
  no_elements = mxGetN(prhs[2]);
  centers = (TPoint3D*)mxGetPr(prhs[2]);

  bft_transducer_set(xdc, no_elements, centers);
}




/*******************************************************************
 * FUNCTION  : mexFunction 
 * ABSTRACT  : Entry function of the interface between Matlab and
 *             the beamforming library
 *******************************************************************/ 
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
   int function_id;
   

   if (!mexNoEntries){
#ifdef DEBUG
       printf("Debug mode\n");
#endif

#ifdef DEBUG_TRACE
       printf("Tracing functions\n");
#endif

#ifdef SHOW_ENTRIES
       printf("Showing some the input arguments\n");
#endif

#ifdef  SPECIAL_CASE
      printf("Special case \n");
      nice(-20);
#endif

       signal(SIGABRT,atAbort);
       mexAtExit(mexBFTExit);
       mexNoEntries ++;
   }
  
   if( nrhs < 1 )
      mexErrMsgTxt("\nmexFunction\nERROR- needed at least one argument.\n");

   function_id = (int)floor(mxGetScalar(prhs[0]) + 0.5);
   switch(function_id){
       case BFT_INIT: bft_init(nlhs, plhs, nrhs, prhs); break;
       case BFT_END: bft_end(nlhs, plhs, nrhs, prhs); break;
       case BFT_PARAM: bft_param(nlhs, plhs, nrhs, prhs); break;
       case BFT_NO_LINES: bft_no_lines(nlhs, plhs, nrhs, prhs); break;
       case BFT_XDC_FREE: call_bft_free_xdc(nlhs, plhs, nrhs, prhs); break;
       case BFT_TRANSDUCER: call_bft_transducer(nlhs, plhs, nrhs, prhs); break;     
       case BFT_CENTER_FOCUS: bft_center_focus(nlhs, plhs, nrhs, prhs); break;
       case BFT_FOCUS: bft_focus(nlhs, plhs, nrhs, prhs); break;
       case BFT_FOCUS_TIMES: bft_focus_times(nlhs, plhs, nrhs, prhs); break;
       case BFT_APODIZATION: bft_apodization(nlhs, plhs, nrhs, prhs); break;
       case BFT_DYNAMIC_FOCUS: bft_dynamic_focus(nlhs, plhs, nrhs, prhs); break;
       case BFT_BEAMFORM: bft_beamform(nlhs, plhs, nrhs, prhs); break;
       case BFT_SUM_IMAGES : bft_sum_images(nlhs, plhs, nrhs, prhs); break;
       case BFT_ADD_IMAGES : bft_add_images(nlhs, plhs, nrhs, prhs); break;
       case BFT_SUM_APODIZATION : bft_sum_apodization(nlhs, plhs, nrhs, prhs); break;
       case BFT_SUB_IMAGES: bft_sub_images(nlhs, plhs, nrhs, prhs); break;
       case BFT_FOCUS_2WAY: bft_focus_2way(nlhs, plhs, nrhs, prhs); break;
       case BFT_FOCUS_PIXEL: bft_focus_pixel(nlhs, plhs, nrhs, prhs); break;
       case BFT_FILTER:  bft_filter(nlhs, plhs, nrhs, prhs); break;
       case BFT_DELAY:   bft_delay(nlhs, plhs, nrhs, prhs); break;
       case BFT_DELAY_FILTER: bft_delay_filter(nlhs, plhs, nrhs, prhs); break;
		 case BFT_XDC_SET: bft_xdc_set(nlhs, plhs, nrhs, prhs); break;
		 
       default: printf("\007 mexFunction :\n");
                printf("Unknown function id. \n");
   }
 
} 
