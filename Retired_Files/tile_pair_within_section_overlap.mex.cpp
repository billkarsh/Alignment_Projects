/*
 * Code for computing correspondences between two images with overlap.
 * This is to be used as a backup when the overlap is too small for SIFT feature point
 * detection.
 *
 * Algorithm is contained in dmesh.cpp/h. This was written by Louis Scheffer,
 * Visiting Scientist, Janelia Farm Research Campus, HMMI.
 *
 * Wrapper for MATLAB written by Shiv Vitaladevuni, JFRC, HHMI
 *
 * v0   12222008  init. code
 */

#include	<stdio.h>

#include	<mex.h>

#include	"tile_pair_within_section_overlap.h"


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("Usage [points_1, points_2] = tile_align_deformable_mesh(image_1, ...\n");
    mexPrintf("\t\timage_2)\n");
    mexPrintf("Find correspondences between image_1 and image_2\n");
    mexPrintf("Inputs:\n");
    mexPrintf("1. image_1      1st image, MxN uint8 matrix\n");
    mexPrintf("2. image_2      2nd image, MxN uint8 matrix\n");
    mexPrintf("Outputs:\n");
    mexPrintf("1: points_1    correspondence points in image_1 2xR double matrix\n");
    mexPrintf("2: points_2    correspondence points in image_2 2xR double matrix\n");
    return;
  }
  if(nrhs!=2){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  if(nlhs!=2){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  const mxArray * image_1_mx = prhs[0];
  const mxArray * image_2_mx = prhs[1];

    const int * sizeImage;
    sizeImage = mxGetDimensions(image_1_mx);
    int width = sizeImage[1], height = sizeImage[0];

  unsigned char * image_1 = (unsigned char *) mxGetPr(image_1_mx);
  unsigned char * image_2 = (unsigned char *) mxGetPr(image_2_mx);

  int n_point;
  double * points_1, * points_2;
  FILE * f_out_log = fopen("dmesh_log.txt", "w");
  printf("######### Entering InSectionOverlap ########\n");
  InSectionOverlap(
    height, width,
    image_1, image_2,
    n_point, points_1, points_2,
    f_out_log);
  printf("######### Exited InSectionOverlap ########\n");
  fclose(f_out_log);

  plhs[0] = mxCreateDoubleMatrix(2, n_point, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(2, n_point, mxREAL);
  mexPrintf("n_point = %d\n", n_point);
  if(n_point<=0)
    return;

  double * p1 = mxGetPr(plhs[0]);
  double * p2 = mxGetPr(plhs[1]);
  for(int i=0; i<n_point; i++){
    p1[2*i] = points_1[2*i+1];
    p1[2*i+1] = points_1[2*i];
    p2[2*i] = points_2[2*i+1];
    p2[2*i+1] = points_2[2*i];
  }

  delete [] points_1;
  delete [] points_2;

  return;
}
