/*
 * Code for aligning to images through piecewise affine transformations.
 * The images may have "folds", in which case the mesh is discontinuous across
 * the folds.
 *
 * Algorithm for deformable mesh based alignment is contained in dmesh.cpp/h. This
 * was written by Louis Scheffer, Visiting Scientist, Janelia Farm Research Campus, HMMI.
 *
 * Wrapper for MATLAB written by Shiv Vitaladevuni, JFRC, HHMI
 *
 * v0   12102008  init. code
 */

#include	<stdio.h>

#include	<mex.h>

#include	"tile_align_deformable_mesh.h"


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("Usage [map_mask, transforms] = tile_align_deformable_mesh(image_1, ...\n");
    mexPrintf("\t\tfold_mask_1, image_2, fold_mask_2)\n");
    mexPrintf("Inputs:\n");
    mexPrintf("1. image_1      1st image, MxN uint8 matrix\n");
    mexPrintf("2. fold_mask_1  fold mask for image_1, 0 for fold, 1,...,C for connected components,MxN uint8 matrix\n");
    mexPrintf("3. image_2      2nd image, MxN uint8 matrix\n");
    mexPrintf("4. fold_mask_2  fold mask for image_2, 0 for fold, 1,...,C for connected components,MxN uint8 matrix\n");
    mexPrintf("5. within_section_flag    set to 1 if both tiles belong to the same section\n");
    mexPrintf("Outputs:\n");
    mexPrintf("1: map_mask    map of transforms, <10 no mapping, 10+i use ith transform, MxN uint16 matrix\n");
    mexPrintf("2: tranforms   array of affine transforms, one for each 10+i in map_mask, 6xT double matrix\n");
    return;
  }
  if(nrhs!=5 && nrhs!=6){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  if(nlhs!=2){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  const mxArray * image_1_mx = prhs[0];
  const mxArray * fold_mask_1_mx = prhs[1];
  const mxArray * image_2_mx = prhs[2];
  const mxArray * fold_mask_2_mx = prhs[3];
  const mxArray * is_within_section_mx = prhs[4];
  char *params;
  if(nrhs>=6) {
    int length = mxGetN(prhs[5])+1;
    params = (char *) calloc(length, sizeof(char));
    mxGetString(prhs[5], params, length);
    mexPrintf("Got option string '%s'\n", params);
    }
  else
    params = (char *)calloc(1, sizeof(char));  // make an empty string

    const int * sizeImage;
    sizeImage = mxGetDimensions(image_1_mx);
    int width = sizeImage[1], height = sizeImage[0];

  unsigned char * image_1 = (unsigned char *) mxGetPr(image_1_mx);
  unsigned char * fold_mask_1 = (unsigned char *) mxGetPr(fold_mask_1_mx);
  unsigned char * image_2 = (unsigned char *) mxGetPr(image_2_mx);
  unsigned char * fold_mask_2 = (unsigned char *) mxGetPr(fold_mask_2_mx);
  int is_within_section = (int) * mxGetPr(is_within_section_mx);

  plhs[0] = mxCreateNumericMatrix(height, width, mxUINT16_CLASS, mxREAL);
  unsigned short int * map_mask = (unsigned short int *) mxGetPr(plhs[0]);
  int n_transform = 0;
  double * transforms = NULL;
  FILE * f_out_log = fopen("dmesh_log.txt", "wt");
  // Configure deformable mesh
  SetDMeshParams(params);
  free(params);

  // Then run the pipeline
  printf("######### Entering PipelineDeformableMap ########\n");
  PipelineDeformableMap(
    height, width,
    image_1, fold_mask_1,
    image_2, fold_mask_2,
    map_mask, n_transform, transforms, stdout, f_out_log);
  printf("######### Exited PipelineDeformableMap ########\n");
  fclose(f_out_log);
  printf("n_transform = %d\n", n_transform);
  plhs[1] = mxCreateDoubleMatrix(6, n_transform, mxREAL);
  double * transforms_out = mxGetPr(plhs[1]);
  {
    int i,j,k;
    k=0;
    for(i=0; i<n_transform; i++){
      printf("Transform %d:\n", i);
      for(j=0; j<6; j++, k++){
        transforms_out[k] = transforms[k];
        printf("%g ", transforms[k]);
      }
      printf("\n");
    }
  }

  if(n_transform>0)
    delete [] transforms;

  return;
}
