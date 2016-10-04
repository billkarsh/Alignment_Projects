// Given a piecewise affine tranformation from an image plane (I) to
// another plane (L) and a label map on L, find for each pixel in I,
// the corresponding label. Label 0 is assigned for undefined pixels.
// If a mask is also provided, then a label is assigned only if the
// transform id in the mask is equal to the piecewise transform id
// for the pixel.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	04072009	init. code
//

#include	<mex.h>
#include	<math.h>
#include	<algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  mexPrintf("START: get_label_piecewise_affine_transformation_mask\n");
    if(nrhs==0){
        if(nlhs==1){
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
            *mxGetPr(plhs[0]) = 1;
            return;
        }
        mexPrintf("Usage mapped_label = get_label_piecewise_affine_transformation_mask(label_map, transform_id_map, transforms, mask);\n");
        mexPrintf("input params:\n");
        mexPrintf("\t1. label map that is to be mapped by the piecewise affine transformation PxQ double\n");
        mexPrintf("\t2. transformation ids for the piecewise affine transform MxN uint32\n");
        mexPrintf("\t3. list of affine transforms 6xT double\n");
    mexPrintf("\t4. default value scalar double\n");
        mexPrintf("\t5. mask with transformation ids MxN uint32 (optional)\n");
        mexPrintf("output:\n");
        mexPrintf("\t1. mapped labels PxQ double\n");
        return;
    }
    if(nrhs<3 || nrhs>6){
        mexErrMsgTxt("Wrong number of inputs\n");
        return;
    }

    int numDim = mxGetNumberOfDimensions(prhs[0]);
    if(numDim>2){
        mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
        return;
    }

  const mxArray * label_map_mx = prhs[0];
  const mxArray * transform_map_mx = prhs[1];
  const mxArray * transforms_mx = prhs[2];

  const mxArray * default_value_mx = NULL;
  if(nrhs>=4)
    default_value_mx = prhs[3];

  bool is_masked = false;
  const mxArray * mask_mx = NULL;
  if(nrhs>=5){
    is_masked = true;
    mask_mx = prhs[4];
  }

    const int * sizeCanvas;
    sizeCanvas = mxGetDimensions(label_map_mx);
    int width_canvas = sizeCanvas[1], height_canvas = sizeCanvas[0];
  mexPrintf("width_canvas:%d, height_canvas:%d\n", width_canvas, height_canvas);
    const int * sizeImage;
    sizeImage = mxGetDimensions(transform_map_mx);
    int width_image = sizeImage[1], height_image = sizeImage[0];
  mexPrintf("width_image:%d, height_image:%d\n", width_image, height_image);

  double * label_map = mxGetPr(label_map_mx);
  unsigned int * transform_map = (unsigned int *) mxGetPr(transform_map_mx);
  double * transforms = mxGetPr(transforms_mx);

  double default_value = 0;
  if(default_value_mx!=NULL)
    default_value = * mxGetPr(default_value_mx);

  unsigned int * mask = NULL;
  if(is_masked)
    mask = (unsigned int *) mxGetPr(mask_mx);

  plhs[0] = mxCreateNumericMatrix(height_image, width_image, mxDOUBLE_CLASS, mxREAL);
  {
    double * mapped_label = mxGetPr(plhs[0]);
    int x, y;
    unsigned int transform_id, * t;
    double * transform;
    t = transform_map;
    int x1, y1;
    for(x=0; x<width_image; x++){
      for(y=0; y<height_image; y++, t++){
        mapped_label[x*height_image+y] = default_value;
        transform_id = *t;
        if(transform_id<=0)
          continue;
        transform_id--; // 0 indexed
        transform = transforms + 6*transform_id; // affine transforms
        x1 = (int)(transform[0]*(double)(x+1) + transform[2]*(double)(y+1) + transform[4] - 1);
        y1 = (int)(transform[1]*(double)(x+1) + transform[3]*(double)(y+1) + transform[5] - 1);
//        mexPrintf("x:%d, y:%d -> x1:%d, y1:%d\n", x, y, x1, y1);
        if(x1<0 || x1>=width_canvas || y1<0 || y1>=height_canvas)
          continue;
        if(is_masked)
          if(mask[x1*height_canvas+y1]!=transform_id+1)
            continue;
        mapped_label[x*height_image+y] = label_map[x1*height_canvas+y1];
      }
    }
  }
  mexPrintf("STOP: get_label_piecewise_affine_transformation_mask\n");
    return;
}
