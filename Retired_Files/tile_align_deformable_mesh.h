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

#ifndef _TILE_ALIGN_DEFORMABLE_MESH_
#define _TILE_ALIGN_DEFORMABLE_MESH_

// For setting parameters for deformable mesh algorithm
void SetDMeshParams(const char *arg);

// The main routine the pipeline should call.
void PipelineDeformableMap(
  int w, int h, // size of all images
  unsigned char *AboveRaster,                // the higher layer
  unsigned char *fold_mask_above,            // 0 for fold, 1=first connected region, 2 = second, etc.
  unsigned char *BelowRaster,                // the lower layer
  unsigned char *fold_mask_below,
  unsigned short int *map_mask,                  // the resulting map.  <10 means no mapping, 10+i means use transform i
  int &Ntrans,                       // how many tranforms result?  (returned value)
  double * &array_of_transforms,     // array of these values.
  FILE *fout,			     // write most stuff here
  FILE *flog);                       // write logging stuff here

#endif
