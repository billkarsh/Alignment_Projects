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

#ifndef _TILE_PAIR_WITHIN_SECTION_OVERLAP_
#define _TILE_PAIR_WITHIN_SECTION_OVERLAP_

// The main routine the pipeline should call.
void InSectionOverlap(
  int w, // width
  int h, // height
  unsigned char *pa, // image 1
  unsigned char *pb, // image 2
  int &Npts, // number of correspondence points found
  double * & apts, // corrdinates of correspondences in image 1
  double * & bpts, // corrdinates of correspondences in image 2
  FILE *flog // log text file pointer
  );
#endif
