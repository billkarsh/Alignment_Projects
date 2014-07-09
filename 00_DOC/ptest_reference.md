# 1_DMesh a.k.a "ptest"

This document covers:

* **ptest.exe**: The main image-pair matcher and point generator.
* **matchparams.txt**: The main parameters controlling ptest.
* **1_PTestx**: Environment for ptest experiments and parameter tuning.

Historically, "dmesh" simply meant deformable mesh, and "ptest", or pair-test was a test application for trying the deformable mesh method to see how well it worked, and of course, for trying a pair of images to see if they were amenable to the method. The experimental application never went away and the original names simply stuck; that's how it is.

We will divide the document into five main sections that follow the workflow through ptest:

1. Parsing the job parameter hierarchy
2. Loading and preprocessing images
3. Finding an approximate match using FFTs
4. Refining the match by mesh optimization
5. Sanity checking the result

## Job Parameters

### General Program Anatomy

Most of my sources have an internal organization that I've adhered to for decades. Reading from top to bottom:

**Data Stuff:**

* includes
	* app-specific headers
	* my library headers
	* 3rd-party headers
	* std library headers
	* using directives
* constants and enums
* pragmas
* macros
* types
* globals
* statics

**Function Stuff:**

* "static" and "private" functions
* functions that call things above
* "public" methods and functions
* main

I try to keep functions pretty simple, main() especially, so you get a flavor of module workflow from the top-level functions (at the bottom of the file, of course). By the way, I almost always order a function's parameters like this:

```
returncode = foo( outputs, inputs );
```

When the returned value is a success indicator, I use the convention `0=failure, !0=success`.

Most main() functions begin with a call to **C_myglobals_::**`SetCmdLine()` which will typically do the following:

* Set parameter defaults.
* Open an output (*.log) file.
* Parse and check command line parameters...
* Turn those into C_myglobals_ member data.

### And now, ptest

The main() function for ptest lives in `dmeshdriver.cpp` and the SetCmdLine() method is in file `CGBL_dmesh`.

#### PWD

ptest must be run within an S-type or D-type subfolder of a temp directory structure. Moreover, the temp directory must contain two files:

* `imageparams.txt`: Contains absolute path to IDB
* `matchparams.txt`: Embodies default job parameters

>Within an alignment pipeline context the temp directory is properly configured by the `mongo.sht` script. To try ptest outside (or alongside) the pipeline context use the wrapper `1_PTestx/ptestx` which first creates a temp directory and then calls ptest. We'll cover ptestx in later sections of this document.

#### Image Precedence

The ptest command line is `ptest <za.ia^zb.ib> [options]`. ptest aligns image A to image B, that is, finds transforms from the A-coordinate system to that of B. The image labels `zk.ik` refer to **layer.tile-id**. **_There is currently no command-line specification of image subregion indices_**. The zk.ik labels are used to:

* Label logs and reports.
* Select context-dependent parameters (za = zb) vs (za > zb).
* Reference IDB images, masks and values.

If no command-line image options are specified, the labels look up paths to images and masks from the IDB. That can be overridden on an item-by-item basis using:

* -ima=path: path to image A (or relative to S- or D-type folder)
* -imb=path: path to image B
* -fma=path: path to foldmask A
* -fmb=path: path to foldmask B

#### Parameter Precedence

>Note that many matchparams variables names explicitly include tag **_SL** or **_XL** to separately handle **same layer** or **cross layer** cases. Note too, that a few parameters are applicable only to the cross layer case yet their names do not have a tag, e.g. **{SCALE, XSCALE, YSCALE, SKEW, MTA}**. I'll claim that's for backward compatibility.

Generally the default parameters are set by the matchparams.txt file. SetCmdLine() compares za to zb and sets the context-dependent globals subgroup `GBL.ctx` using the corresponding _SL or _XL variants. Subsequent code need not retest which is in effect.

The file defaults can be overridden by command-line options. This is a very straightforward 1:1 substitution for most parameters available on the command line (see SetCmdLine() code or 1_Ptestx/pgo.sht).

It's sometime useful to apply a known deformation (always applied to image A) before aligning, say, to describe a tilted image plane. A complicated precedence hierarchy has evolved to allow a variety of specifications of **_strictly affine_** deformation. Here's the precedence from highest to lowest:

1. -Tdfm specified on command line
2. {-SCALE, -XSCALE, -YSCALE, -SKEW} 'factors' on command line (SL or XL)
3. {SCALE, XSCALE, YSCALE, SKEW} 'factors' in matchparams.txt (XL ONLY)
4. TAB2DFM in matchparams (extracts deformation from Tab)
5. Unity (default).

To further clarify this, the fullest description possible is a complete affine Tdfm, which can only be given on the command line. If you supply that it trumps anything else.

Short of that, you may have one or more known scale factors or skew, in which case, the given combination {SCALE, XSCALE, YSCALE, SKEW} of factors is composed into an effective Tdfm using `Taffine::ComposeDfm()` which sets `Tdfm = Scale x (Yskew x Xskew)`.

Even if you haven't given any explicit deformation via Tdfm or factors, there is always a starting guess Tab which can either be implicit (taken from IDB transforms) or explicitly given on the command line. Internally, we decompose it: `Tab = Trgd x Tdfm`, where Trgd contains only effective rotation and translation and Tdfm is a pure deformation matrix. You can elect to retain the Tdfm part of Tab by setting `matchparams::TAB2DFM=Y`.

A final precedence twist is that you can override the effective angle of Tab by setting option `-CTR=angle`; very useful for debugging.

## Image Loading & Preprocessing

### Your Own Preprocessing

Of course, this is content-based matching and the pipeline has a few tricks available to improve the signal to background but the results will be better still if you use the best images you can and fix what you can in advance.

* **Software Lens Correction**: If you know there are repeatable geometric distortions in the images (due to the imaging optics, say) you might try this tool: [Distortion_Correction](http://fiji.sc/Distortion_Correction).

* **Sparse Content**: This is primarily a problem in fluorescence data where the features or the labeling may be sparse. Usually you'll have several color channels so obviously you want to align on the densest of them. In some cases it may be useful to superpose multiple channels to increase feature density. Finally, you might apply a **histogram equalization filter** to boost low intensity levels.

### CPixPair::Load()

After `dmeshdriver::main()` calls SetCmdLine() to sort out parameters it creates a `0_GEN/CPixPair` object 'px' and calls `px.Load()`. This loads and optionally filters both A and B. We'll point out notable constraints and highlights here.

#### Dimensions

Immediately upon loading the A and B rasters we check the dimensions, that is, we require wa == wb and ha == hb (though width and height are independent of each other). **_As a general rule, then, all images in a stack must have the same dimensions_**. The natural full dimensions are stored in the px object as `{wf, hf}`.

We find that reducing the image scale by a factor of two or four doesn't significantly degrade overall accuracy, but saves substantially on memory footprint and compute time. Therefore, we iteratively downsample the images (divide each dimension by 2) until each dimension is <= `MAX1DPIX = 2048`. Most subsequent work will be done at this scaled down size `{ws, hs}` with integral scale factor `scl = wf/ws`. Therefore, another data requirement is that **_each dimension must be divisible by two or four (if larger than 2048)_**.

#### Background Subtraction

After raster loading comes intensity flattening which we do by projecting the intensity onto low order Legendre polynomials (independently for X- & Y-axes) to compute those coefficients (one pass) and then subtracting out those components (second pass). The zero order poly is a constant, so models a detector pedestal value. The first order is a line so these can model stage tilt. The second order is a half-cosine which models intensity fall-off toward the edges (vignetting). Higher orders sometimes seem promising but have dubious physical interpretation.

* Orders up to 6 are supported.
* Select highest order with `matchparams::PXBRO`.

#### Resin Masking

For EM images, use `matchparams::PXRESMSK=Y` to enable the generation of mask images that discriminate tissue (mask value 1) from resin (0). We make masks with utility function `0_GEN/Maths::ResinMask8()`. The algorithm is entirely heuristic but has proven very effective on EM images from different labs/scopes here at Janelia:

1. Downsample by 8X or so.
2. Apply Sobel edge filter 3 times to fill space with edges.
3. Smooth with a median filter.
4. Threshold at some 8-bit level (100 or so) to make mask.
5. Reascale mask to image dimensions.

There's a Boolean parameter for same-layer Y/N. If same-layer we use a smaller median kernel and a lower threshold. These choices tend to keep a bit more fine detail in the mask which can be very helpful in finding something to match on the periphery of fly brains where there are thin strands of tissue remaining from excising the brain. In the same-layer case, these tiny features actually exist in both images so make genuine matches. If not same-layer tiny features are unlikely to match well so we prefer to smooth them away.

>Note: This function additionally zeroes mask pixels where the image intensity is either zero or saturated at 255, which are unlikely to be real tissue values.

Getting back to px.Load(), the first thing we do is make resin masks with `same-layer = false` in order to get an estimate of how much image area is real tissue (sum all 1-pixels). We simply refuse to match these images if both of them are below 15% tissue. Remember the rule of conservatism in alignment by features: **_It is far worse to include wrong matches than to miss a few good ones...one only needs sufficiently many good points for a solution._**

>Note: If you're reading along in the code, you'll see that if this tissue fraction test fails, we log a message with format `FAIL: Reason`. It is our fairly uniform practice to tag key messages with `FAIL` in all-caps for easy grep searching.

#### Image Pointers

Class CPixPair publishes pointers to raster images, named with tags `aln` or `vfy` according to whether we recommend that image be used for determining alignment (with enhancement filters applied) or for qualifying the result and measuring resulting correlation.

Pointers are used because the scaling and filtering operations are optional or dynamically determined so in many cases, for example, the filtered and natural version may be the same underlying storage location.

#### Spatial Filtering Experiments

Although I've played repeatedly with increasing correlation using spatial filters like **DoG**, I can't say I was ever impressed with any brand of that. However, CPixPair is a ready home for filtering experiments should you feel so compelled.

You'll also notice some code in here for activating a built-in flavor of 'software lens support'. That was a simple experiment in infrastructure; only taken as far as per-camera affine corrections. Affine correction in particular isn't very potent considering all the other ways of specifying affine correction in ptest (**pre-tweak affines** is another to be covered shortly). The method of lens correction we **_do_** advocate is externally generating new corrected images for alignment before the pipeline ever begins. By this means all stages of the pipeline enjoy correction without needing to reapply it.

In summary, while resin filtering is a really big win, we enable neither magic kernels nor the crude lens 'feature'.

#### Intensity Normalization

Take note that, whatever combination of operations is applied to the pixels, we always call `0_GEN/Maths::Normalize()` afterward. This function first calculates the current mean and standard deviation of a value array and then rescales the array to `mean = 0, stdev = 1`. What we do with pixels for the most part is to pair-up A-pixels with B-pixels and then measure their correlation coefficient (R) which is essentially some flavor of **SUM<sub>i=0,N </sub>( A<sub>i</sub>B<sub>i</sub> )**. The idea behind the normalization is that middle tones are common and boring whereas extreme values (big contrasts) are more rare and interesting. Hence, we down-weight average values by mapping them to zero so their contribution to R is minimal.

>Note: You'll see in `0_GEN/Maths` that we have several alternative ways of 'normalizing' wherein we prefer central values or the tails, and so on, groping for ways to enhance signal-to-background in R sums. This has not been fruitful in our hands.

### Foldmask Loading & Preprocessing

Let's pick up the flow again inside `dmeshdriver::main()`. After px.Load() returns, we use the scaling factor px.scl (result of the 2048 size cap) to adjust some of the matchparams geometric parameters (which are always expressed in that file in units of full-size image pixels). Don't worry, LIMXY is appropriately scale-adjusted where it's used.

>The next thing you see in main() is an optional clause to call `InSectionOverlap()` and then exit. This is an archaic experiment to try (for two images having no relative rotation, i.e., in the same layer) the full course of FFT approximate displacement followed by mesh optimization. I'll probably retire this at some point because it's not super instructive; you can safely ignore this. It has always been difficult to choose whether to scrub away history and make it pretty for you, or to keep snippets of the development path as reminders of how Lou and I got here.

After scaling adjustments, action resumes in `CalcTransforms()`; a static function higher up in this file. We're almost ready to start calculating. The last thing we have to do is to prepare the foldmasks, one for each of the two images.

#### What's a Foldmask?

A serial section may suffer damage {tears, folds, other} that divide it into two or more distinct subregions, and these subregions may have moved relative to one other by being picked up from the slicer differently or adhering to the substrate differently. Consequently, each subregion would need its own coordinate transform. A foldmask is a metaimage of the same XY dimensions (full scale) as the original. The 8-bit (unsigned) foldmask pixels identify which subregion the original pixels belong to, using a simple 4-neighbor connectivity rule.

Special foldmask value zero identifies any damage region that should be ignored in all alignment calculations. All pixels having value one are 4-way connected to tissue region one, but have no connection to pixels of value two or any other, and so on for foldmask values two through 255. The rule that, **_distinct subregions within a section are not connected_**, is what allows them to float free of one another and obtain their own solution transforms. Of course when aligning a subregion to an adjacent layer, any pairing of an A-subregion onto a B-subregion is permissible.

There is important missing machinery in our pipeline surrounding foldmasks. Folds, tears and the like are defects that occur at the level of whole layers (montages), that is, the whole montage might be broken into some handful of distinct subregions. Typically a layer is tiled by a dozen to several thousand image tiles. Montage subregion three may appear in hundreds of images and it would be extremely valuable if foldmask value three referred to the same connected montage region in all images for this layer. However, although we have an automated tool called `1_Tiny` that does a fair job of calculating a foldmask for an image, it assigns numbers to found subregions within a given image independently of all other images. There is no existing process that coordinates the numbering across image tiles. Although we have a degree of foldmask "awareness" in most parts of the pipeline, the ad hoc numbering issue is a huge gap that needs development to exploit folmasks fully. I never got to that before I had to move on to other projects.

#### Using Foldmasks to Remove Resin

Returning now to the loading of foldmasks...although labeling connected subregions is imperfect, ptest nevertheless associates a foldmask with every image whether or not the no-folds `-nf` option is invoked. This is handled by `0_GEN/FoldMask::GetFoldMask()`. If using real foldmasks such as those generated by the tiny program, the masks will be loaded from disk files stored in the IDB. If the `-nf` option is specified, we simply allocate a buffer and fill it with ones signifying that the entire image is 4-way connected as a single region.

Now comes a beautiful trick. Every image now has a mask where zero means the pixel is not used and non-zero means it is used. It is entirely consistent with this labeling scheme, to take the logical AND (product) of this mask with any other binary mask wherein zero and one similarly signify not/used. That's exactly how we use the resin masking option to knock out resin pixels. If `matchparams::PXRESMSK=Y` then px.Load() calculates its resmsk members making them non-empty and they are sent as parameters to GetFoldMask.

#### Rectangular Cropping Regions

The same mask intersection method can be used to remove the margins of images; useful because the extreme periphery of images may suffer especially large aberrations. To invoke this mechanism, one creates a cropping file which is a simple text file having a line per referenced camera index with these fields (whitespace separated):

```
cam-id x0 y0 dx dy
```

Remember that in your initial layout metadata file you supply a zero-based camera id number for every image (up to four cameras). Next, before you execute your `dbgo.sht` script to make an IDB, edit the script to include option `-crop=mycroprectfile.txt`.

When dbgo.sht runs it copies your file into the top level of the IDB, giving it the standardized name `crop.txt`. Subsequently, any pipeline component like ptest can check for and implement your rectangles.

## Aside: I/O and Loops Over Subregions

Now that the parameter, image and foldmask data are loaded, CalcTransforms() calls `PipelineDeformableMap()` to do all the work, both thumbnail matching and mesh optimization. Find this function at the bottom of file `1_DMesh/dmesh.cpp`. This function's inputs are easy to locate:

* func-param : const CPixPair &px : all the image data
* func-param : const uint8 *foldmaskA
* func-param : const uint8 *foldmaskB
* func-param : FILE* flog : handle to output text log `.../temp/za/S(D)x_y/pair_za.ia^zb.ib.log`.
* global : GBL : pointer to global parameters class via `#include "CGBL_dmesh.h"`

## Approximate Thumbnail Matching Using FFTs

_fin_

