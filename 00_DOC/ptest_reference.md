# 1_DMesh a.k.a "ptest"

This document covers:

* **1_DMesh**: The project that **builds ptest.exe**.
* **matchparams.txt**: The main parameters controlling ptest.
* **1_Ptestx**: Environment for ptest experiments and parameter tuning.

Historically, "dmesh" simply meant deformable mesh, and "ptest", or pair-test was a test application for trying the deformable mesh method to see how well it worked, and of course, for trying a pair of images to see if they were amenable to the method. The experimental application never went away and the original names simply stuck; that's how it is.

The document format is pretty simple. I'll run through the algorithm without up-front definition of terms; I think you'll catch on as we go. Following that are several appendix-like sections to add more detail where I think it might be helpful.


* Outline of the Algorithm

Appendices:

* [General Source Code Anatomy](ptest_reference.md#general-source-code-anatomy)
* [Input Parameter Handling](ptest_reference.md#input-parameter-handling)
* [Output (Exhaustive List)](ptest_reference.md#output-exhaustive-list)
* [Image Loading and Preprocessing](ptest_reference.md#image-loading-and-preprocessing)
* [Foldmasks and Connected Regions](ptest_reference.md#foldmasks-and-connected-regions)
* [Approximate FFT Matching](ptest_reference.md#approximate-fft-matching)
* [Notes on Mesh Creation](ptest_reference.md#notes-on-mesh-creation)
* [Ptestx Environment (Try Folder)](ptest_reference.md#ptestx-environment-try-folder)

#### Disclaimer

This software is presented, **_as is_**, in the hopes that it may be useful to you--it has certainly allowed us to achieve very good quality alignment quickly and with only modest effort. Nevertheless, the software had been under active development right up to the time of its publication, hence, inevitably exposes a variety of flaws: {old experiments, deprecated parameters, incomplete logic and no doubt incorrect logic in some areas}. This (small amount of) baggage complicates the code by its presence, but the overtly wrong stuff can be bypassed by suitable parameter choices. On the whole I am comfortable claiming that the utility of this code far outweighs any embarrassment I may suffer.

In particular, I advise the following:

1. Foldmask handling, as I will discuss later in this document, should work quite differently than it does at present. What it does now is peculiar because it reflects shifting ideas over time. I have retained significant amounts of foldmask handling infrastructure because I fully believe these will ultimately play a role in the real thing or at least demonstrate how key pieces might work...But all that machinery is doing little useful work at the moment. Still, what it does now is effectively correct if you turn off foldmasks using at least one of:
	* `matchparams::FLD=N`
	* command line option `-nf`
2. Sanity checking a match by assessing the Fourier spectrum in a difference image may be both poorly grounded and wrongly coded. This is performed by `0_GEN/Metrics::FourierMatch()`. Rather, you should bypass that entirely and instead use the earth mover metric by setting:
	* `matchparams::EMM=Y`.

## Outline of the Algorithm

In a nutshell, it's just this.

1. Use FFT operations to get rough overlap of A & B.
2. Describe A within overlap using deformable mesh of contiguous triangles.
3. Use gradient descent to deform A relative to B, hence, refine A-B correlation.
4. Sanity check results.
5. Report results: paired A-B point matches + diagnostics.

Once more with details...

### How it Works (Ignoring Foldmasks)

For simplicity, we first do a walk through without foldmasks, which means we assume montages (hence all their tiles) are 4-way connected throughout; no tears dividing tissue into pieces that need separate transforms. This is a selectable operation mode `-nf` that works well and is handled in a logically consistent manner throughout the pipeline. It's how we do things in practice.

>I'll revisit foldmasks several times, but basically, we always convert images to point clouds and calculate on those. Enabling foldmasks only affects how we select the pixels that compose a cloud, not how we use clouds thereafter.

### Working Directory

Ptest is always launched from within an 'S' or 'D' type subfolder of a temp directory structure. So its `pwd` will look like `.../temp/za/Sx_y/`. This allows ptest to use relative paths to find and place non-image files. Of particular note, the temp directory must contain two files:

* `imageparams.txt`: Contains absolute path to IDB
* `matchparams.txt`: ptest default job parameters

>Within an alignment pipeline context the temp directory is properly configured by the `mongo.sht` script. To try ptest outside (or alongside) the pipeline context use the wrapper `1_Ptestx/ptestx` which first creates a temp directory and then calls ptest. We'll cover ptestx in later sections of this document.

### Input

Ptest gets as input:

* Two image labels "za.ia^zb.ib" on the command line
* An estimated affine Tab either from the IDB or the command line
* Parameters that tweak behavior: matchparams.txt + command line overrides.

### Output

The output is principally this:

* pts.same(down)

There are many other useful and detailed diagnostic outputs but really the point pairs are what we feed forward down the pipeline. Each ptest job will append to its file some handful of entries like this:

```
CPOINT2 1.12-1 548.347679 716.475949 0.11-1 392.363415 37.427421
CPOINT2 1.12-1 906.183423 788.011835 0.11-1 749.536275 120.037851
CPOINT2 1.12-1 531.969897 905.889984 0.11-1 370.319378 227.175693
...
```

On each line the LHS point is the centroid of a triangle (in a mesh) in A-coords and the RHS point results from operating the associated affine on the LHS point, transforming it into B-coords. There's just one mesh between the two images, so there aren't A-centroids **and** B-centroids.

### Step 1 - Approximate FFT Match (Tab => tr_guess)

Again, ptest gets as input a rough guess "Tab" as a starting point. For a variety of purposes we need a better guess than that at how images A & B overlap before we begin mesh operations, but there are two really important issues:

1. The deformable mesh carves the intersection area into small patches and must be able to measure meaningful correlation changes for them. This becomes more unreliable as patches become smaller. We want patches to be pretty small so we get a good density of triangles, hence, density of reported point pairs to work with in the solver. We can go pretty small as long as we start the mesh off close to a local minimum, so the pressure is on to get a good approximate transform. In addition, you should select rational mesh sizing parameters according to what you know about the contrast and feature density in your data:
	* `matchparams::MNL` : sets nominal triangle edge length scale (same or cross layer)
	* `matchparams::MTA` : sets nominal triangle area (cross layer only).
2. In a moment we'll describe the mesh optimization procedure which expresses the cross correlation of two image patches as a function of mesh vertices (control points), and calculates analytical derivatives of that to decide which direction to step the vertices a la gradient descent. The point here is that correlation isn't generally a smooth function of position, but for typical serial section images (typical feature sizes) it is approximately smooth within a ball of some 20 pixels radius about a true match (that's all we ask of it). So again, a good starting guess is crucial.

I'll just reiterate here that this idea of progressive refinement is central to the success of the pipeline. Looking now at the whole pipeline, we build montages first because that's pretty reliable without global context. Then we assemble montages into a crude stack with the strip matcher which works well because it uses large amounts of real estate (context) to make matches. Then montages are refined block by block. Then blocks are refined tile by tile. Inside ptest now, the tiles are first coarsely aligned with the largest available areas, and then refined on small mesh triangles. Working from large scales to small is what regularizes and stabilizes the process from end to end.

Let's begin. Using the currently recommended mode of angle handling `matchparams::MODE=M or N` we follow these steps:

1. Crop the images to be just a little larger than the predicted overlap region. That's governed by the initial guess Tab and `matchparams::XYCONF`, a [0..1] user supplied estimate of how accurate Tab is. Cropping improves FFT performance. At the same time, the pixels from each image raster are re-expressed in a new point-cloud format consisting of two paired vectors (standard library container for ordered items). One vector contains only image coordinate points, the other contains the intensity values at those points. This method of description accommodates free-form region shapes and allows separately manipulating intensity and geometry.
2. Make thumbnails, that is, down-scale everything by a power of two (usually 4 or 8) for FFT efficiency `matchparams::THMDEC`.
3. Find the best FFT correlation peak within a given disc. That means we create a correlation image wherein pixel position encodes relative displacements between A and B, and the intensity is essentially a correlation coefficient. We apply a binary mask to that describing which parts of the correlation image we should search for the correlation peak. The mask is primarily a disc centered at the nominal displacement obtained by guess Tab, with a radius given by the user supplied estimate of maximal Tab error, `matchparams::LIMXY`. The mask also excludes area in accordance with matchparams specifications for minimum required 1D and 2D image overlap. The peak must exceed the user supplied correlation threshold `matchparams::RTRSH`. If it doesn't we can (optionally via `matchparams::PRETWEAK`) try to automatically apply small amounts of distortion to A and try again. If we can't get over the threshold we quit, reporting no solution for this image pair.
4. Otherwise, we recalculate the final affine at full scale and double-check that it lies within the LIMXY expected disc. We log how far the solution is from the disc center. Large displacements (even though they may be within LIMXY) often indicate bad initial stage coordinates or pathological block-block matches.

The final affine is tabulated in file `.../temp/za/S(D)x-y/ThmPair_za^zb.txt` and is returned to the caller `PipelineDeformableMap()`.

### Step 2 - Point Cloud Intersection (From tr_guess)

From the estimated transform "tr_guess" we again represent the intersection as two paired vectors {one of points and one of values} and we construct these by "pushing_back" all A-points that map via tr_guess into B. We check that the area exceeds the user supplied minimum mesh area `matchparams::MMA`. We then compute the bounding box for the point cloud which will provide a frame in which to erect a mesh.

### Step 3 - Single Triangle Mesh Optimizer (tr_guess => Updated tr_guess)

#### Building a Mesh

Triangle meshes are described by two associated vectors. The first `vector<Point>` lists all the vertices in the mesh (also termed "control points"). The second `vector<triangle>` lists triangle structures whose essential data are really just a triplet of integer indices into the control points list.

Any point in the plane can be expressed as a linear combination of three other non-collinear reference points; a representation called [barycentric coordinates](http://en.wikipedia.org/wiki/Barycentric_coordinate_system). Given a mesh with N control points, we can re-express all the A-points (A-pixels), giving each a vector of N values wherein all values are zero except for the three multipliers corresponding to the vertices of the triangle that contains the given A-point. In this way, we make all of the A-points functions of all N control points. The A-values are just as they were before, and they remain constants throughout the mesh calculations. We've merely changed the list of A coordinate points to a list of constant multipliers, plus a list of control points that will now be treated as variables. By the way, of course we only need three control points to describe all the points in a plane whether they lie inside or outside a triangle. What we've done here is to partition the points among several triangles, always choosing the nearest (enclosing) triangle for each point.

#### Mesh Deformation

Using this representation, as we drag one or more control points this way or that way, the A-pixels follow along. Note too that the continuity of the A point cloud is preserved under control point displacement. We've constructed a deformable mesh. Note that in what follows, the B-pixels will remain in a plain old raster image that remains constant. Only control points (triangle vertices) are variables. That is, the only thing that changes is how the draggable A-pixels line up against the fixed B-pixels.

#### Mesh Optimization

The mesh optimizer measures the correlation of A and B for a given set of control point locations, and it decides how to move control points to get a better correlation (gradient descent). It does this repeatedly until no further improvement is obtained. The correlation measure is just:

**R = SUM<sub>i=0,N</sub>( A<sub>i</sub>B<sub>i</sub> ), N = number of A-pixels**

Peculiar though it seems at first, it is the B-values that are variables in this expression. To evaluate the i-th term of the expression, we use the current floating-point location of A<sub>i</sub> (expressed in B-coords) to interpolate among the fixed B-raster values, thereby obtaining a B<sub>i</sub>-value that depends upon the control points. This is of course multiplied by the always fixed value of A<sub>i</sub>.

The R expression can be differentiated analytically with respect to each of the control point coordinates to get a gradient vector, which we normalize, and then multiply by a step size (initialized at 10 pixels). The step size is decreased as convergence approaches.

#### Baby Steps: Try Just One Triangle

Before erecting a full mesh with lots of small triangles, we first make a mesh with one triangle, simply inscribed in the intersection bounding box. The assignment of points to triangles is done in function `0_GEN/Geometry::BestTriangle()`. The rule is to use the enclosing triangle, or if not enclosed, use the nearest. In this case all points get the same triangle assignment and will move together in the mesh optimization. This rationale for this exercise is twofold. (1) The result is used to improve the FFT-derived affine match. (2) If a reasonable result can't be found using all the real estate we have, then it's unlikely to work in the small triangle case.

Two matchparams override mesh behavior:

* `OPT_SL=N` : In same-layer cases, the optimizer is bypassed, that is, initial control points are not moved. This had application in a rare fluorescence data situation. In instances of spatially sparse bright extended objects the mesh optimizer tended to shrink A-objects to get them to fit wholly inside their B counterparts.
* `ONE=Y` : Bypasses full mesh creation and optimization; the single triangle mesh is used instead.

>Note: An unfortunate naming convention has stuck. In several parts of code, parameters and log output, the single triangle mesh is termed the "**affine**" mesh, as opposed to the multi-triangle mesh which is termed "**deformable**".

### Step 4 - Full Mesh Optimizer (tr_guess => Mesh of {Triangles, Centroids, Affines})

Next we create a multi-triangle mesh within the intersection area. The original scheme Lou had was quite complex. It would surround the point cloud with a boundary-hugging polygon and proceed to carve out triangles from the border inward until the remaining area was exhausted. It was quite error-prone, especially with sparse or discontiguous point clouds. That method lives in file `1_DMesh/CreateMesh` along with the simpler brute force method that supplants it. The new method is discussed in appendix: [Notes on Mesh Creation](ptest_reference.md#notes-on-mesh-creation).

The mesh is sent to the optimizer: `1_DMesh/ImproveMesh::ImproveMesh()`, which does the following:

1. Convert A-cloud to barycentric coordinates.
2. Remember starting control points `orig` and make copy `cpts`.
3. Optimize `cpts` using `0_GEN/Correlation::ImproveControlPts()`.
4. Compare `cpts` to `orig` to sanity-check triangle area changes, if OK...
5. Use `cpts` and `orig` to calculate an affine for each triangle.
6. Compute the optimized triangle centroids.

#### Area Changes

Step (4), area checking, works like this. Loop over the triangles in the mesh and collect the largest individual triangle area change, and compute the total mesh area change.

Next test these against user-supplied maxima:

* `matchparams::TMC` : individual **T**riangle **M**ax **C**hange.
* `matchparams::TSC` : **T**riangle **S**ummed **C**hange.

#### Triangle To Affine

In step (5) the procedure to convert control point changes to an affine is this:

* Let there be a right triangle in global space with vertices { (0,0), (1,0), (0,1) }.
* Let there be a mapping T to local coordinates wherein the vertices are labeled { a, b, c }.
* Affine mapping T can be written directly as:

```
[ [ b.x-a.x c.x-a.x a.x ]
  [ b.y-a.y c.y-a.y a.y ] ]
```

* Construct such a mapping `To` from global space to `orig` vertices.
* Construct another `Tc` from global space to `cpts` vertices.
* Then the affine mapping from `orig` to `cpts` is `[Tc] X [To-inverse]`.

#### Affine Index Raster

You will notice that `dmeshdriver::main()` expects `CalcTransforms()` to compute a 16-bit image called `rmap`. The pixels in this image are either zero, if not mapped from A to B, or are `10 + index of best affine to map with`. You can view this image (basically a depiction of the triangles) by using ptestx and setting `matchparams::GBL.mch.WMT=Y`. The image is then written to `.../temp/za/ia/zb.ib.map.tif`.

The rmap image is also used to create RGB-composite `comp.png` and B-coverage `registered.png` (see [Output (Exhaustive List)](ptest_reference.md#output-exhaustive-list)).

### Step 5 - QA Metrics

After the optimizer returns (in `1_DMesh/RegionToRegionMap::RegionToRegionMap()`) we calculate a variety of quality scores for the resulting A-B match. Most metrics are measured (and logged) on the mesh as a whole (labeled "DEF"), and on the individual triangles (labeled "TRI"):

* Simple mean square difference =  **1/N X SUM<sub>i=0,N</sub>( A<sub>i</sub> - B<sub>i</sub> )<sup>^2</sup>**.
* **_Lou's clever approximation to an_** [approximate earth mover's metric](http://ttic.uchicago.edu/~ssameer/Research/Papers/WEMD_CVPR08.pdf).
	* Require per-triangle EMM < `matchparams::GBL.mch.EMT`.
	* Require whole-mesh EMM < `1.10 X matchparams::GBL.mch.EMT`.
* Fraction of "yellow" overlap pixels = **1/N X (count A-B values within 25%)**.
	* Require yellow >= `matchparams::GBL.mch.FYL`.
* Self-consistency: Calculate maximum spread over triangles of {angle, affine components}.
	* Require spread in effective angle <= `matchparams::GBL.mch.LDA` radians.
	* Require spread in matrix-like components <= `matchparams::GBL.mch.LDR` radians.
	* Require spread in translation components <= `matchparams::GBL.mch.LDC` pixels.
* Close to Original: Average mesh affines, calculate difference from tr_guess.
	* Require difference in matrix-like components <= `2 X matchparams::GBL.mch.LDR` radians.
	* Require difference in translation components <= `2 X matchparams::GBL.mch.LDC` pixels.
* Close to User translation: User provides command-line expected translations `-XYexp=dx,dy`.
	* Require difference of **_any_** triangle translation <= `matchparams::GBL.mch.DXY`.

## <a name="general-source-code-anatomy"></a>General Source Code Anatomy

Most of my sources have the following internal organization, reading from top to bottom:

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

I try to keep functions pretty simple, main() especially. You should be able to glean module workflow by reading the bottom-most functions.

I almost always order a function's parameters like this:

```
returncode = foo( outputs, inputs );
```

When the returned value is a success indicator, I use the convention `0=failure, !0=success`.

Most main() functions begin with a call to **C_myglobals_::**`SetCmdLine()` which will typically do the following:

* Set parameter defaults.
* Open an output (*.log) file named after the application.
* Parse and check command line parameters...
* Turn those into C_myglobals_ member data.

## <a name="input-parameter-handling"></a>Input Parameter Handling

Find ptest `main()` in file: `dmeshdriver.cpp`. The first task is to call the parameter parser: `CGBL_dmesh::SetCmdLine()`.

Generally, the default behavior is:

* Default parameter values come from `matchparams.txt`
* Images and foldmasks are fetched from from `IDB`

These are supplemented and/or selectively overridden by a few command line options.

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

## <a name="output-exhaustive-list"></a>Output (Exhaustive List)</a>

#### Internal

The main() function for ptest first prepares the input data {parameters, images, foldmasks}, then calls the workhorse function `PipelineDeformableMap()` which itself writes the principal pts.same(down) files, and, it returns to the caller several items for diagnostic report generation. These internally consumed "outputs" are:

* func-param : `int NTrans` : count of ALL 'triangles', hence, a-to-b affines found
* func-param : `double* tr_array` : the NTrans affines (6 doubles each)
* func-param : `uint16* map_mask` : A-image where each pixel indexes which listed affine to use

#### External

Here's the roster of all externally visible ptest outputs:

**Always**

* file : `.../temp/za/S(D)x_y/pair_za.ia^zb.ib.log` : crucial status and diagnostic data
* file : `.../temp/za/S(D)x-y/ThmPair_za^zb.txt` : tabulates results of thumbnail FFT matching
* file : `.../temp/za/S(D)x-y/pts.same(down)` : centroids of found triangles (A and B-coords)

**Ptestx Diagnostic Options**

* file : `.../temp/za/ia/zb.ib.map.tif` : rendered map_mask (iff `matchparams::WMT=Y`)
* file : `.../temp/za/ia/zb.ib.tf.txt` : the NTrans affines (iff `matchparams::WTT=Y`)
* ---
* file : `.../temp/za/S(D)x-y/comp.png` : RGB composite A on B result (iff `-v` option)
* file : `.../temp/za/S(D)x-y/registered.png` : accumulate all A mappings to any B (iff `-v` option)
* file : `.../temp/za/S(D)x-y/qual.tif` : intensity encodes A onto B match quality (iff `-v -heatmap` options)
* ---
* file : `.../temp/za/S(D)x-y/thmA_i.tif` : input raster to A-FFT (iff `-dbgcor` option)
* file : `.../temp/za/S(D)x-y/thmB_i.tif` : input raster to B-FFT (iff `-dbgcor` option)
* file : `.../temp/za/S(D)x-y/corr_A_i.tif` : correlation peak mask (iff `-dbgcor` option)
* file : `.../temp/za/S(D)x-y/corr_R_i.tif` : correlation Pearson-R (iff `-dbgcor` option)
* file : `.../temp/za/S(D)x-y/corr_S_i.tif` : correlation Wetzel spectrally modified (iff `-dbgcor` option)

## <a name="image-loading-and-preprocessing"></a>Image Loading and Preprocessing

### Your Own Preprocessing

Of course, this is content-based matching and the pipeline has a few tricks available to improve the signal to background but the results will be better still if you use the best images you can and fix what you can in advance.

* **Software Lens Correction**: If you know there are repeatable geometric distortions in the images (due to the imaging optics, say) you might try this tool: [Distortion_Correction](http://fiji.sc/Distortion_Correction). This is especially useful if the distortion is nonlinear.

* **Sparse Content**: This is primarily a problem in fluorescence data where the features or the labeling may be sparse. Usually you'll have several color channels so obviously you want to align on the densest of them. In some cases it may be useful to superpose multiple channels to increase feature density. Finally, you might apply a **histogram equalization filter** to boost low intensity levels.

### CPixPair::Load()

After `dmeshdriver::main()` calls SetCmdLine() to sort out parameters it creates a `0_GEN/CPixPair` object 'px' and calls `px.Load()`. This loads and optionally filters both A and B. We'll point out notable constraints and highlights here.

#### Dimensions

Immediately upon loading the A and B rasters we check the dimensions, that is, we require wa == wb and ha == hb (though width and height are independent of each other). **_As a general rule, then, all images in a stack must have the same dimensions_**. The natural full dimensions are stored in the px object as `{wf, hf}`.

We find that reducing the image scale by a factor of two or four doesn't significantly degrade overall accuracy, but saves substantially on memory footprint and compute time. Therefore, we iteratively downsample the images (divide each dimension by 2) until each dimension is <= `MAX1DPIX = 2048`. Most subsequent work will be done at this scaled down size `{ws, hs}` with integral scale factor `scl = wf/ws`. Therefore, another data requirement is that **_each dimension must be divisible by two or four (if larger than 2048)_**.

After main() calls px.Load() the resulting scale factor px.scl is used to adjust matchparams values having dimensionality. Matchparams file values are at the full image size. You may notice that LIMXY is absent here. That does get appropriate scale adjustment at the point it's used.

#### Background Subtraction

After raster loading comes intensity flattening which we do by projecting the intensity onto low order Legendre polynomials (independently for X- & Y-axes) to compute those coefficients (one pass) and then subtracting out those components (second pass). The zero order poly is a constant, so models a detector pedestal value. The first order is a line so these can model stage tilt. The second order is a half-cosine which models intensity fall-off toward the edges (vignetting). Higher orders sometimes seem promising but have dubious physical interpretation.

* Orders up to 6 are supported.
* Select highest order with `matchparams::PXBRO`.

#### Resin Masking

For EM images, use `matchparams::PXRESMSK=Y` to enable the generation of mask images that discriminate tissue (mask value 1) from resin (0). Resin is nasty for correlation based feature finding because we always normalize images first. Things that are very dark or very white like dirt are clipped off, yet resin produces middle values. The signal can easily be overwhelmed by resin as it matches to anything pretty well. However, the resmask procedure usefully identifies resin and we can exclude those areas. We make masks with utility function `0_GEN/Maths::ResinMask8()`. The algorithm is entirely heuristic but has proven very effective on EM images from different labs/scopes here at Janelia:

1. Downsample by 8X or so.
2. Apply Sobel edge filter 3 times to fill space with edges.
3. Smooth with a median filter.
4. Threshold at some 8-bit level (100 or so) to make mask.
5. Rescale mask to image dimensions.

There's a Boolean parameter for same-layer Y/N. If same-layer we use a smaller median kernel and a lower threshold. These choices tend to keep a bit more fine detail in the mask which can be very helpful in finding something to match on the periphery of fly brains where there are thin strands of tissue remaining from excising the brain. In the same-layer case, these tiny features actually exist in both images so make genuine matches. If not same-layer tiny features are unlikely to match well so we prefer to smooth them away.

>Note: This function additionally zeroes mask pixels where the image intensity is either zero or saturated at 255, which are unlikely to be real tissue values.

Getting back to px.Load(), the first thing we do is make resin masks with `same-layer = false` in order to get an estimate of how much image area is real tissue (sum all 1-pixels). We simply refuse to match these images if both of them are below 15% tissue. Remember the rule of conservatism in alignment by features: **_It is far worse to include wrong matches than to miss a few good ones...one only needs sufficiently many good points for a solution._**

>Note: If you're reading along in the code, you'll see that if this tissue fraction test fails, we log a message with format `FAIL: Reason`. It is our fairly uniform practice to tag key messages with `FAIL` in all-caps for easy grep searching.

#### Image Pointers

Class CPixPair publishes pointers to raster images, named with tags `aln` or `vfy` according to whether we recommend that image be used for determining alignment (with enhancement filters applied) or for qualifying the result and measuring resulting correlation.

Pointers are used because the scaling and filtering operations are optional or dynamically determined so the same underlying storage may back multiple instances. For example, if scaling isn't needed, then the scaled image and the full-size image are the same thing.

#### Spatial Filtering Experiments

Although I've played repeatedly with increasing correlation using spatial filters like **DoG**, I can't say I was ever impressed with any brand of that. However, CPixPair is a ready home for filtering experiments should you feel so compelled.

You'll also notice some code in here for activating a built-in flavor of 'software lens support'. That was a simple experiment in infrastructure; only taken as far as per-camera affine corrections. Affine correction in particular isn't very potent considering all the other ways of specifying affine correction in ptest (**pre-tweak affines** is another to be covered shortly). The method of lens correction we **_do_** advocate is externally generating new corrected images for alignment before the pipeline ever begins. By this means all stages of the pipeline enjoy correction without needing to reapply it.

In summary, while resin filtering is a really big win, we enable neither magic kernels nor the crude lens 'feature'.

#### Intensity Normalization

Take note that, whatever combination of operations is applied to the pixels, we always call `0_GEN/Maths::Normalize()` afterward. This function first calculates the current mean and standard deviation of a value array and then rescales the array to `mean = 0, stdev = 1`. What we do with pixels for the most part is to pair-up A-pixels with B-pixels and then measure their correlation coefficient (R) which is essentially some flavor of **SUM<sub>i=0,N </sub>( A<sub>i</sub>B<sub>i</sub> )**. The idea behind the normalization is that middle tones are common and boring whereas extreme values (big contrasts) are more rare and interesting. Hence, we down-weight average values by mapping them to zero so their contribution to R is minimal.

>Note: You'll see in `0_GEN/Maths` that we have several alternative ways of 'normalizing' wherein we prefer central values or the tails, and so on, groping for ways to enhance signal-to-background in R sums. This has not been fruitful in our hands.

## <a name="foldmasks-and-connected-regions"></a>Foldmasks and Connected Regions

After `dmeshdriver::main()` loads images, it calls local function `CalcTransforms()` which first loads and prepares the foldmasks, one for each of the two images, and then calls worker function `PipelineDeformableMap()`.

### What's a Foldmask?

A serial section may suffer damage {tears, folds, other} that divide it into two or more distinct subregions, and these subregions may have moved relative to one other by being picked up from the slicer differently or adhering to the substrate differently. Consequently, each subregion would need its own coordinate transform. A foldmask is a metaimage of the same XY dimensions (full scale) as the original. The 8-bit (unsigned) foldmask pixels identify which (montage) subregion the original pixels belong to, using a simple 4-neighbor connectivity rule.

Special foldmask value zero identifies any damage region that should be ignored in all alignment calculations. All pixels having value one are 4-way connected to tissue region one, but have no connection to pixels of value two or any other, and so on for foldmask values two through 255. The rule that, **_distinct subregions within a section are not connected_**, is what allows them to float free of one another and obtain their own solution transforms. Of course when aligning a subregion to an adjacent layer, any pairing of an A-subregion onto a B-subregion is permissible.

There is important missing machinery in our pipeline surrounding foldmasks. Folds, tears and the like are defects that occur at the level of whole layers (montages), that is, the whole montage might be broken into some handful of distinct subregions. Typically a layer is tiled by a dozen to several thousand image tiles. For example, montage subregion three may appear in hundreds of images and it would be extremely valuable if foldmask value three referred to that same connected montage region in all images for this layer. However, although we have an automated tool called `1_Tiny` that does a fair job of calculating a foldmask for an individual image, it has no knowledge of montages. So in each image it just labels each subregion according to the order found. There is no existing mechanism that coordinates the numbering across image tiles. Although we have a degree of foldmask "awareness" in most parts of the pipeline, the ad hoc numbering issue is a huge gap that needs development to exploit foldmasks properly.

These are some of Bill's thoughts on what's needed:

* We can make pretty good montages without fold awareness using the `-nf` option. So that should probably still be done first as we do today. Then we need a tool to resolve the montages into regions and a means to map those masks onto individual image-sized foldmasks.
* By the way, I've tried running `1_Tiny` on montages but tiny really needs full scale and good contrast pixels to work properly.
* The z.id.rgn labels in the pts.same files at this point would be wrong because they have rgn=1. We might want a scheme to edit the pts.same files: assign real labels and remove entries with (argn != brgn) or (argn x brgn = 0). Alternatively, we could do montaging again **with** foldmasks using ptest's current scheme of looping all A-regions against all B-regions, but skipping pairings of bad or unlike labels. Of course previous pts.same files need to be cleared or named differently, which entails an adjustment to the solver's catalog scheme.
* A weakness of redoing montaging with foldmasks is that some tiles may contain only very small portions of subregion 3, and we often fail to get connection points on small patches. That's also an issue for cross-layer ptest jobs. The original tiles are not necessarily optimally placed for tile-tile matching of subregions. Stephan Saalfeld has suggested dynamically drawing virtual tile content from an underlying true montage. Alternatively, one might just paint new true tiles in a manner akin to external lens correction.
* With properly labeled folmasks in hand, we can still use strip-strip alignment to get coarse rigids as we do now, but the block-block matching should be modified so that each block contains purely label-i from one montage and purely label-j from the other. Now, as we find cross-layer block affines we write `make.down` entries for ptest that name all labels `ptest za.ia-ra^zb.ib-rb -Tab=affine`. At the same time ptest needs a mode wherein it processes specified subregions instead of nested loops over regions. Rather, nested region loops move to the block matcher.

### Using Foldmasks to Remove Resin

Returning now to the loading of foldmasks...although labeling connected subregions is imperfect, ptest nevertheless associates a foldmask with every image whether or not the no-folds `-nf` option is invoked. This is handled by `0_GEN/FoldMask::GetFoldMask()`. If using real foldmasks such as those generated by the tiny program, the masks will be loaded from disk files stored in the IDB. If the `-nf` option is specified, we simply allocate a buffer and fill it with ones signifying that the entire image is 4-way connected as a single region.

Now comes a beautiful trick. Every image now has a mask where zero means the pixel is not used and non-zero means it is used. It is entirely consistent with this labeling scheme, to take the logical AND (product) of this mask with any other binary mask wherein zero and one similarly signify not/used. That's exactly how we use the resin masking option to knock out resin pixels. If `matchparams::PXRESMSK=Y` then px.Load() calculates its resmsk members making them non-empty and they are sent as parameters to GetFoldMask().

### <a name="rectangular-cropping-regions"></a>Rectangular Cropping Regions

The same mask intersection method can be used to remove the margins of images; useful because the extreme periphery of images may suffer especially large aberrations. To invoke this mechanism, one creates a cropping file which is a simple text file having a line per referenced camera index with these fields (whitespace separated):

```
cam-id x0 y0 dx dy
```

Remember that in your initial layout metadata file you supply a zero-based camera id number for every image (up to four cameras). Next, before you execute your `dbgo.sht` script to make an IDB, edit the script to include option `-crop=mycroprectfile.txt`.

When dbgo.sht runs it copies your file into the top level of the IDB, giving it the standardized name `crop.txt`. Subsequently, any pipeline component like ptest can check for and implement your rectangles.

### Connected Regions

To recap, foldmasks are **images** whose pixels are region labels [0..255], and this is a very useful and compact means of enumerating subregions and of pairing images with masks. Nevertheless, whole images are not directly useful for match calculations because subregions must be treated independently. Rather, we need to separate the images into their component subregions. Each pairing of a subregion from A against a subregion from B may yield a unique transform.

Subregions may have irregular boundaries and, within an image, may even appear as discontiguous patches (though they ought to be connected somewhere in the containing montage). We can accommodate any shape by describing subregions as two associated (1-1) lists of pixels: one of pixel coordinates (Points) and one of pixel values.

Function `0_GEN/FoldMask::ConnRgnsFromFoldMask()` separates a foldmask into a list of its component **_connected subregions_** with the following structure (see `0_GEN/FoldMask.h`):

```
class ConnRegion {

// Note: All coordinates here are scaled.

public:
    vector<Point>  pts;     // pixels within the region
    double         dx, dy;  // deltas to line up with image2    <== ignore this
    IBox           B;       // bounding box in original image   <== ignore this
    int            id;      // this region's mask value
// other uninteresting members
};
```

Fields {dx, dy, B} are archaic but still serve some older projects; sorry for that. Really, a region is just a set of image points and their common label.

## <a name="approximate-fft-matching"></a>Approximate FFT Matching

### Correlation Calculators

File `0_GEN/Correlation` contains everything that directly calls any FFT library, including filter convolution and all our correlation calculators. The file also contains core parts of the mesh optimizer. There are three main correlation calculators "correlators", each of which take identical parameter lists but internally find peaks using different strategies. The correlators take point clouds as input, which are painted into zeroed rasters for the FFT libraries. The rasters are presized to accommodate all possible displacements between A & B, so the correlation image suffers no aliasing effects. The raster dimensions are increased to the nearest powers of two for algorithm efficiency.

#### R Method, R-Image

From the beginning of my tenure, finding a reliable peak had been an active area of development. The original scheme, embodied by `CorrImagesR()`, calculates Pearson's R, which normalizes the usual cross correlation by the standard deviations of A and B within the overlap region. This method almost guarantees R values between 1 and -1, making the interpretation of R attractively portable across many contexts. Unfortunately, these normalizing factors can get very noisy where displacements are large and overlap area is very small. In consequence, the large displacement periphery of the correlation images is often peppered with high-R tiny peaks. **_Note that we re-map correlation image pixels such that zero displacement (top-left A & B corners coincident) is at the center coordinate of an "R" image._** These large distractor peaks make it tricky to identify the true peak.

#### A-Image

One of the early remedies that is still in force for all calculators is to simply avoid looking for peaks in the edge areas of any correlation image. This is expressed as the requirement that the 1-D and the 2-D overlaps exceed user specified thresholds. Each prospective XY-displacement is tested by user-provided validity callbacks, and the results compose the boolean "A-image"; a mask for the R-image. There are two callbacks:

* `LegalRgn` : caller's choice, but typically, test the 1-D and 2-D overlap bounding box dimensions.
* `LegalCnt` : caller's choice, but typically, test the overlap area's non-zero pixel count.

There is a subtle difference between the two since the point cloud can be sparse in some regions. These callbacks use parameters:

* `matchparams::OLAP1D`
* `matchparams::OLAP2D`.

A later addition to the A-image was correlator parameters {`Ox, Oy, Rx, Ry`} that specify (iff R >= 0) a disk for the expected peak. This disk is logically AND-ed with the edge mask described above.

#### F Method, F-Image

During some early era I put considerable effort on filtering the R-image in the spatial domain, giving rise to `CorrImagesF()`. I tried several heuristic methods to improve reliability. The current implementation applies a 3x3 LoG filter to the R image to sharpen candidate peaks, resulting in the **F-image**. Each candidate in the F-image is surrounded by its own annular zone called its **guard band**. The inner radius is where the F-value has fallen to 50% of that candidate's peak F-value. The outer radius of the band is nominally 3X the inner radius. A candidate is rejected if it is too close to another, in particular, if any other peak is in its guard band AND the encroaching peak is higher than: (`matchparams::NBMAXHT`) X (candidate's own peak height). The idea is that true peaks are well isolated from any rivals. This method can get true peaks where R fails, but the reverse is also true, so neither is a winner. Of course I've tried combining the two, but the best weighting between them is decidedly content-dependent.

The F-type calculator internally creates {F-image, R-image, A-image}. F and A are used to find the peak within the valid area, but the correlation we report is interpolated from the R-image since, again, Pearson R has nicer normalization properties.

#### S Method, S-Image

In April of 2014 Stephan Saalfeld organized an [Alignment Hackathon](https://github.com/saalfeldlab/alignment-hackathon-2014-04/wiki) here at Janelia. I had the good fortune to talk with Art Wetzel (Pittsburgh Supercomputing Center) about his current research interests and learned that he was developing a new approach to alignment that relies heavily on FFT-based correlation over many length scales. I can't divulge much of this unpublished work, but Art convinced me I should try manipulating the power spectra of the FFTs used in correlation to emphasize feature size scales of relevance, hence, "S method." Art's simple trick, employed here, boosts signal to noise from 10% to 100% or more over the corresponding R-image depending upon content. `CorrImagesS()` is now our de facto standard. Unfortunately, the natural value range is very variable and I employ a hack to drive the reported values closer to [-1..1] so that one still specifies thresholds in matchparams and scriptparams files. This method internally generates {S-image, R-image, A-image}.

#### Ptestx

An absolutely invaluable diagnostic and tuning facility is to run wrapper function `ptestx` with command line option `-dbgcor`. The correlator outputs the following diagnostics (in PWD S- or D- folder), and then ptest quits:

* The cropped/filtered thumbnails used in correlation {`thmA.tif, thmB.tif`}.
* Internal {`corr_R.tif, corr_A.tif`} always.
* Internal {`corr_S.tif, corr_F.tif`} according to correlator type.
* `Table of largest peaks` in the pair_xxx.log file.

### Angle Sweep Modes

#### Don't Do This

To really appreciate the legacy matching machinery in here you need a brief history of the development path. Also understand that Lou Scheffer, my esteemed predecessor, had only intended to align just one data set, the result of which would keep the world busy for ten years. Expediency was the rule. It would never be used again by Lou or anybody else... until I asked for it and it was generously provided. I did not come to this code from a background in alignment, so having a great deal of respect for Lou, I always assumed he did it this way for good reason. My changes were slow and incremental and involved repeated testing to be sure I was improving rather than degrading the master's work. Over time I became expert and in hindsight many things are obvious and many are just plain silly, I know that.

The original ptest I got from Lou didn't have any awareness of whether its two images were from the same or different layers. All pairs were subjected to an angle scan from -45 to +45 degrees. At each angle step one ran the correlator and collected {dx, dy, R}. Next one employed an elaborate analysis to determine the peak R of the sweep, hence, the correct relative angle and displacement (dx, dy) for the image pair. An elaborate analysis of the R vs angle curve was needed because the original R correlation method was very unreliable, often returning a false peak and making the angle curve very jagged. Ad hoc rules were needed like "the true peak must be higher than the points on either side, and those must be higher than their neighbors" and the like, hopefully discriminating noise from signal. These never worked very well.

An improvement I had created was to do the sweep and collect {dx, dy, R} as before, but then fit dx-vs-angle and dy-vs-angle to lines (over small domains) and require the linear correlation coefficients to exceed a preset value (0.7). This rule of peak invariance under small rotations worked sometimes but not all the time.

Another scheme I tried was that each S- or D-folder contained a file `ThmPair_za^zb.txt` that would collect sweep results from each pair job {errorcode, angle, R, transform}. **_These files continue to be written regardless of selected angle mode_**. After enough entries existed in the file then the remaining jobs could read the file and use the median angle as a concensus estimator for the central angle of their scan, AND, one could then use a narrower scan range. By the way, all ptest jobs in an S- or D-folder are run from a make file on a single machine and in that condition we can use mutexes to arbitrate file access to ThmPair files (similar idea for pts.same(down) files). Note too that the jobs composing an S- or D-folder come from a block, that is, a common local neighborhood, so it is reasonable that such jobs should produce a shared angle.

At some point a revelation struck. When the images were in the same layer, although there was considerable variance in the reported best angles, it looked like that angle was usually zero. Could we skip angle scans altogether in the same layer? This was very exciting, and at first, experiment worthy but non-obvious. Ultimately this was a big win, as you can imagine. It was the beginning of including labels on the command line and splitting many of the matchparams into same-layer and cross-layer values.

Still, for the cross-layer case, there was only limited success with these ad hoc sweeping techniques. Admittedly, the invention of resin masking came fairly recently and I am certain that angle sweeping would look much better with resin masking enabled, but I have never tested that. Rather, I am going to insist that sweeping two images inside of ptest will always be inferior to our modern preferred approach. The modern way says that if the images are in the same layer then the angle is zero, end of story. For the cross-layer case, the modern way is to sweep much larger blocks of tissue against each other for enhanced signal to noise, and then report the best angle to ptest.

In summary, these are the older `matchparams::MODE` selectors you should **not** use:

* `E` : Approx matcher just returns initial Tab.
* `F` : Just read ThmPair file and fetch existing result for this label.
* `C` : Angle sweep using parameter `CTR` as the sweep center angle.
* `Z` : Angle sweep using zero as center angle.
* `Y` : Angle sweep, center angle dynamically determined as median of ThmPair table.

The above sweep modes primarily set how the center angle is determined. The range of the angle sweep is set by these matchparams:

* `HFANGDN` : half range denovo (degrees), that is, when ThmPair entry count is too low for median.
* `HFANGPR` : half range when prior angles (median) is available.

#### Do This Instead

Ptest should not hunt angles. Rather it will use zero for montages, and it will extract an angle from the block-derived Tab in the cross-layer case.

* `MODE_SL=M` : Force Tab to angle zero.
* `MODE_XL=N` : Extract angle from Tab.

These modes not only specify the angle between the images, they also paint the correlation A-image with a disc specifying where the translation peak is. Tab sets the disc origin. The radius of the disc is specified by `matchparams::LIMXY`.

#### Caveat

There remains a critical failure mode with the disc methods (modes M and N). The promise and the premise of these modes is that ptest does not, itself, have the proper perspective to judge the validity of the peak searching disc. The correlator always reports **_something_** from the disc...right or wrong, there is always a maximal value in there. So, if your inputs are wrong ptest can't tell. It would then use a bad starting guess for the mesh and could quite possibly report bogus points. To avoid that every effort should be made to use sensible values for parameters, here `LIMXY`, and whenever applicable, diligently check intermediate pipeline results because at each stage what comes out is only as good as what goes in.

### <a name="tweaks"></a>Tweaks

#### Pre-tweaks

Using command line parameters you can specify a deformation "Tdfm" to be applied to A before matching it to B, if you already have one in mind. More often however, A and B are slightly different due to dehydration or focus disparities that are not known in advance. There is a built-in mechanism called `Pretweaks()` that dynamically generates deformations. Here's how it's optionally engaged:

* Set `matchparams::PRETWEAK=Y` to enable the feature.
* If using an angle sweep mode
	* Sweep is done normally first
	* If best R is below RTRSH, Pretweaks() is called
	* Sweep is tried again one more time
* If using disc mode
	* Pretweaks() is called unconditionally
	* Disc method applied

To see how Pretweaks() works let's start with the correlator which acts on two point clouds. We never call the correlator directly, instead we use a helper class `0_GEN/CThmScan`. In particular, method `RFromAngle()` transforms the A-cloud and then calls the correlator. The transform is the affine composition `R X Tdfm X Tptwk`, where R is a pure rotation matrix made from RFromAngle's `deg parameter`, Tdfm is a class member initialized from the user's Tdfm, and Tptwk is another member initially set to unity, but modifiable via the Pretweaks() method.

Pretweaks starts with a list of all five affine "deformations" {total scale, X-scale, Y-scale, X-skew, Y-skew}, each one of which is describable by a scalar parameter. At the given constant center angle (passed to RFromAngle()) it applies each type of deformation using eleven parameter values centered around one (for scale) or around zero (for skew). It tabulates {type, parameter, R} for all 55 possibilities. If the best of these (highest R) is better than before, that deformation multiplies the existing Tptwk, and the winning type is struck from the type list. Next the procedure is repeated with the list of four remaining types, and so on. Any arbitrary small deformation can be composed this way.

#### Post-tweaks

Another parameter, `matchparams::TWEAKS`, activates the same mechanism, but this time, Pretweaks() is called after a best angle and translation are found. The idea is to fine tune the final affine before sending it on to the mesh calculator. Pretweaks() is fairly time consuming and the mesh apparatus subsumes all affine transformations, so I recommend keeping this off until proven necessary.

## <a name="notes-on-mesh-creation"></a>Notes on Mesh Creation

My current mesh building method starts by erecting a regular grid of rectangles over the intersection bounding box. The initial grid cell size is given by your "mesh nominal length" `matchparams::MNL`. **_This value is used for same- and cross-layer cases._** We then adjust the dimensions down slightly so the rectangle mesh roughly circumscribes the bounding box. After one more step, these rectangles will each be split along the diagonal to become two triangles, but we have one more sizing step first.

**_In the cross-layer case_** we consult your "minimum triangle area" parameter `matchparams::MTA`. We gradually reduce the triangle count while the average mesh triangle area remains lower than MTA. Using MTA only for cross-layer cases allows you some freedom to separately size the two cases.

Finally, to improve calculation efficiency, we get a count of points falling in each triangle, which I term the occupancy. All triangles having occupancy lower than 30% of their area are culled from the mesh, and any control points that become unreferenced are removed.

You can see key results from mesh building in the job's pair_xxx.log file. That paragraph looks something like this:

```
---- Building mesh - deformable ----
Lx Dx Nx:  1649    74.00  22
Ly Dy Ny:  1649    74.00  22

STAT: From 2722500 pts, got 968 triangles, 529 control points.
Timer: MeshCreate took 10.380 seconds.
```

* Lx, Ly are the dimensions of the (**scaled by px.scl**) intersection bounding box.
* Dx, Dy are the scaled rectangular grid cell dimensions.
* Nx, Ny are the counts of rectangular grid cells in each dimension.

### Tips on Mesh Sizing

In the cross-layer example above the user had encountered memory paging problems and slow performance due to suboptimal choices for the mesh sizing parameters. In the solver, we're only going to get one transform per image tile, so the 968 triangles (that is, point pairs) being generated here is excessive. It's instructive to look at this case in detail.

The matchparams were set to defaults from an unrelated data set, remember, these are specified in full scale pixels:

* MNL=300
* MTA=10000

The images were each 6600 X 6600 at full scale, and in this block-face example cross-layer images overlap each other almost completely. Let's first check that the mesh algorithm worked as advertised. Since 6600 is bigger than the linear size limit `MAX1DPIX = 2048` in CPixPair::Load(), each dimension would be halved twice, down to 1650 X 1650 which is consistent with the {Lx, Ly} reported above. Next, 6600 / MNL = 6600 / 300 = 22, which is the the number of grid cells created in each dimension. Finally, the average triangle in this mesh is 74 X 74 / 2 = 2738, which is greater than MTA/(scl X scl) = 10000/(4 X 4) = 625.

Here's what one should do. Ten point pairs along an edge in the same-layer case is adequate in most cases but these are huge images so I might go up to 16 points per edge. So in each dimension we only need 16 triangles, or 8 rectangles, so we should set MNL = 6600 / 8 = 825.

To adjust sizing for the cross-layer cases, note that the solver does not apply any explicit weighting factor for same-layer vs cross-layer point-pairs. Each point-pair constraint equation gets equal weight. If one wanted to apply a weight to a constraint equation, one would do that by multiplying all terms of that equation by a weight. Implicitly then, the same-layer and cross-layer contributions are differentially weighted according to how many constraints of each type we use. We can write the ratio of same- to cross-layer connections in a succinct way like this.

```
E = average number of same-layer edge connections.
F = average number of face connections to one cross layer.
R = desired balance of same/cross connections.

R = (4 X E) / (2 X F), so,

F = 2 X E / R.
```

For block-face imaging, the cross-layer matching is very good and we may be inclined to set R = 1. Then in this example we want each cross-layer pair to use 2 X 16 = 32 triangles, so we set MTA = 6600 X 6600 / 32 = 1361250, but we can make that a round number like 1200000. You are free to decide R, hence, MTA according to your beliefs about match quality.

## <a name="ptestx-environment-try-folder"></a>Ptestx Environment (Try Folder)

Ptestx is a very simple helper tool for ptest. All it really does is create a `temp` folder structure, copy two needed input files there (matchparams and imageparams), change directory into the S- or D-folder and then call ptest with the remaining command line parameters. The parameters pertaining specifically to ptestx are:

* `-idb=path` : path to your idb folder.
* `-d=temp` : ptestx creates local "temp" folder with this name.
* `-clr` : delete existing local temp folder first.
* `-prm=path` : copy this matchparams.txt file into local temp.

I regularly use ptestx to tune parameters in matchparams before committing batch jobs to the cluster. Here's what I do. Following the general workflow prescribed in [**user_quick_reference**](user_quick_reference.md):

### Set Up a Try Folder

After **"Create Image Database and 'temp0' Work Folder"** and before "**Extract Same-Layer Point-Pairs**" create a folder called "try" here:

* prms
	* matchparams.txt
	* scriptparams.txt
* mylayout.txt
* topgo.sht
* idb <= exists after running dbgo.sht
* temp0 <= exists after running mongo.sht
* try
	* pgo.sht <= copy here from `1_Ptestx`
	* matchparams.txt <= copy here from prms folder

As shown, copy the example script `1_Ptestx/pgo.sht` into your try folder. Also copy the `matchparams.txt` file from your starting prms folder into the try folder so you can play with a local copy without disturbing "the real thing".

The strategy is to edit the local `try/matchparams.txt` and apply it to several image pairs until you like the results, and only then make it the working copy. Note that, at this point, the true working copy exists in two places, so you'll want to update both with your changes:

* Copy back to prms folder for posterity (master copy for this project).
* Copy to temp0 folder because mongo.sht places the **_real_** working copy there.

### Test Same-layer Case

The reason I had you run both dbgo.sht and mongo.sht first is that these create the infrastructure that makes it super easy to call ptestx. In particular, mongo.sht determines which images likely overlap in the montages and then writes make.same files listing these pairs. So, open one of the make.same files and copy the label data as indicated by the [[]] brackets in this snippet (grab the -nf flag too):

```
... etc ...
4/1.10.map.tif:
	ptest 1.4^1.10 -nf ${EXTRA}

4/1.8.map.tif:
	ptest [[1.4^1.8 -nf]] ${EXTRA}

3/1.9.map.tif:
	ptest 1.3^1.9 -nf ${EXTRA}
... etc ...
```

Now edit your local pgo.sht to read like this:

```
ptestx 1.4^1.8 -nf -clr -d=temp -prm=matchparams.txt -v
```

When you run pgo.sht, you'll get a new temp folder inside your try folder. Dig down into `try/temp/1/S0_0` to find the pair_xxx.log file. If it ran to completion the S0_0 folder will also contain pts.same and comp.tif showing A in green, B in red and yellow where they overlap. If the job failed to complete, the log will indicate why (look for a line labeled by "FAIL:"). If the issue appears to be in the initial FFT stage, edit the pgo.sht script like this:

```
ptestx 1.4^1.8 -nf -clr -d=temp -prm=matchparams.txt -CTR=0 -dbgcor
```

Now you'll get the thumbnail images the correlator uses (see [Output (Exhaustive List)](ptest_reference.md#output-exhaustive-list)).

To make a change, edit your local copy of matchparams: `try/matchparams.txt` and run pgo.sht again.

At this stage you can only tune parameters for montaging because we don't yet know which tiles pair up in the cross-layer case. Nevertheless, after copying any changes back to the working matchparams you are ready to resume the main workflow by changing directory into temp0 and running `ssub.sht i j`.

### Test Cross_layer Case

Returning to the main workflow, after **"Run Block-Block Alignment"**, the D- folders are created and each contains a make.down script which tells us which tiles pair together and the recommended -Tab parameter between them. So, before submitting `dsub.sht i j` to the cluster, we want to try some cross-layer cases to make sure the "_XL" matchparams make sense. This time we open one or more make.down files and copy out a test case (as indicated by [[]]):

```
... etc ...
12/0.2.map.tif:
	ptest 1.12^0.2 -Tab=1.003811,-0.027988,-1280.633680,0.023455,1.007369,-681.869390 -nf ${EXTRA}

12/0.1.map.tif:
	ptest [[1.12^0.1 -Tab=1.004147,-0.025540,-1279.945138,0.023444,1.007651,176.454982 -nf]] ${EXTRA}

11/0.11.map.tif:
	ptest 1.11^0.11 -Tab=1.003358,-0.027931,-155.192769,0.023796,1.006780,177.440927 -nf ${EXTRA}
... etc ...
```

Edit pgo.sht to look like this:

```
ptestx 1.12^0.1 -Tab=1.004147,-0.025540,-1279.945138,0.023444,1.007651,176.454982 -nf -clr -d=temp -prm=matchparams.txt -v
```

Carry on editing local matchparams and trying image pairs until confident. Then copy the matchparams back to the prms folder, and back to the temp0 folder. Resume workflow by changing directory into temp0 and running `dsub.sht i j`.

_fin_

