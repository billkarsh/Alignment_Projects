# Aligner Methods Overview

### Basic Concepts: Regions & Connecting Points

Our world is composed of two basic building-blocks:

* 2-D __regions__ (each having an affine)
* __point-pairs__ (connecting the regions)

Regions tile (cover) a space and they must have non-empty overlaps (intersections) with other regions. What the aligner fundamentally does is the following two step process:

1. Match putatively identical content between pairs of images, and tabulate matched point-pairs (aka correspondence points, features, or just 'points' for brevity).
2. Use the point-pairs to solve for 'best-fitting' linear mappings from local image coordinates to a shared common global space (a montage or stack).

Regions are labeled to tell them apart and to link them to meta-data such as which points belong to them. Many labeling schemes are possible; ours however, uses a physically motivated triple of index numbers:

`z.id-rgn`

* __z__: 0-based positive integer labeling a section (aka layer, montage).
* __id__: 0-based positive integer __*unique within a layer*__, labeling an image (aka tile).
* __rgn__: 1-based positive integer labeling a subregion within an image.

We carry the rgn-index to allow images to be subdivided by folds, tears and the like into subregions that can get their own transforms. Rgn index __zero__ signifies a fold or other masked area to be ignored. Otherwise the rgn-index must be in range [1..255].

Note that each image can be paired with an 8-bit mask image of equal dimensions called a **'foldmask'** whose pixels are the rgn-indices just described.

>Note too, that effective management of foldmasks and rgn differentiation is only about half finished in the present version and is not recommended for use on critical data. Rather, choose the `-nf` (**_no folds_**) option wherever applicable in parameter files and command lines. The `-nf` mode does work as advertised.

#### Points File Content

As an alignment progresses it generates files of type `pts.same` and `pts.down` with the tabulated points-pairs. Here are some lines from a pts.down file, which collects matches from different layers.

```
CPOINT2 1.12-1 548.347679 716.475949 0.11-1 392.363415 37.427421
CPOINT2 1.12-1 906.183423 788.011835 0.11-1 749.536275 120.037851
CPOINT2 1.12-1 531.969897 905.889984 0.11-1 370.319378 227.175693
...
```

For example, line 1 states point (x,y)~(548,716) in image `A` with z.id-rgn = 1.12-1 matches point (392,37) in image `B` 0.11-1.

There is always an A-image, listed first, being matched to a B-image. The coordinates are local to the respective images. Global coordinates will emerge at the end from the final solve.

Note too that between a given image pair there are multiple point-pairs reported so that we can construct non-degenerate linear systems to solve for 6-parameter affines, and if desired, 8-parameter homographies.

#### LSQ Solver

The input data needed to solve for the desired affines are these:

1. A specification of which regions to solve for; usually a whole number of layers, so specified as a range `z=i,j` of layer indices. These are the dependent variables.
2. Starting guesses for each of the transforms, which will be iteratively improved. We'll discuss where appropriate guesses come from in a later section.
3. The connecting point-pairs, which are treated as constants.

The solving algorithm is very simple. It looks like this:

##### For each 'iteration' of the solver...

**_The Data_**

- We have a list `Xsrc` of input transform elements from the previous iteration (or the starting guesses if pass zero). These remain constant through the current iteration.
- We have a list `Xdst` of output elements that we will calculate on this iteration.
- For every transform/region we have lists of all those point-pairs that have either an A- or B-member in that region. The point-pairs are stored in structures like this:

```
struct pointpair {
    int   a, b;      // which regions/transforms
    Point pa, pb;    // matched coordinate pairs
}
```

**_The Algorithm_**

For each of the `Xdst` destination transforms (one at a time)...

1. Loop over the point-pairs belonging to this transform, which we will call the a-member.
2. For each point-pair add an equation/row to an overdetermined matrix `[A]` and column vector `[B]` of the form: `Ta(pa) = Tb(pb)`. This expresses that the same content should map to the same place in the global space.
3. After adding all equations, create the corresponding normal equations: `[At][A][X]=[At][B]`, where [At] signifies A-transpose.
4. Solve for X which are the elements of this `Xdst`.
5. Do (1..4) for each of the `Xdst`.

At the end of each iteration, swap (pointers to) `Xsrc` <=> `Xdst`.

By the way, all of that can be (is) done in a massively parallel way. Solving is usually the quickest part of the process.

Finally, the last Xdst are written to disk as 'the solution', along with some metrics about the residual error distribution, how 'regular' the solution looks, and if we had to kick out image tiles because their transforms became implausibly distorted.

### Systematic Workflow to Make Point-Pairs

To feed the end-game process, just described, that derives the final global transforms from points that link regions (images) to each other we employ the following key devices:

1. A pipeline component (`ptest.exe` built by project `1_DMesh`) that takes three principal **inputs={images A, B, and approximate transform Tab from A to B}** and produces **output=tabulated points**.

2. A workflow (a sequence of bash scripts and binaries) designed to efficiently determine which images should be paired with each other (along with estimated Tab guesses) so that we can feed these pair-wise ptest jobs to a compute cluster in a parallelized way.

>The ptest program is such an important and complex workhorse that we devote a separate reference to it: [**ptest_reference**](ptest_reference.md). For the remainder of the workflow discussion let us just assume a magical black box taking image pairs to points tables.

#### Main Workflow Steps

Thus far we introduced a two-step process cartoon--

1. Create points.
2. Solve all the transforms.

We now replace that with a more accurate but still cryptic cartoon--

1. Create image meta-database **'IDB'**.
2. Create organized workspace **'temp'** and required bash scripts.
3. Pair **same-layer** images using stage coords and extract `pts.same` files.
4. Use `pts.same` files to solve for independent layer montages.
5. Create **'cross_wkspc'** and scripts to drive cross-layer work.
6. Render montages at low resolution and coarsely align using selected strips.
7. User reviews, fixes and blesses resulting low resolution stack.
8. Coarse stack re-expanded and saved as a **'scaffold'** of starting guesses for final solve.
9. Subdivide layers into grids of rendered blocks and refine alignment.
10. Pair cross-layer images using block affines and extract `pts.down` files.
11. Solve final stack denovo using all `pts.same` and `pts.down` files.

For concreteness, the actual sequence of commands is given here: [**user_quick_reference**](user_quick_reference.md).

Each step is manually initiated by a designated bash script. Following each step one runs a designated report-generating bash script and reviews the recommended sanity checking data. Each step thereby proceeds from a state of reasonable confidence, so there is a reasonable chance that a successful outcome will emerge, even in a system with many millions of images.

> Note: Our bash scripts have file type `*.sht` which stands for **shell text**. I map this file type to a text editor so when I double-click on it (say, in Konqueror) I can review what I am about to do and make last minute parameter changes. I spend more time reading scripts than executing them.

The main workflow steps get more coverage in the following sections.

#### IDB (A Structured Repository for Raw Data)

An alignment session begins with a user provided file of information about the images. It can be a text file or TrakEM2 XML file. A text layout line has the following fields:

>Z tileID a00 a01 a02 a10 a11 a12 col row cam full_path

* Z is 0-based layer index.
	* All the images for a layer should be grouped together in a layout file.
	* Listed Z's should never decrease.
	* Z's need not be contiguous: we will automatically match across gaps.
* TileID is 0-based and is unique within its layer.
	* `Tip: ID's can be anything, but assigning, say, col*1000+row, makes navigating easier`.
* Components of affine with 2D order like this:
	* `.[ [a00 a01 a02]`
	* `.. [a10 a11 a12] ]`.
	* Equiv 1D labeling like this: [a0 a1 a2 a3 a4 a5] = [a00 a01 a02 a10 a11 a12].
	* Most often, [a0 a1 a3 a4] = identity [1 0 0 1] and [a2 a5] = image [left top].
* col and row are 0-based col and row indices for human debugging.
	* If unavailable, these should be set each to -999.
* cam is 0-based camera identifier for debugging.
	* In a one-camera setup, use cam = 0.
* Full absolute path to image.
	* Images may be accessed from a variety of working directories so must be absolute.
	* Images can be TIF (8,16) PNG (8,16) MRC(16).

Script `dbgo.sht (project 1_MakeIDB)` reformats the user layout into a new IDB folder/file organization for efficient access by pipeline stages. The IDB device also has these features:

* Scripts/destinations for autogenerating 8-bit png images from mrc files.
* Scripts/destinations for autogenerating associated foldmasks.
* Designated location for foldmask metadata (counts of rgns in each image).
* _`Project 1_Tiny` builds pipeline component for mrc and foldmask operations._
* File `0_GEN/PipeFiles.cpp` provides API to access these data using {z,id} key-pair.

#### Workspace 'temp' (Universal Filing Scheme for Machines and Humans)

The aligner produces copious output: scripts, tables, logs and so on. To make navigating these things easier for everybody we have created an organizational scheme that is assumed throughout the pipeline:

* Top-level project folder/
	* topmost scripts.
	* prms/ (alignment parameter files)
	* idb/
	* temp0/ (version 0 workspace)
		* copies of current matchparams.txt and scriptparams.txt
		* same-layer scripts
		* down scripts
		* 0/ (everything for layer 0)
			* // ------------------------------------------------
			* // Iff `scriptparams::createauxdirs=Y`
			* 0/ (auxdir for debug info from matching z.id=0.0 onto others)
			* 1/ (auxdir for debug info from matching z.id=0.1 onto others)
			* etc.
			* // ------------------------------------------------
			* S0_0/ (workspace for same-layer block (x,y)=(0,0))
			* S0_1/ (ditto for block (0,1))
			* etc.
			* montage/ (workspace for solving layer 0 montage)
			* D0_0/ (workspace for down-layer block (x,y)=(0,0))
			* D0_1/ (ditto for block (0,1))
			* etc.
		* 1/ (ditto for layer 1)
		* n/ (ditto for layer n)
		* cross_wkspc/ (strip and block intermediate alignment)
		* stack/ (workspace for solving final stack)

The folders become filled with logs and output data corresponding to the pipeline stage currently running.

Workspace folder name, e.g. 'temp0', is just a local tradition and you can set any name you want in script `mongo.sht (project 1_MakeMontages)`.

#### Same-Layer Work (Montaging)

Script `mongo.sht (project 1_MakeMontages)` not only creates a new 'temp0' workspace, it also uses the stage coordinates embedded in your initial layout file to determine which image pairs have sufficient overlap to list in a make.same script.

Internally, it actually creates a single long list of such pairs for each whole montage. Then, the pairs are partitioned into a grid of blocks according to `scriptparams::montageblocksize` and the top-left corner of each pair's A-member. An Sx_y subfolder is created for each grid cell along with a make.same file for the pair jobs in that block. It's just a divide and conquer scheme for parallelizing the work load. The assignment of jobs to this or that grid cell is not important. In fact, if any grid cells look underpopulated we reassign their members to neighbor cells to getter a better balance.

Typical make.same 'make rules' look like this:

```
4/1.8.map.tif:
	ptest 1.4^1.8 -nf ${EXTRA}

3/1.9.map.tif:
	ptest 1.3^1.9 -nf ${EXTRA}
```

Note that makefiles like this will indeed be given to the make command and in the past we may have actually used features of the make system such as checking the build date of a target like '4/1.8.map.tif'. However, these days, all we really want to do with the make facility is simply list the ptest jobs to be executed, and use make's -j parameter to say for example: `make -f make.same -j 8`, that is, on a given cluster node, run this list of jobs, running no more than eight jobs at a time. Again, it's a device for load balancing and resource management; nothing more complicated than that.

What you should note about the rules (jobs) in a make.same file is that there is a requisite pair of images (always with same layer index) and the only other parameter is -nf signifying that foldmasks are not being used. Conspicuously absent is any parameter of form [-Tab=a0,a1,a2,a3,a4,a5] to give a starting hint for A=>B. Rather, the absence of that hint indicates that ptest should look up the initial affines from the IDB and create Tab=[Tb-inv][Ta].

Montaging is highly reliable. If a pair of images in the same layer overlap and if that overlap contains any kind of features, then those features should actually match in both pictures. It doesn't matter if the features themselves are excellent quality tissue or are dirt on the sample; we still have two pictures of the same thing. Moreover, we can montage any layer in isolation without reference to other layers.

We make all the montages ahead of any cross-layer work so that we can render them (and parts thereof) to deduce how the layers fit together. Once we understand how the layers match up, we will be able to run make.down jobs and get high quality pts.down data. **_In the end, the final solve uses all the pts.same and pts.down data to solve for all the transforms denovo_**. This is often a point of confusion, but the montage transforms only assist with rough cross-layer positioning. **_Montage transforms are not used in the final solve_**. No claim is made that this is fundamentally necessary or optimal. It is simply an accurate description of how the present method works.

#### Cross-Layer Work (Strip Alignment)

Cross-layer work has two phases:

1. **Strip-strip**: to get a coarse global alignment and a solving scaffold.
2. **Block-block**: to improve the former, and **write make.down files** (tile-pair lists).

>Note: All of the work described in this section is done at reduced resolution, usually between 4X and 50X or so. This emphasizes larger features like blood vessels that change slowly across sections. It also allows cheaper and faster computation.

**_In the rest of this section we'll discuss the strip-strip phase._**

Very often in serial section work there are layers that got turned around at some point in the slice pick-up or grid positioning processes. A key benefit of strip-strip alignment is to unburden the user from (nearly all) manual prealignment of layers. Here's an outline of the steps:

1. For every pair of adjacent montages (A onto B below).
2. Use the `1_Scapeops` component to:
	* Paint all of B.
	* Paint a horizontal strip through the middle of A.
	* Paint a vertical strip through the middle of B.
	* Make a wide angle sweep of strip A against B and log the best rigid transform.
3. Now use the rigids to compose a low resolution stack with 1 image/layer (as TraKEM2).
4. Have the user open, inspect and amend that result.
5. Re-expand the saved result, transferring the rigids to all tiles => scaffold.

Of course, for the top-most layer, we only tell `scapeops.exe` to make a montage for layer B (there is no layer A above).

We use strips for the angle sweeps to save computation time. Crossed strips have the potential to find overlap over a wide range of translation and rotation. On occasion this fails to find a valid transform:

* Damaged or missing tissue right smack in the intersecting part of the strip.
* Irregular montage shape such that strips simply do not cross.
* Extremely low contrast data with very low signal/noise.
* Irregular image intensity patterns making strong artifactual features.

A real strength of this method is that a few misses do not matter because the user will visually check all the results and can quickly fix any errant layers using the `Align/Align Using Manual Landmarks` tool in TrakEM2. In an hour one can review around a thousand layers. Again, this is all done at low resolution with one flat image per layer, so the review process is very manageable.

> Note: Blockface methods still benefit from software tweaking, but in that case you wouldn't want to waste time on a wide angle sweep looking for a strip angle you know is essentially zero. We've made the strip angle search parameters variables. In this case you may want to set `scriptparams::stripsweepspan=0` degrees and set `scriptparams::stripsweepstep=1` degree. This will skip the angle sweep and simply calculate the best translational components for zero degrees. The block-block refinement will further tweak that result in the usual way.

#### Cross-Layer Work (Block Alignment)

Whereas strip alignment roughly got the layers in the right place, what we want to get out of block-block alignment is sufficient improvement of the alignment in a grid of smallish neighborhoods to establish which tiles really pair with each other, and moreover, what their estimated Tab transforms are (this will assist ptest in selecting the correct peak from their correlation image).

You may ask why, after getting a strip alignment, we need to refine that any further. The answer is that for really big montages, that is, with linear dimensions in the hundreds of tiles, although strip alignment gets their centers pretty close, there is a large lever arm out to the corners. We could easily be off by a whole tile. That only gets worse if the layers are not identical in shape, as may happen from folds, tears, or just plain differential dehydration from beam exposure or slice thickness variation.

To get a good alignment result each tile needs to be well connected to the whole, that is, each tile needs an _adequate_ count of match points in all directions. Sometimes, at some point in the sample handling and/or imaging pipeline, images can go missing...sometime big chunks may go missing. That means that sometimes, there is either an insufficient quantity of tissue in the adjacent (below) layer or the tissue is of marginal/poor quality.

The block matcher `1_Cross_ThisBlock` attempts to match each block to successively lower (farther away) tissue until it gets adequate coverage of the block with reasonable quality tissue. Each block may match piecewise to several other layers below. You can control several behaviors of the block matcher using `scriptparams.txt` parameters:

* The farthest away one should look.
* The threshold for adequate coverage.
* The threshold for minimally useable correlation.
* The threshold for _nominal_ correlation, meaning "If it's this good we're done."

As the block matcher continues to look at farther tissue it keeps everything it has already seen sorted according to correlation value. When determining if it has sampled enough to meet the coverage threshold it considers the highest quality tissue first, even if farther away. An exception is that the immediately adjacent layer has priority as long as its quality is at least 90% that of the overall best. In our experience, using only the best quality matches, without _requiring_ links to adjacent tissue, produces too much jitter.

Again, each block-block job is running in its corresponding `Dx_y` subfolder and should produce two key types of output file there:

* `make.down` file listing ptest image pairs with Tab guesses.
* One or more files named `ThmPair_Za^Zb.txt`.

The ThmPair files are tables that get filled by ptest (when you later run `dsub.sht`). Each line shows what ptest determined as the best affine between a given A-tile/rgn against a B-tile/rgn, along with an error code, a correlation and gross angle. These are results from the FFT matching of thumbnails. Their original purpose is a bit archaic: it allowed all of the ptest jobs running in a given block (on one machine) to share their results with one another and thereby reach a consensus about the likely block angle. That's now superceded by the block matcher.

Nevertheless, these ThmPair files serve a new important function. The block matcher creates one of these files for each new layer that this block connects to and the layer indices compose the names of these files. This is the clue that the solver will use to determine, for every layer, the other layers to which it connects. That catalog of layerwise connectivity helps the solver to partition a large stack onto multiple machines so that each machine loads point data for an inner range zi of layers it is solving for, and, for an outer range zo of layers that those depend upon.

After the block-block phase completes (and is sanity checked) the user returns to the top level of the temp folder and runs the script `dsub.sht` which crawls into all the `Dx_y` subfolders and runs the make.down files within. Upon completion, all of the pts.same and pts.down files are in hand. It only remains to solve the full stack.

### Solving the Stack

_fin_

