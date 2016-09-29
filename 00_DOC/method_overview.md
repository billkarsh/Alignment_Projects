# Aligner Methods Overview

**Note:**

A practical example with additional discussion is available here: [Alignment_Tutorial](https://github.com/billkarsh/Alignment_Tutorial).

**Headings:**

* [Basic Concepts: Regions and Connecting Points](method_overview.md#basic-concepts-regions-and-connecting-points)
	* [Points File Content](method_overview.md#points-file-content)
	* [LSQ Solver](method_overview.md#lsq-solver)
* [Systematic Workflow to Make Point-Pairs](method_overview.md#systematic-workflow-to-make-point-pairs)
	* [Main Workflow Steps](method_overview.md#main-workflow-steps)
	* [IDB (A Structured Repository for Raw Data)](method_overview.md#idb-a-structured-repository-for-raw-data)
	* [Workspace 'temp' (Universal Filing Scheme for Machines and Humans)](method_overview.md#workspace-temp-universal-filing-scheme-for-machines-and-humans)
	* [Same-Layer Work (Montaging)](method_overview.md#same-layer-work-montaging)
	* [Cross-Layer Work (Strip Alignment)](method_overview.md#cross-layer-work-strip-alignment)
	* [Cross-Layer Work (Block Alignment)](method_overview.md#cross-layer-work-block-alignment)
* [Iterative LSQ Solver](method_overview.md#iterative-lsq-solver)
	* [LSQi Interface Part](method_overview.md#lsqi-interface-part)
	* [LSQw Worker](method_overview.md#lsqw-worker)
	* [Manual Convergence and Start/Stop Operation](method_overview.md#manual-convergence-and-Start-stop-operation)
	* [Packed Storage Formats](method_overview.md#packed-storage-formats)
	* [Solving Details and Parameters](method_overview.md#solving-details-and-parameters)
	* [Solution Viewer](method_overview.md#solution-viewer)
	* [Error Viewer](method_overview.md#error-viewer)


## <a name="basic-concepts-regions-and-connecting-points"></a>Basic Concepts: Regions and Connecting Points

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

### <a name="points-file-content"></a>Points File Content

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

### <a name="lsq-solver"></a>LSQ Solver

The input data needed to solve for the desired affines are these:

1. A specification of which regions to solve for; usually a whole number of layers, so specified as a range `z=i,j` of layer indices. These are the dependent variables.
2. Starting guesses for each of the transforms, which will be iteratively improved. We'll discuss where appropriate guesses come from in a later section.
3. The connecting point-pairs, which are treated as constants.

The solving algorithm is very simple. It looks like this:

#### For each 'iteration' of the solver...

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

## <a name="systematic-workflow-to-make-point-pairs"></a>Systematic Workflow to Make Point-Pairs

To feed the end-game process, just described, that derives the final global transforms from points that link regions (images) to each other we employ the following key devices:

1. A pipeline component (`ptest.exe` built by project `1_DMesh`) that takes three principal **inputs={images A, B, and approximate transform Tab from A to B}** and produces **output=tabulated points**.

2. A workflow (a sequence of bash scripts and binaries) designed to efficiently determine which images should be paired with each other (along with estimated Tab guesses) so that we can feed these pair-wise ptest jobs to a compute cluster in a parallelized way.

>The ptest program is such an important and complex workhorse that we devote a separate reference to it: [**ptest_reference**](ptest_reference.md). For the remainder of the workflow discussion let us just assume a magical black box taking image pairs to points tables.

### <a name="main-workflow-steps"></a>Main Workflow Steps

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

### <a name="idb-a-structured-repository-for-raw-data"></a>IDB (A Structured Repository for Raw Data)

An alignment session begins with a user provided file of information about the images. It can be a text file or TrakEM2 XML file. A text layout line has the following fields:

```
Z tileID a00 a01 a02 a10 a11 a12 col row cam full_path
```

* Z is 0-based layer index.
	* All the images for a layer should be grouped together in a layout file.
	* Listed Z's should never decrease.
	* Z's need not be contiguous: we will automatically match across gaps.
	* Z need not start at 0 (set desired `-z=i,j` range in topgo.sht).
* TileID is a non-negative signed 32-bit integer and is unique within its layer.
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
* Full rooted absolute path to image.
	* Images may be accessed from a variety of working directories so must be absolute.
	* Images can be TIF (8,16) PNG (8,16) MRC(16).
* Spaces are space or tab.

Script `dbgo.sht (project 1_MakeIDB)` reformats the user layout into a new IDB folder/file organization for efficient access by pipeline stages. The IDB device also has these features:

* Scripts/destinations for autogenerating 8-bit png images from mrc files.
* Scripts/destinations for autogenerating associated foldmasks.
* Designated location for foldmask metadata (counts of rgns in each image).
* _`Project 1_Tiny` builds pipeline component for mrc and foldmask operations._
* File `0_GEN/PipeFiles.cpp` provides API to access these data using {z,id} key-pair.

### <a name="workspace-temp-universal-filing-scheme-for-machines-and-humans"></a>Workspace 'temp' (Universal Filing Scheme for Machines and Humans)

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

### <a name="same-layer-work-montaging"></a>Same-Layer Work (Montaging)

Script `mongo.sht (project 1_MakeMontages)` not only creates a new 'temp0' workspace, it also uses the stage coordinates embedded in your initial layout file to determine which image pairs have sufficient overlap to list in a make.same script.

Internally, it actually creates a single long list of such pairs for each whole montage. Then, the pairs are partitioned into a grid of blocks according to `scriptparams::montageblocksize` and the top-left corner of each pair's A-member. An Sx_y subfolder is created for each grid cell along with a make.same file for the pair jobs in that block. It's just a divide and conquer scheme for parallelizing the work load. The assignment of jobs to this or that grid cell is not important. In fact, if any grid cells look underpopulated we reassign their members to neighbor cells to getter a better balance.

Typical make.same 'make rules' look like this:

```
4/1.8.map.tif:
	ptest >>pts.same 2>pair_1.4^1.8.log 1.4^1.8 -nf ${EXTRA}

3/1.9.map.tif:
	ptest >>pts.same 2>pair_1.3^1.9.log 1.3^1.9 -nf ${EXTRA}
```

Note that makefiles like this will indeed be given to the make command and in the past we may have actually used features of the make system such as checking the build date of a target like '4/1.8.map.tif'. However, these days, all we really want to do with the make facility is simply list the ptest jobs to be executed, and use make's -j parameter to say for example: `make -f make.same -j 8`, that is, on a given cluster node, run this list of jobs, running no more than eight jobs at a time. Again, it's a device for load balancing and resource management; nothing more complicated than that.

What you should note about the rules (jobs) in a make.same file is that there is a requisite pair of images (always with same layer index) and the only other parameter is -nf signifying that foldmasks are not being used. Conspicuously absent is any parameter of form [-Tab=a0,a1,a2,a3,a4,a5] to give a starting hint for A=>B. Rather, the absence of that hint indicates that ptest should look up the initial affines from the IDB and create Tab=[Tb-inv][Ta].

Montaging is highly reliable. If a pair of images in the same layer overlap and if that overlap contains any kind of features, then those features should actually match in both pictures. It doesn't matter if the features themselves are excellent quality tissue or are dirt on the sample; we still have two pictures of the same thing. Moreover, we can montage any layer in isolation without reference to other layers.

We make all the montages ahead of any cross-layer work so that we can render them (and parts thereof) to deduce how the layers fit together. Once we understand how the layers match up, we will be able to run make.down jobs and get high quality pts.down data. **_In the end, the final solve uses all the pts.same and pts.down data to solve for all the transforms denovo_**. This is often a point of confusion, but the montage transforms only assist with rough cross-layer positioning. **_Montage transforms are not used in the final solve_**. No claim is made that this is fundamentally necessary or optimal. It is simply an accurate description of how the present method works.

### <a name="cross-layer-work-strip-alignment"></a>Cross-Layer Work (Strip Alignment)

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

### <a name="cross-layer-work-block-alignment"></a>Cross-Layer Work (Block Alignment)

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

## <a name="iterative-lsq-solver"></a>Iterative LSQ Solver

The solver comprises two parts:

* `1_LSQi` user interface part
* `1_LSQw` worker part

### <a name="lsqi-interface-part"></a>LSQi Interface Part

LSQi primarily decides how to distribute the work onto one or many cluster nodes and then launches the job(s).

LSQi first crawls the temp folder hierarchy looking for S-folders, D-folders and the ThmPair files within D-folders. From this scan it builds a catalog file that tabulates for each layer, the maximum x and y grid cell indices for S- and D-types and the list of all other layers that this layer connects to (below it). Workers use the grid extent info to figure out which pts.same and pts.down files they need to load. LSQi, itself, uses the layer connectivity data to determine how the problem separates into blocks of whole layers.

Essentially, LSQi compares the number of layers you're solving for (zi=inner range) with parameter `-zpernode`. If zi fits within zpernode it's a single machine job and the next decision is whether to launch the work in-process on the current machine `(maxthreads=1` or `-local` option set) or to resubmit the job to the cluster to get the desired slot count.

If zi exceeds zpernode then `nwks` machines will be used. A custom Grid Engine environment is needed to reserve the required cluster resources. We've provided a kit (`00_impi3`) to help the system administrator set up the needed scripts. Here's the [ReadMe](../00_impi3/ReadMe.md) from the kit that explains how LSQi sets up and launches the MPI-based cluster job.

### <a name="lsqw-worker"></a>LSQw Worker

#### Unified Multi/Single Processing Workflow

LSQw always gets called with command line argument `-nwks=k` and this signals the code whether to start MPI and run in a multiprocessor fashion (nwks>1) or not. In all cases, the same source code does the work. This is managed in a simple way. Workers are numbered (global `wkid`) from 0 to n-1 and regardless of worker count, nwks, there is always wkid=0 who is implicitly designated the coordinator of all workers. The multi-worker and single worker flow are unified like this:

```
1. All workers (including 0) calculate your own share of the total
2. If( nwks > 1 )
		Wait for everyone else to catch up
3. If( I'm not zero )
		Send my results to zero
   Else if( nwks > 1 )
		Recv and consolidate results from others
4. Continue
```

All MPI operations are abstracted/wrapped by files `1_LSQw::lsq_MPI.{cpp,h}`. As it turns out, very few MPI functions are needed to efficiently handle all data sharing:

* MPIWaitForOthers(): A barrier (fence) that blocks until all workers reach same point
* MPISend(): Synchronously send a byte array to specified worker
* MPIRecv(): Synchronously recv a byte array from specified worker

#### Bag of Calculators

LSQw can be thought of as a collection of calculators. Several of them have the aforementioned structure/workflow, namely {calculate my part, consolidate at wkid=0}:

* Bounds: Calculate bounding rectangle for solution
* Dropout: Tabulate tiles lost for various reasons
* Error: Tabulate point mapping errors, RMS, largest deviants
* Magnitude: Calculate largest deviations of effective scale from unity
* Split: Resolve/separate solution into its disconnected islands (if any)

#### Untwist Calculator

Another calculator `Untwist` improves the accuracy of an input scaffold. Remember the scaffold is generally created from rough alignment of scaled-down montages and that is done before block alignment, hence, before down connection points are known. However, by the time that the final stack solve is performed, all of the points are in hand and we can use them to find improved rigids between the layers. Applying the rigids to the whole stack involves cumulative products. The total correction Rk applied to layer k is the product of foregoing pairwise corrections (rj):

`Rk = rk*rk-1*...*r1`.

In this case, pairwise transforms are communicated among workers by writing the rigids to wkid-labeled disk files. In a subsequent phase, each worker parses the files from other workers with wkid lower than its own.

#### Solve Calculator

The solving itself is really just another calculator, but it has a unique data sharing structure. Again, each worker is responsible for solving the transforms for its own "inner" range `zi=i,j` of layers. The zi are non-overlapping. Each layer is solved by just one worker. However, each worker must carry points and transforms for an "outer" range of layers `zo=p,q` that the inner layers depend upon. At startup each worker reads the `ranges.txt` file that LSQi writes to determine the index ranges of layers it must exchange with the worker to its left or right.

After each solve iteration a worker has updated its own zi layers, but it must get updated zo layers from its left/right neighbors. This is done in an orderly and simple manner by the `Updt()` method of the solution class: `XArray`. Here, 'odd' and 'even' refer to the worker wkids:

```
1. Odds send left,   evens recv from right
2. Odds send right,  evens recv from left
3. Evens send left,  odds recv from right
4. Evens send right, odds recv from left
```

### <a name="manual-convergence-and-Start-stop-operation"></a>Manual Convergence and Start/Stop Operation

The solver does not itself contain an automatic convergence mechanism. Rather, the vision for how the solver would be used is this:

```
1. Run for n1-thousand iters taking original scaffold 'X0' to solution X1
2. Spot check errors and smoothness. If not satisfied...
3. (Optionally modify parameters like -Wr or -Etol)...
4. Run n2-thousand iters more taking X1 to X2.
5. Back to (2): Evaluate and repeat as necessary.
```

The solver's `-mode=XXX` option specifies which operations to perform (which calculators to invoke) on this pass:

* `-mode=catalog`: Just calculate the block/layer connectivity file `lsqcache/catalog.txt`.
* `-mode=eval`: Calculate characteristics of solution X: {bounds, magnitudes, errors, drops}.
* `-mode=A2A`: Solve for new affines from given affines (Xprev = `-prior=XXX`).
* `-mode=A2H`: Solve for homographies from given affines.
* `-mode=H2H`: Solve for homographies from given homographies.
* `-mode=split`: Resolve X into its separate islands (typically the final operation).

>Note that any of the solve modes {A2A, A2H, H2H} implicitly also run `eval` upon completion.

### <a name="packed-storage-formats"></a>Packed Storage Formats

The start/stop/eval and internal iterative workflow suggested that care be taken to use efficient data representations to minimize both disk I/O and MPI data exchange. Compact data also facilitate tackling huge problems.

#### Binary Point Data

The correspondence points in pts.same and pts.down files need to be loaded for solve or evaluate operations and they are in human readable text formats in these files because that enables one to debug and possibly edit these critical data. However, they are never modified by the solver.

To improve I/O performance, after LSQw has launched and determined its layer range it looks for a cached binary form of its required points data. The file has name pattern: `lsqcache/pnts_wkid_zi_zj.bin`. If the file does not exist it is created from the equivalent text-based points files.

> Note: If you edit a text-based points file, remember to delete the cached binary file!

By default LSQw will look for the cached {catalog, points} data in a local subdirectory (of PWD) with name `lsqcache`. You can override that path with `-cache=altpath/lsqcache`. I do that frequently in the following usage scenario:

1. Solve **_normally_** in standard folder `temp/stack` with suggested `runlsq.sht`.
2. For comparison, create folder `temp/stack2` with copy of that runlsq.sht...
3. Edit `temp/stack2/runlsq.sht`
	* to make experimental parameter changes,
	* **and** to set `-lsqcache=../stack/lsqcache` to avoid remaking these data.

#### Solution Folders 'X'

The pipeline and solver generally represent collections of transforms as a folder with name pattern: `X_<A,H>_<BIN,TXT>[_optional-tag]`. For example, the cross-layer work begins by assembling (using `gathermons.sht`) all of the individual layer montages into a folder named `X_A_BIN_mons`.

* **X**: always present, signifies solution (AX=B)
* **A, H**: signifies affines (6 doubles) or homographies (8 doubles)
* **BIN, TXT**: signifies binary or text data
* Parsers that open these folders ignore anything following {BIN,TXT}

##### Binary Folders

Binary folders like `X_A_BIN` contain transform **and** flags files like this:

* `X_A_2814.bin`: All affines for layer 2814.
* `F_2814.bin`: All flag fields for those affines.

All transforms get an 8-bit flag field with these bits:
```
enum RgnFlags {
	fbDead		= 0x01,  // transform invalid
	fbRead		= 0x02,  // no starting guess exists
	fbPnts		= 0x04,  // too few pnts to solve
	fbKill		= 0x08,  // ill behaved...dropped
	fbCutd		= 0x10   // cross-layer links cut
};
```

Most data within the pipeline and solver are organized at the highest level by layer number. Within a layer, transforms, flags and other solver data are indexed by their (id,rgn) doublet using a standard C++ `map<int,int>` device, one map per layer. The left-hand tile-id is associated with a cumulative count of subregions over all tile-id values lower than that. All indices are zero-based throughout the code. The schematic index calculation is then:

`I(id,r) = map::find(id)->second + r - 1;`

The map is calculated by function `0_GEN/PipeFiles.cpp::IDBGetIDRgnMap( map, idb, z )` which scans tabulated counts of regions (in `fm.same` files) within the IDB.

##### Text Folders

Text versions of solution folders, like `X_A_TXT`, contain only transform files like:

* `X_A_2814.txt`: Affines for layer 2814.

There is a file per layer, of course. Within a file, the lines are labeled affines, so read: `id rgn A0 A1 ...`. In this case the flags are implicit. A listed transform exists with no attached issues, so gets flag value zero. All unlisted transforms are missing in action, so get flag value 0x03.

### <a name="solving-details-and-parameters"></a>Solving Details and Parameters

The schematic we gave earlier for solving is mostly true. There are a few minor modifications to describe in this section. To recap, the repetitive scheme is basically this:

1. Hold all tforms and points fixed except one tform (Ta).
2. Get an updated value for Ta by creating matrix equations with Ta(pa)=Tb(pb), for all point pairs having member pa belonging to Ta domain. This new Ta does not yet overwrite the previous Ta.
3. Repeat the single tform solve for each tforms in turn.
4. Now the collected new tforms replace the previous ones.
5. Updated tforms (and flags) are exchanged with neighbor workers.
6. Repeat this process for given number of iters.

One of the most difficult parts of the pipeline to get right was the block-block matching. In part that's because the extant data had many problems including highly variable intensity from tile to tile and large amounts of beam damage. This would on occasion make a bad guess (Tab) that would be propagated in the make.down file to the ptest jobs for all tiles in the block and that would make bad points which causes catastrophic results in the solver as you can imagine.

>Notes:

>1. Ptest does not have a reliable independent way to validate the Tab guess it gets as input. Again, Tab combined with the matchparams::LIX parameter create a correlation peak search mask and the correlator **_will_** return the peak value from that disk, correct or not. If the wrong thumbnail starting transform is calculated, the ensuing mesh work will probably be wrong yet points may be reported. There are opportunities in ptest to trap a bad result of course, like the correlation thresholds or the limits of mesh deformation and so on, but these tests are often too permissive. Sometimes you will make the tests permissive to get solutions for perfectly good cases that would be missed otherwise. As always, balancing parameters is an art. The methods are robust only if the data are of good quality.

>2. It can be maddening to track down which tile pairs are producing offending points because effects can be quite long range. Sometimes too, an errant tile may go unstable and be removed from the result only after it has already propagated mayhem to neighbors. Now, as you try to view the result the true offender isn't even in there.

In the next few sections we look at a few tricks to keep the inherently oscillatory solving process from getting caught in a loop, and to identify cancers at an early stage and remove them surgically.

#### Reducing Complexity

The solver is really not that complicated. The whole solving apparatus is contained in the single smallish file `1_LSQw/lsq_Solve.cpp`. Moreover, it's much simpler than that once you understand that much of the code is simply duplicated six times over to optimize addressing performance for the three fundamental solving modes. More specifically the three main solving loops `{_A2A, _A2H, _H2H}` are all copy/paste clones of each other with some macro names changed and some values of '6' exchanged for '8'. That gets us a 3X complexity reduction.

Next, each main solver has a matched companion, these are `{Cut_A2A, Cut_A2H, Cut_H2H}`. What happens in each main solver is that we first try to use all the points that connect this tile to others. If we incur a degenerate matrix or if the resulting transform looks excessively distorted our immediate suspicion is that the down connecting points may be bad due to a bad block-block match. Before giving up hope for this tform we call Cut_XXX which re-solves in the identical way, but using only in-layer connections, effectively "cutting" the tile's cross connections. That allows a tile to remain joined to its in-layer neighbors and move with them. It's not a perfect solution, but this way we can keep a tile in the set, marking it as "cutd" (cut downs) instead of writing it off as "killed". So you only need to understand _A2A to understand its five other copy/paste variants.

#### Solving Data

The global data we need are managed by `lsq_Globals` and consist of two main things:

1. The master list of all points `vector<CorrPnt> vC;` loaded from the binary points files.
2. The master list of all resident regions `vector<Rgns> vR;`

vC is a giant single zero-based list of these structures:

```
class CorrPnt {
public:
	Point   p1, p2;
	int     z1, z2,  // index into vR
            i1, i2;  // index into vR[iz]
	union {
	int		used;
	struct {
	uint16	r1, r2;
	};
	};
// methods not shown
};
```

Basically there are two linked points, each in their resp. local image coords, and each is associated with a region indexed by a pair of indices (z, i). Here, index i is not a tile-id but the result of mapping i = I(id, rgn) as described above.

Crucially, there is also a boolean `used` flag that marks this point as used Y/N. As solving progresses we may determine certain points to have too large an error to be included.

The regions (**_regions, not affines_**) are grouped first by zero-based layer index, so everything for a layer is addressed by `vR[iz]`. Within a layer, the zero-based single index is used to get to a region. The Rgns region structure carries these data:

```
class Rgns {
// The rgns for given layer
// indexed by 0-based 'idx0'
public:
	vector<vector<int> >  pts;      // ea rgn's pts
	vector<uint8>         flag;     // rgn flags
	map<int,int>          m;        // map id -> idx0
	int                   nr,       // num rgns
                          z;        // common z
// other members, methods not shown
};
```

Each region has a list `pts` of indices into the master list vC of point-pairs. All vC[i] having either member p1 or p2 belonging to this region get listed (initially). So each vC[i] is referenced by two regions. As solving progresses we may determine some points to be bad and we mark `vC[i].used = false`, moreover, we may remove that point from the region's pts list for efficiency.

The flags described earlier {fbRead, fbKill, ...} describe, hence, live with, the regions rather than the affines or homographies. Yes, they are exported with solutions to clarify whether a solution is valid or not.

Each layer gets a mapping i = I(id, rgn) as described earlier. Each entry in the map is the cumulative count of regions for all tile-id (id) below this one. But that leaves open what the total count of regions is (including the highest tile-id) so we need to also carry `nr` the total count. That's used for sizing arrays and for loop limits.

Finally, all indexing is converted to a zero-based scheme when LSQw initializes. The `z` member, however, is the real-world layer number which is needed for IDB lookups or reports.

The other important data for solving are the source and destination lists of affines (or homographies) which are XArrays like this:

```
class XArray {
public:
	int						NE;
	vector<vector<double> >	X;
public:
	void Resize( int ne );
	void Load( const char *path );
	void Save() const;
private:
	bool Send( int zlo, int zhi, int XorF, int toLorR );
	bool Recv( int zlo, int zhi, int XorF, int fmLorR );
public:
	bool Updt();
};
```

* `NE` is the number of elements (6 or 8).
* 'X' is addressed first by zero-based layer `X[iz]` then by combined region index (i). The start of the elements for tform (i) is `X[iz][NE*i]`.

#### Solving Algorithm

At the highest level the solving function is this:

```
void Solve( XArray &Xsrc, XArray &Xdst, int iters )
{
// set some mode-dependent parameters
// select mode proc, like _A2A
// determine how many threads `nthr` to use

// Iterate
	Xs = &Xsrc;
	Xd = &Xdst;
	for( pass = 0; pass < iters; ++pass ) {
		Do1Pass( proc );
		// swap Xs<->Xd
		XArray	*Xt = Xs; Xs = Xd; Xd = Xt;
	}
}
```

Here's the essence of Do1Pass. Note that as we solve we may decide to mark points as not used and we may decide to mark regions too, but the data are shared by many threads, all working on their own assigned regions (affines) so we can't really change the flag data while the threads are running. Rather, each thread keeps a private list of marks `vthr[ithr]` where it pushes edit requests. After all threads complete, the lists are consolidated and the edits are executed by `UpdateFlags()`. Finally data are exchanged with neighbor workers:

```
static void Do1Pass( EZThreadproc proc )
{
// multithreaded phase
	vthr.clear();
	vthr.resize( nthr );

	if( !EZThreads( proc, nthr, 1, "Solveproc" ) )
		exit( 42 );

// single-threaded phase
	UpdateFlags();
	vthr.clear();

// synchronize
	Xd->Updt();
}
```

Finally we look at a proc, for example `_A2A`. This is (like all my threads) a Posix pthread function. This a good time to mention that all threaded work in the code base is done using my `0_GEN/EZThreads` package. The general design pattern I use to assign N workers to M tasks is this:

* Reduce N until N <= M.
* Each thread/worker gets a zero-based id = j.
* Worker j works on tasks {j, j+N, j+N+N, ...} until j >= M.

_A2A uses a helper class `Todo Q` to do indexing. What's the first region I should work on? What's the next? Here the tasks are the regions and they straddle layers, so advancing to the next region is done modulo the previous vR[iz].nr.

Read the code itself along with the highlight commentary here:

```
static void* _A2A( void* ithr )
{
// Get the first region (Q.iz, Q.ir)
do { // for each region

// Make shorthand pointers and locals for my points list vp, point count np

// If too few points to solve (np < 3) mark this region bad and skip it.

// Create regularizer solution 'rgd', either translation or rigid. As we solve for
// affines we also solve for the more constrained rgd solution that, for example,
// can not shrink or skew. The output solution will be a component-wise admixture
// of the affine and rgd in the user-specified  proportion. For example, option
// -Wr=R,0.001 selects a rigid with weight 0.001.

// Get ready to build the normal equation matrix 'AX=B' where the sought transform
// components are X. A & B are termed the LHS and RHS matrices. Each point-pair
// makes a new equation Ta(pa) = Tb(pb) which is a new row in a virtual overdetermined
// matrix (more rows/points than variables). The AddConstraint_Quick() function builds
// LHS & RHS one equation at a time by calculating the product of A-transpose with
// the overdetermined A & B. When we 'Solve' RHS gets replaced by the solution X,
// so we simply point RHS into the correct portion of Xdst.

// For each point...

// Do some preliminary sanity checking

	// Skip if a !used point
	
	// Pick one of two functionally identical clauses according to whether p1 is
	// in the target A region or p2 is inside A. Another 2X complexity reduction!!
	
	// Employ index caching to do the fewest array accesses.
	// Make sure the 'other' point is in a valid region, or push an edit request.
	// Pair the appropriate tform with its point and apply to make global points.
	// Note that Ta and Tb are both the previous tforms, not yet the new solution
	// for Ta.
	
	// Check the error: the distance between these global points. If the error exceeds
	// command line option -Etol then do not use this point in the equations (but do
	// not mark it as problematic). IMPORTANT: As we iterate, the whole system will
	// tend to oscillate. Waves of change will wash back and forth through the volume.
	// We do not want to reject any points until the waves have died down a bit,
	// that is, until some number of iterations have completed. We call that ad hoc
	// parameter 'editdelay'.

// OK, we have a qualified point-pair...

// Next, we'll apply an oscillation damping measure. Let's call the global points
// from the previous iteration A = Ta(pa) and B = Tb(pb). Whereas we would normally
// solve for a new Ta'(pa) like this:
// Ta'(pa) = B, we instead do this:
// Ta'(pa) = B' = Wb*B + (1-Wb)*A, for Wb < 1.
// This is legal because prior to convergence the Ta and Tb are approximate, and
// as we get close to convergence, of course this modification tends toward zero.
// Its damping effect comes from effectively reducing the differences that the
// new Ta need to accommodate. Convergence is dramatically quicker.

// Now use the local pa and modified global B to add equations both for the
// regularizer and for the affine, and count the number of used points (++nu).

// Done with point loop

// ----------------------------
// Now solve and evaluate

// If the number of equations (points) was too small, KILL this region.

// Solve it...
// If the solution is pathological, try resolving with cut down points.
// Otherwise, solve for the regularizer and form the admixture.

// If enough iters have been done (somewhat damped down) then measure the
// 'squareness' of the transform, which means: Apply the tform to a right
// angle and report the deviation from 90 degrees, really, |sin(90-angle)|.
// If too deviant, again try resolving with cut down points.
// Otherwise,
// Check if the number of points used (nu) is fewer than the number listed (np)
// and if so then cull unused points from our list (we're the only thread
// currently using that list so we can do it now).

    } while (Q.Next());
}
```

### <a name="solution-viewer"></a>Solution Viewer

The pipeline and solver work primarily with binary versions of solution data, and the final output from the solve is binary, e.g., `X_A_BIN`. When you want to inspect these data use the provided viewer tool `1_XView`. You can convert IDB or X_folder data to various forms of text or make a TrakEM2 xml file from it.

### <a name="error-viewer"></a>Error Viewer

The solver creates binary tables of final point errors and reports these in output folder `Error`. Use the provided viewer tool `1_EView` to histogram these data. The output is a text file that can be imported into Excel or other plot program. The parameters let you specify the bin width and maximum error size. There is an overflow row. The columns in the output are:

* Error
* Counts for all points
* Counts for in-layer only
* Counts for cross-layer only

Moreover, you can histogram two Error folders at once, placing their respective result columns into one file for ready comparison in Excel. This is probably the most meaningful use of the error distribution: measuring if it has changed as a result of some intended improvement.

Remember that outlier errors may cast doubt on the validity of a solution transform, but may just as easily cast doubt on the validity of the points found by ptest.

_fin_

