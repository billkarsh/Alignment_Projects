
=============
General Notes
=============

// Note that printf with %.21f writes doubles to full precision
// such that scanf with %lf recovers the identical value.

- Submit brief jobs immed using " qsub -l short=true".

- Interactive: qlogin -l excl=true

- Can only occupy 2500 slots at once. Can only have 20,000 jobs
in the queue at once. Cost nominally $0.05/slot/hour.

- See ip addresses for shares: > mount.
- See host name: > hostname.

- ABComposite shows (BLUE pix) there are missed pixels in mappings
from A to B. This is an off by one error in some decision about
what's the interior.

- Fold masks are not yet generated appropriately for our images.

- Script scr should read the overlap requirements from Thumbparams.

- Try Gaussian blur instead of DoG, called from PixPair. This should
enhance larger features.

- MOS needs better name handling (that is, ween off of string labels
and onto idb z.id:rgn labels).

- LSQ should update connectivity (used flags) as inliers update.
- LSQ inlier removal may require update of connectivity.
- LSQ weight CPOINT by Q,R?
- LSQ take list of user CPOINT with stiffness?
- LSQ user points with long-range point matching!
- LSQ fit homography components on good data and force on user data.

- Tell what is actual overlap from running script on xml result,
but disable writing files. Not sure what the full thought was here.

- Handle folds and splits for optical.

- Foldmasks: Tiny:
- Results highly sensitive to sdev, which has a switch TINYSTAT to
specially calculate sum/n instead of sum/(n-1).
- Don't yet understand use of D parameter in tiny::ImageToFoldMap,
so my versions of the conn regions routines are currently commented out.

- There had been a deadlock problem in ptest waiting for mutex for
Thumbs table, but this seems to have abated with new cluster.

- Davi suggests that average values of homography tilt params
{g,h} might be forced for images that didn't get placed by lsq.

=================
Conceptual Issues
=================

Use #36, compare to sk sampl3
Use #30, compare to Munley - no grip, 36 much better

kart priority: 31, 36, 32, 38, ?{39, 34, 30, 33, 35}
Best yet: 21.183 in #36 (9/9/13)
21.292 in #32 (10/28/13)
21.060 in #31 (11/04/13)
Ribtect3 42R3
Valhalla Armadillo M
Sparco SPK-5 L
Alpinestars Bionic Rib S/L

- Agenda_101612.pptx good vision of big plan.

------------------------------------------------
Saalfeld Notes:
---------------
- Resume chat Wed Dec 11.

- Hackathon alert: sometime in April.

- I suggest we can manually extend cuts through fold masks but
recover correct connectivity in the cut zone by adding corr.
point pairs into the mix wherein the (x,y) coords are the same
for A::B, as are the {z,id} labels, but the rgn numbers differ.
For example: CPOINT2 z.id:1 x0 y0 z.id:2 x0 y0. Essentially, this
describes a single point labelled two ways.

- In dmesh, when finding corr points across layers, use even more
scale reduction (like 10 instead of ~2 as now). This is faster,
and focuses attention on larger features which are more stable.
It's a more appropriate size scale. The scale factor should be a
parameter that is decided by user mainly by slice thickness.

- In doing block-block matching, SS advocate that all good layers
be included because this effectively averages over more data so
reduces noise. It also adds longer range rigidity.

- SS suggest main places to extend current pipeline are:

- Incorporate use of software lenses per camera. These only need be
employed in CPixPair as images are loaded, and in the final rendering
phase, which is outside my scope.

- The model that is fit to the corr points should become a smooth
interpolation function like a spline rather than an affine. Details
of this remain vague to me right now.

- Yes, we need to break up huge connected regions into small patches
for computation, but the units we work with need not be natural
images. Rather, we should tile conn regions flexibly so that each
pixel is surrounded by the most possible context. So for example, if
a natural image contains only a small spike of protruding region that
offers too little area to match...then slide the "image" boundary
over to include more conn region and use that in dmesh. This ability
to paint image pairs (to be input to dmesh) from a large source conn
region also supports use of very large images (such as may come from
FIB SEM). Is this really needed? Lower priority than other features.
------------------------------------------------

- CLEAN OUT OLD DATA {Bocklab, Tomo}.

- -DTINYSTAT in tiny foldmask project sets a flag for maths.cpp
but that is compiled in separate archive. Need to revisit when
do foldmasks.

- Consider scatter/gather (vectored) IO for solution sync or
point loading.

- Unknown how large a block SuperLU can do.
- Unknown if SuperLU can work with more processors.
- Unknown how SuperLU breaks, but it looks like MPI breaks.
- Unknown how to glue blocks together in principle.
- Probly. need to lsq the scaffold first, then use that for blocks.
- Need way to draw montage foldmasks and assign to image-wise fm's.

- If fm's used, then the montage point-pairs can't be used in LSQ
because they're not yet labelled by region. So real workflow will
look more like:
	- make montages without foldmasks
	- calculate, draw, check, assign montage foldmasks
	- map montage foldmasks to tile foldmasks
	- montage again with full foldmasks to get point labels
	- now do all the cross layer work

- Need to visualize the block results to look for block errors.
- Bad block signs: ?{low R, low point cnt, some T consistency rule}?
- Need blocks to connect to more distant layers.
- Need blocks to connect to more proximal layers.
- Need study correlation drop-off with layer dZ.
- Need more threads all over.
- Need scape painter to observe fm's (and cropping).
- Need blocks to observe fm's.
- Need ptest to take params about which rgns go with this Tab hint.
- Need some documentation.
- Need methods paper, at some point.
- Need scipts to qsub in batches and stay within 20000 limit.

- Use microCT to make a scaffold. That is, align EM vs microCT to
make TAffineTable.txt external scaffold file to lsq.

- What is current status regarding tools for Davi to smooth intensity
distributions among neighbor tiles?

=============
Memory Issues
=============

- Presently, SuperLU working only by turning off the iterative
refinement option. This yield results identical to the original
single-processor SLU, which it turn out, never implemented
iterative refinement.

- Got result for 52 layer case using -nproc=32 and -nproc=96,
though in the latter case there were many logged instances of
"reg_mr Cannot allocate memory". Need to understand what this
means. Also note that these messages appeared only for the first
solve pass, not at all for the second, which means it is indep.
of the memory needs of the application; rather something else
is consuming memory.

- Need better debugging for SuperLU. Currently, extra folder
tomo/Libraries/SuperLU_DIST_3.3 is there to track changes I've made
for debugging in the deeper tomo/Libraries/SUPERLU/SuperLu_Dist_3.3.
I've also created working subfolder SRC_BAK to track debugging
changes that are different from the functionality changes tracked
in Bill_Changes folder. Need to be more selective about what gets
written and when. Probably should arrange to turn debugging lines
on dynamically instead of through static compile-time define.

-----------------------------------
Some web discussion:

The argument after -nn reflects the number of nodes the job is
running on, which on Vayu, should be set to ncpus/8.
You may find that jobs running across large numbers of CPUs fail
with an error similar to

v634:7702:6f8aa940: 136873455 us(136873455 us!!!):
reg_mr Cannot allocate memory

When COMSOL is running in distributed mode, it tends to send very
large messages, which can be larger than the amount of memory MPI
is able to allocate. In order to remedy this, the TCP/IP protocol
can be used for communications. This is done by adding the line

setenv I_MPI_FABRICS 'shm:tcp'

to your batch script. Please note that TCP/IP is slower than the
default communication protocol, and should only be used if errors
are encountered.
-----------------------------------


