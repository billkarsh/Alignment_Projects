
=============
General Notes
=============

- Note that printf with %.21f writes doubles to full precision
such that scanf with %lf recovers the identical value.

- Submit brief jobs immed using "qsub -l short=true".
- Interactive: "qlogin -l excl=true".

- Can only occupy 2500 slots at once. Can only have 20,000 jobs
in the queue at once. Cost nominally $0.05/slot/hour.

- See ip addresses for shares: > mount.
- See host name: > hostname.

- Log into any node, e.g.: > ssh h04u23.
- Observe resource usage: > top.

- ABComposite shows (BLUE pix) there are missed pixels in mappings
from A to B. This is an off by one error in some decision about
what's the interior.

- MOS needs better name handling (that is, ween off of string labels
and onto idb z.id-rgn labels).

LSQ maybes:
- LSQ weight CPOINT by Q,R?
- LSQ take list of user CPOINT with stiffness?
- LSQ user points with long-range point matching!
- LSQ fit homography components on good data and force on user data.
- Davi suggests fit average homography {g,h}, force onto dropouts.

Foldmasks: Tiny:
- Handle folds and tears for optical.
- Fold masks are not yet generated appropriately for our images.
- Results highly sensitive to sdev, which has a switch TINYSTAT to
specially calculate sum/n instead of sum/(n-1).
- Don't yet understand use of D parameter in tiny::ImageToFoldMap,
so my versions of the conn regions routines are currently commented out.


=================
Conceptual Issues
=================

kart priority: 31, 36, 32, 38, ?{39, 34, 30, 33, 35}
Best yet:
21.183 in #36 (9/9/13)
21.292 in #32 (10/28/13)
21.060 in #31 (11/04/13)
21.188 in #31 (02/04/14)
Ribtect3 42R3
Valhalla Armadillo M
Sparco SPK-5 L
Alpinestars Bionic Rib S/L

- Agenda_101612.pptx good vision of big plan.

------------------------------------------------
Saalfeld Notes:
---------------
- I suggest we can manually extend cuts through fold masks but
recover correct connectivity in the cut zone by adding corr.
point pairs into the mix wherein the (x,y) coords are the same
for A::B, as are the {z,id} labels, but the rgn numbers differ.
For example: CPOINT2 z.id-1 x0 y0 z.id-2 x0 y0. Essentially, this
describes a single point labelled two ways.

- In dmesh, when finding corr points across layers, use even more
scale reduction (like 10 instead of ~2 as now). This is faster,
and focuses attention on larger features which are more stable.
It's a more appropriate size scale. The scale factor should be a
parameter that is decided by user mainly by slice thickness...
The counter argument is that low res matching is done at the block
level and what we want from image pairs is high res correspondence
points, so use highest res.

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
- Need more threads all over.
- Need scape painter to observe fm's (and cropping).
- Need blocks to observe fm's.
- Need ptest to take params about which rgns go with this Tab hint.
- Need some documentation.
- Need methods paper, at some point.
- Need scipts to qsub in batches and stay within 20000 limit.

- LSQw needs option -freset to clear flag history from prior.

- Use microCT to make a scaffold. That is, align EM vs microCT to
make external scaffold to lsq.

- What is current status regarding tools for Davi to smooth intensity
distributions among neighbor tiles?


