
=============
General Notes
=============

- Submit brief jobs immed using " qsub -l short=true".

- Interactive: qlogin -l excl=true

- ABComposite shows (BLUE pix) there are missed pixels in mappings
from A to B. This is an off by one error in some decision about
what's the interior.

- Fold masks are not yet generated appropriately for our images.

- Correct order of operations should be:
	1) Make montages (thumb + mesh) without foldmasks
	2) Make foldmasks with global numbering
	3) Make cross plane (thumb and mesh)

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

kart priority: 32, 31, 30, ?{38, 36, 34, 35}

- Agenda_101612.pptx good vision of big plan.

- CLEAN OUT OLD DATA {Bocklab, Tomo}.

- Unknown how large a block SUperLU can do.
- Unknown if SuperLU can work with more processors.
- Unknown how SuperLU breaks, but it looks like MPI breaks.
- Unknown how to glue blocks together in principle.
- Probly. need to lsq the scaffold first, then use that for blocks.
- Need way to draw montage foldmasks and assign to image-wise fm's.
- Need to visualize the block results to look for block errors.
- Bad block signs: ?{low R, low point cnt}?
- Need blocks to connect to more distant layers.
- Need study correlation drop-off with layer dZ.
- Need scape painter to observe fm's (and cropping).
- Need blocks to observe fm's.
- Need ptest to take params about which rgns go with this Tab hint.
- Need some documentation.
- Need methods paper, at some point.

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


