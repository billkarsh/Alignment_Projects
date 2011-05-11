
- Submit brief jobs immed using " qsub -l short=true".

- Interactive: qlogin -l excl=true

- ABComposite shows (BLUE pix) there are missed pixels in mappings
from A to B. This is an off by one error in some decision about
what's the interior.

- If it should turn out that there are scale diffs between A and B
then the Interpolate and Distribute operations need scale args.

- Fold masks are not yet generated appropriately for our images.

- Correct order of operations should be:
	1) Make montages (thumb + mesh) without foldmasks
	2) Make foldmasks with global numbering
	3) Make cross plane (thumb and mesh)

- Script scr should read the overlap requirements from Thumbparams.

- Jobs issue from montage centers first, and work outward, but should
amount of overlap also be used to order jobs?

- Try Gaussian blur instead of DoG, called from PixPair. This should
enhance larger features.

- Whole pipeline should use common DB files mapping {z,id}<->name.
- MOS needs better name handling.

- LSQ should update connectivity (used flags) as inliers update.
- LSQ inlier removal may require update of connectivity.
- LSQ weight CPOINT by Q,R?
- LSQ take list of user CPOINT with stiffness?
- LSQ user points with long-range point matching!

- Test if we can make stacks by rigid thumbs alone--enough CPOINT?

- Tell what is actual overlap from running script on xml result,
but disable writing files.

- Lou requests that for large images, it may be a good idea to
check the final correlation from deformable mesh using the full
scale images instead of the 'scale' images. This is a pain because
I downsample in place and would have to reload originals.

- Handle folds and splits for optical.
- Explore Davi aligner.


- Foldmasks: Tiny:
- Results highly sensitive to std, which has a switch TINYSTAT to
specially calculate sum/n instead of sum/(n-1).
- Don't yet understand use of D parameter in tiny::ImageToFoldMap,
so my versions of the conn regions routines are currently commented out.

- Perhaps write only fm.png, not fm.tif.

- May be an issue that now, not all tf.txt are produced in symmetric
way so lou pipeline may not find one. How to fix that?

- Big angles in larva data.

- Git.


