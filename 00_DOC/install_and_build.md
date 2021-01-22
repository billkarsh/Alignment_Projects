
# Install & Build Alignment_Projects

Repo URL: <https://github.com/billkarsh/Alignment_Projects>

1. Click 'Fork' on my github page
2. `> git clone <'clone url'>`

### Required 3rd Party Support Libraries

* libtiff
* libpng
* zlib
* One of: {Intel MKL, fftw}

### Set Your Local Path Data

- Decide on (mkdir) directory to hold all your build output binaries.
- Add that to your bash_profile/PATH environment variable.

- Copy file **00_ENV/aln_makefile_std_defs** to some directory of yours that is outside of the git folder, so updates to git will not overwrite your custom path data.

- Edit your copy of aln_makefile_std_defs...
    * Set PATH_SRCCODE to your local copy of 'Alignment_Projects'.
    * Set PATH_OUT to your choice from step (1).
    * If not at Janelia, then you need to point the PATH_XXX variables to your local copies of library and include directories for {libtiff, libpng, zlib}. You will also need to do that for the fftw package, unless you can use Intel MKL support for FFT operations. MKL is preferred because the code is fully re-entrant allowing more efficient multithreading.
    * Follow other instructions in aln_makefile_std_defs to link against libraries of your choice.


- All our make files have this include statement: `include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs`. So you must edit your bash_profile and define ALN_LOCAL_MAKE_PATH:
`export ALN_LOCAL_MAKE_PATH=/groups/~location_of_my_std_defs`. *(Example bash_profile edits below)*.

### Customize Cluster 'qsub' Dispatch Calls

- Copy files **00_ENV/QSUB_{1NODE,MNODE}.sht** to a directory defined in your PATH environment variable. These define how our software should submit 'qsub' commands to your local cluster OS. You should only need to select the appropriate cluster type in **QSUB_1NODE.sht**, and possibly edit the corresponding **QSUB_1NODE_XXX.sht** to make adjustments for your policies. You should not need to edit the projects that call these cluster dispatch files.

- MPI is not required unless you need to solve very large stacks on multiple machines. To enable that your system administrator must set up your Grid Engine with the requisite "impi3" parallel environment. Folder **00_impi3** contains resources for that. This is invoked from the QSUB_MNODE.sht script.

### Edit .bash_profile

Ours looks like this:

```
# ---------------------------------------------------------
# Alignment specific definitions
#

# Set path for local builds
export ALN_LOCAL_MAKE_PATH=/groups/apig/tomo

# Set universal mrc image file margin
export MRC_TRIM=12

# Append paths to needed binaries
PATH=$PATH:$HOME/bin
PATH=$PATH:/groups/apig/tomo/lou_stuff/binaries
PATH=$PATH:/groups/apig/tomo/commands
export PATH

# --- prior to 10/06/2016 ---
# Set place to look for libpng14 library
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64
#export LD_LIBRARY_PATH

# Disable Intel MKL internal threading
export MKL_NUM_THREADS=1

#
# End alignment definitions
# ---------------------------------------------------------
```

### Edit .bash_rc

Compiling against MPI libraries requires environment symbol definitions.

Here's what that looks like at Janelia:

```
# ---------------------------------------------------------
# .bashrc
#

# User specific aliases and functions

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi


# Use Intel compiler for MPI
if [ -f /usr/local/INTEL2016.sh ]; then
     . /usr/local/INTEL2016.sh
fi
export I_MPI_DEBUG=0
export I_MPI_FABRICS=shm:ofa
export I_MPI_FALLBACK=disable
# --- prior to 10/06/2016 ---
# export I_MPI_RENDEZVOUS_RDMA_WRITE=1
# export I_MPI_DEVICE=rdssm:OpenIB-iwarp
# export I_MPI_FALLBACK_DEVICE=0
# export I_MPI_USE_DYNAMIC_CONNECTIONS=0
# --- older ---
# export I_MPI_OFA_USE_XRC=1
# export I_MPI_OFA_DYNAMIC_QPS=1

#
# End bashrc
# ---------------------------------------------------------
```

Here's the contents of our **INTEL2016.sh** setup script:

```
. /usr/local/intel-2016/bin/compilervars.sh intel64
. /usr/local/intel-2016/impi/5.1.3.181/bin64/mpivars.sh intel64
. /usr/local/INTEL-2016/advisor_xe/advixe-vars.sh quiet
export INTEL_LICENSE_FILE="nnnnn@flexlm.int.janelia.org"
```

### Build

Project labeling scheme:

| Folder | Use |
| ---- | ---- |
| 00_XXX | Infrastructure |
| 0_XXX | Common libraries |
| 1_XXX | Main alignment pipeline |
| 2_XXX | Handy utilities |
| macros | Example scripts |
| Other | Bill Development stuff |

Good idea to build the code on the target machine. Log on to a cluster node and build there.

Useful targets for the make command:

```
> make all1     ;build all 1_XXX
> make all2     ;build all 2_XXX
> make          ;build 1_ & 2_ (same as 'make all')
> make clean    ;remove all {.o,.a} but keep .exe
```

Standard way to make everything:

```
> cd 00_ENV
> make
```

### Script type '.sht'

Scripts generally use bash shell (.sh), but I choose to name my files (.sht) for 'bash text'. In my linux GUI (KDE) I have mapped double-clicks on (.sht) files to KEDIT.exe because I more often read and edit scripts than run them (so I can check what I'm about to do).

_fin_

### Docker build instructions

```
docker build -t $(USER)\alignment_projects .
```

then you should have an ubuntu image with MKL installed and compiled binaries at \usr\local\bin

