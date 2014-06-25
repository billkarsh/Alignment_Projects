# "impi3" Grid Engine Parallel Environment

### Licensing

Certain materials in this folder fall under the "Sun Industry Standards Source License" which is now "retired, free and open source" as described in the wiki article:

<http://en.wikipedia.org/wiki/Sun_Industry_Standards_Source_License.desktop>.

### Implementation Notes

Rob Lines is the name of our cluster administrator at Janelia Farm. He put together this set of scripts to implement the grid engine parallel environment "-pe impi3 s" that our solver uses to start mpi jobs running on several whole cluster nodes at once. (The 's' parameter is given as slots = n_whole_nodes * slots_per_node).

The materials in this folder should enable a grid engine admin to similarly implement this operation mode.

Here's how our solver uses this:

If solver interface LSQi see's that option `-zpernode` is < layer range, it will spread the work onto multiple whole machines, the count of which is "nwks" as described in the following sections.

#### Script mpigo.sht

Lsqi first writes a bash script "mpigo.sht" with content like this:

```
#!/bin/sh

tail -n +2 sge.txt > hosts.txt

mpirun -perhost 1 -n 4 -machinefile hosts.txt lsqw -nwks=4 -temp=/nobackup/bock/karsh/fly400_v1/temp0 -cache=/nobackup/bock/karsh/fly400_v1/temp0/stack/lsqcache -prior=../cross_wkspc/X_A_TXT -mode=A2A -Wr=R,0 -Etol=500 -iters=10000 -splitmin=1000 -maxthreads=16 -untwist
```

The script's "tail" command opens file `sge.txt` and strips off the first line, saving the resulting list of host names as file `hosts.txt`.

The next command calls mpirun with its custom parameters:

```
"-perhost 1":               place just one communicator on each node.
"-n 4":                     number of nodes is 4 (= worker count nwks in this example).
"-machinefile hosts.txt":   list of reserved hosts.
```

Everything else is just the command line for our solver worker lsqw.

#### impi3 Node Reservation

After writing the mpigo.sht script detailed above, lsqi tells the scheduler to reserves the nodes like this:

`qsub -N lsqmpi -cwd -V -b y -o sge.txt -pe impi3 64 ./mpigo.sht`

Notes:

1. lsqi directs standard out to a file named `sge.txt` so mpigo.sht can find it.
2. lsqi invokes the parallel environment `-pe impi3 s`, with s = [nwks * scriptparams.txt::slotspernode]. In our case: nwks=4, slotspernode=16.
3. When mpigo.sht is called to execute, sge.txt will have already been written to the PWD. Here is the typical content:

```
-unique -catch_rsh /var/spool/uge814/h01u07/active_jobs/9447715.1/pe_hostfile
h01u07
h05u11
h06u13
h05u16
```


