Parallel_Programming_HW2 - Mid_Point
===========================
This document is used to show the detailed description of the Parallel_Programming_HW2 project.

See `https://github.com/EnzeXu/Parallel_Programming_HW2` or
```shell
$ git clone https://github.com/EnzeXu/Parallel_Programming_HW2.git
```

Compile and Run:
```shell
$ module load python/3.8.13
$ gcc mid_point.c -o a -fopenmp -lm
$ srun -p small --pty -N 1 -n 22 --ntasks-per-socket=22 --mem=80gb --time=0-10:00:00 /bin/bash
$ ./a 1000000000 $P_MAX$
```
