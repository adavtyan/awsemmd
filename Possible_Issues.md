# [Issue #1](https://code.google.com/p/awsemmd/issues/detail?id=#1) #

If you see an error like the one shown bellow and you did make the dummy MPI library in directory LAMMPSDIR/src/STUBS/ try to find a _`*`.a_ file in there. There is a chance that it has a different name (_libmpi\_stubs.a_ for example). Renaming it to _libmpi.a_ should fix the issue.

```
g++: ../STUBS/libmpi.a: No such file or directory
make[1]: *** [../lmp_serial] Error 1
make: *** [serial] Error 2
```