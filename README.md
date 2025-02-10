## Build

1. Set `PETSC_DIR` and `PETSC_ARCH` appropriately for your petsc installation (such that `$PETSC_DIR/$PETSC_ARCH` is the path to the installation directory).
2. `make`


## Running

The `settings.yaml` has most of the default options necessary.

Run using (for example):

```
mpirun -n 4 ./exec -options_file setting.yaml
```

Size of the problem is controlled by `-dm_plex_box_faces` flag. So a 100x100x100 mesh would be:

```
mpirun -n 4 ./exec -options_file setting.yaml -dm_plex_box_faces 100,100,100
```

Order of the elements is set via `-degree`. So quadratic elements would have `-degree 2`

The time taken to write the file is printed out at the end (via the `VecView` event in PETSc).

To have a non-zero solution written to the file, use flag `-project_solution`.
