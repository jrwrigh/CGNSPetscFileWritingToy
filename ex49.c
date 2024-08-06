static char help[] = "Tests dof numberings for external integrators such as LibCEED.\n\n";

#include <petscdmplex.h>
#include <petscds.h>

typedef struct {
  PetscBool useFE;
  PetscInt  degree;
  PetscBool closure_tensor;
} AppCtx;

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscFunctionBeginUser;
  options->useFE          = PETSC_TRUE;
  options->degree     = 1;
  options->closure_tensor = PETSC_FALSE;
  PetscOptionsBegin(comm, "", "Dof Ordering Options", "DMPLEX");
  PetscCall(PetscOptionsBool("-use_fe", "Use FE or FV discretization", "ex49.c", options->useFE, &options->useFE, NULL));
  PetscCall(PetscOptionsInt("-degree", "Degree of FEM discretization", "ex49.c", options->degree, &options->degree, NULL));
  PetscCall(PetscOptionsBool("-closure_tensor", "Use DMPlexSetClosurePermutationTensor()", "ex49.c", options->closure_tensor, &options->closure_tensor, NULL));
  PetscOptionsEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CreateMesh(MPI_Comm comm, AppCtx *user, DM *dm)
{
  PetscFunctionBeginUser;
  PetscCall(DMCreate(comm, dm));
  PetscCall(DMSetType(*dm, DMPLEX));
  PetscCall(DMSetFromOptions(*dm));
  PetscCall(DMSetApplicationContext(*dm, user));
  PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetupDiscretization(DM dm, AppCtx *user)
{
  DM       cdm = dm;
  PetscInt dim;
  MPI_Comm comm = PetscObjectComm((PetscObject)dm);

  PetscFunctionBeginUser;
  PetscCall(DMGetDimension(dm, &dim));

  PetscFE        fe, fe_face;

  PetscCall(PetscFECreateLagrange(comm, dim, 5, PETSC_FALSE, user->degree, user->degree, &fe));
  PetscCall(PetscFEGetHeightSubspace(fe, 1, &fe_face));
  PetscCall(DMAddField(dm, NULL, (PetscObject)fe));
  PetscCall(PetscFEDestroy(&fe));
  PetscCall(DMCreateDS(dm));
  PetscCall(DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL));
  while (cdm) {
    PetscCall(DMCopyDisc(dm, cdm));
    PetscCall(DMGetCoarseDM(cdm, &cdm));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}


int main(int argc, char **argv)
{
  DM       dm;
  AppCtx   user;
  PetscInt depth;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  PetscCall(ProcessOptions(PETSC_COMM_WORLD, &user));
  PetscCall(CreateMesh(PETSC_COMM_WORLD, &user, &dm));
  PetscCall(SetupDiscretization(dm, &user));
  {
    Vec vec;

    PetscCall(DMCreateLocalVector(dm, &vec));
    PetscCall(VecViewFromOptions(vec, NULL, "-vec_view"));
    PetscCall(VecDestroy(&vec));
  }
  PetscCall(DMPlexGetDepth(dm, &depth));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}

/*TEST

  test:
    suffix: 0
    requires: triangle
    args: -dm_refine 1 -petscspace_degree 1 -dm_view -offsets_view

  test:
    suffix: 1
    args: -dm_plex_simplex 0 -dm_plex_box_bd periodic,none -dm_plex_box_faces 3,3 -dm_sparse_localize 0 -petscspace_degree 1 \
          -dm_view -offsets_view

  test:
    suffix: cg_2d
    args: -dm_plex_simplex 0 -dm_plex_box_bd none,none -dm_plex_box_faces 3,3 -petscspace_degree 1 \
          -dm_view -offsets_view

  test:
    suffix: 1d_sfc
    args: -dm_plex_simplex 0 -dm_plex_dim 1 -dm_plex_shape zbox -dm_plex_box_faces 3 1 -dm_view -coord_ltog_view

  test:
    suffix: 2d_sfc
    nsize: 2
    args: -dm_plex_simplex 0 -dm_plex_dim 2 -dm_plex_shape zbox -dm_plex_box_faces 4,3 -dm_distribute 0 -petscspace_degree 1 -dm_view

  test:
    suffix: 2d_sfc_periodic
    nsize: 2
    args: -dm_plex_simplex 0 -dm_plex_dim 2 -dm_plex_shape zbox -dm_plex_box_faces 4,3 -dm_distribute 0 -petscspace_degree 1 -dm_plex_box_bd periodic,none -dm_view ::ascii_info_detail

  testset:
    args: -dm_plex_simplex 0 -dm_plex_dim 2 -dm_plex_shape zbox -dm_plex_box_faces 3,2 -petscspace_degree 1 -dm_view ::ascii_info_detail -closure_tensor
    nsize: 2
    test:
      suffix: 2d_sfc_periodic_stranded
      args: -dm_distribute 0 -dm_plex_box_bd none,periodic
    test:
      suffix: 2d_sfc_periodic_stranded_dist
      args: -dm_distribute 1 -petscpartitioner_type simple -dm_plex_box_bd none,periodic
    test:
      suffix: 2d_sfc_biperiodic_stranded
      args: -dm_distribute 0 -dm_plex_box_bd periodic,periodic
    test:
      suffix: 2d_sfc_biperiodic_stranded_dist
      args: -dm_distribute 1 -petscpartitioner_type simple -dm_plex_box_bd periodic,periodic

  test:
    suffix: fv_0
    requires: triangle
    args: -dm_refine 1 -use_fe 0 -dm_view -offsets_view

TEST*/
