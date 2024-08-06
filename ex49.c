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

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  PetscCall(PetscLogDefaultBegin());  // So we can use PetscLogEventGetPerfInfo without -log_view
  PetscCall(ProcessOptions(PETSC_COMM_WORLD, &user));
  PetscCall(CreateMesh(PETSC_COMM_WORLD, &user, &dm));
  PetscCall(SetupDiscretization(dm, &user));
  {
    Vec vec;

    PetscCall(DMCreateLocalVector(dm, &vec));
    PetscCall(VecViewFromOptions(vec, NULL, "-vec_view"));
    PetscCall(VecDestroy(&vec));
  }

  { // Get time taken to write file
    PetscLogEvent      vec_view_log;
    PetscEventPerfInfo perf_info;

    PetscCall(PetscLogEventGetId("VecView", &vec_view_log));
    PetscCall(PetscLogEventGetPerfInfo(PETSC_DETERMINE, vec_view_log, &perf_info));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\nTime taken for vec writing (sec): %g\n", perf_info.time));
  }

  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
