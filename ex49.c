static char help[] = "Tests dof numberings for external integrators such as LibCEED.\n\n";

#include <petscdmplex.h>
#include <petscds.h>

typedef struct {
  PetscBool useFE;
  PetscInt  degree;
  PetscBool closure_tensor;
  PetscBool project_solution;
} AppCtx;

#define M_PI 3.1415926535897932384626433827950288

static PetscErrorCode project_function(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u,
                                       void *ctx) {
  PetscReal x_tot = 0;

  PetscFunctionBeginUser;
  for (PetscInt d = 0; d < dim; d++) x_tot += PetscSqr(x[d]) ;
  x_tot = sqrt(x_tot);
  for (PetscInt c = 0; c < Nc; c++) {
    PetscScalar value = sin(2.3 * M_PI * x_tot);
    if (PetscAbsScalar(value) < 100 * PETSC_MACHINE_EPSILON) value = 0.;
    u[c] = value;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

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
  PetscCall(PetscOptionsBool("-project_solution", "Project solution onto mesh", "ex49.c", options->project_solution, &options->project_solution, NULL));
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
    if (user.project_solution) {
      PetscErrorCode (*funcs)(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u,
                              void *ctx) = {project_function};
      PetscCall(DMProjectFunctionLocal(dm, 0, &funcs, NULL, INSERT_VALUES, vec));
    }
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
