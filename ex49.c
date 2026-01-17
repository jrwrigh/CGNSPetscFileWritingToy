static char help[] = "Tests dof numberings for external integrators such as LibCEED.\n\n";

#include <petscdmplex.h>
#include <petscds.h>
#include <petscsf.h>

typedef struct {
  PetscBool useFE;
  PetscInt  degree;
  PetscInt  num_comps;
  PetscBool closure_tensor;
  PetscBool project_solution;
  PetscBool part_balance;
  PetscBool create_mat;
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
  options->degree         = 1;
  options->closure_tensor = PETSC_FALSE;
  options->num_comps      = 5;
  options->part_balance   = PETSC_TRUE;
  PetscOptionsBegin(comm, "", "Dof Ordering Options", "DMPLEX");
  PetscCall(PetscOptionsBool("-use_fe", "Use FE or FV discretization", "ex49.c", options->useFE, &options->useFE, NULL));
  PetscCall(PetscOptionsInt("-degree", "Degree of FEM discretization", "ex49.c", options->degree, &options->degree, NULL));
  PetscCall(PetscOptionsInt("-num_comps", "Number of components", "ex49.c", options->num_comps, &options->num_comps, NULL));
  PetscCall(PetscOptionsBool("-closure_tensor", "Use DMPlexSetClosurePermutationTensor()", "ex49.c", options->closure_tensor, &options->closure_tensor, NULL));
  PetscCall(PetscOptionsBool("-project_solution", "Project solution onto mesh", "ex49.c", options->project_solution, &options->project_solution, NULL));
  PetscCall(PetscOptionsBool("-set_partition_balance", "Set partition balance", "DMPlexSetPartitionBalance", options->part_balance, &options->part_balance, NULL));
  PetscCall(PetscOptionsBool("-create_mat", "Create Mat", "DMCreateMatrix", options->create_mat, &options->create_mat, NULL));
  PetscOptionsEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CreateMesh(MPI_Comm comm, AppCtx *user, DM *dm)
{
  PetscFunctionBeginUser;
  PetscCall(DMCreate(comm, dm));
  PetscCall(DMSetType(*dm, DMPLEX));
  PetscCall(DMPlexSetPartitionBalance(*dm, user->part_balance));
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

  PetscCall(PetscFECreateLagrange(comm, dim, user->num_comps, PETSC_FALSE, user->degree, user->degree, &fe));
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

// @brief Get information about a DM's local vector
PetscErrorCode DMGetLocalVectorInfo(DM dm, PetscInt *local_size, PetscInt *global_size, VecType *vec_type) {
  Vec V_loc;

  PetscFunctionBeginUser;
  PetscCall(DMGetLocalVector(dm, &V_loc));
  if (local_size) PetscCall(VecGetLocalSize(V_loc, local_size));
  if (global_size) PetscCall(VecGetSize(V_loc, global_size));
  if (vec_type) PetscCall(VecGetType(V_loc, vec_type));
  PetscCall(DMRestoreLocalVector(dm, &V_loc));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// @brief Get information about a DM's global vector
PetscErrorCode DMGetGlobalVectorInfo(DM dm, PetscInt *local_size, PetscInt *global_size, VecType *vec_type) {
  Vec V;

  PetscFunctionBeginUser;
  PetscCall(DMGetGlobalVector(dm, &V));
  if (local_size) PetscCall(VecGetLocalSize(V, local_size));
  if (global_size) PetscCall(VecGetSize(V, global_size));
  if (vec_type) PetscCall(VecGetType(V, vec_type));
  PetscCall(DMRestoreGlobalVector(dm, &V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DMPrintVecSizes(DM dm, AppCtx *user) {
  MPI_Comm comm = PetscObjectComm((PetscObject)dm);
  PetscMPIInt rank, comm_size;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCallMPI(MPI_Comm_size(comm, &comm_size));

  // Mesh
  const PetscInt num_comp_q = user->num_comps;
  PetscInt       owned_dofs, local_dofs, global_dofs;
  PetscCall(DMGetGlobalVectorInfo(dm, &owned_dofs, &global_dofs, NULL));
  PetscCall(DMGetLocalVectorInfo(dm, &local_dofs, NULL, NULL));
  PetscCall(PetscPrintf(comm, "\n"));
  PetscCall(PetscPrintf(comm, "Partition:                             (min,max,median,max/median)\n"));
  {
    PetscInt *gather_buffer = NULL;
    PetscInt  part_owned_dofs[3], part_local_dofs[3], part_boundary_dofs[3], part_neighbors[3];
    PetscInt  median_index = comm_size % 2 ? comm_size / 2 : comm_size / 2 - 1;
    if (!rank) PetscCall(PetscMalloc1(comm_size, &gather_buffer));

    PetscCall(PetscPrintf( comm, "  Total Global Dofs                  : %" PetscInt_FMT "\n", global_dofs));

    PetscCallMPI(MPI_Gather(&owned_dofs, 1, MPIU_INT, gather_buffer, 1, MPIU_INT, 0, comm));
    if (!rank) {
      PetscCall(PetscSortInt(comm_size, gather_buffer));
      part_owned_dofs[0]             = gather_buffer[0];              // min
      part_owned_dofs[1]             = gather_buffer[comm_size - 1];  // max
      part_owned_dofs[2]             = gather_buffer[median_index];   // median
      PetscReal part_owned_dof_ratio = (PetscReal)part_owned_dofs[1] / (PetscReal)part_owned_dofs[2];
      PetscCall(PetscPrintf(
          comm, "  Global Vector %" PetscInt_FMT "-DoF nodes          : %" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ", %f\n", num_comp_q,
          part_owned_dofs[0] / num_comp_q, part_owned_dofs[1] / num_comp_q, part_owned_dofs[2] / num_comp_q, part_owned_dof_ratio));
    }

    PetscCallMPI(MPI_Gather(&local_dofs, 1, MPIU_INT, gather_buffer, 1, MPIU_INT, 0, comm));
    if (!rank) {
      PetscCall(PetscSortInt(comm_size, gather_buffer));
      part_local_dofs[0]             = gather_buffer[0];              // min
      part_local_dofs[1]             = gather_buffer[comm_size - 1];  // max
      part_local_dofs[2]             = gather_buffer[median_index];   // median
      PetscReal part_local_dof_ratio = (PetscReal)part_local_dofs[1] / (PetscReal)part_local_dofs[2];
      PetscCall(PetscPrintf(
          comm, "  Local Vector %" PetscInt_FMT "-DoF nodes           : %" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ", %f\n", num_comp_q,
          part_local_dofs[0] / num_comp_q, part_local_dofs[1] / num_comp_q, part_local_dofs[2] / num_comp_q, part_local_dof_ratio));
    }

    if (comm_size != 1) {
      PetscInt num_remote_roots_total = 0, num_remote_leaves_total = 0, num_ghost_interface_ranks = 0, num_owned_interface_ranks = 0;
      {
        PetscSF            sf;
        PetscMPIInt        nrranks, niranks;
        const PetscInt    *roffset, *rmine, *rremote, *ioffset, *irootloc;
        const PetscMPIInt *rranks, *iranks;
        PetscCall(DMGetSectionSF(dm, &sf));
        PetscCall(PetscSFSetUp(sf));
        PetscCall(PetscSFGetRootRanks(sf, &nrranks, &rranks, &roffset, &rmine, &rremote));
        PetscCall(PetscSFGetLeafRanks(sf, &niranks, &iranks, &ioffset, &irootloc));
        for (PetscInt i = 0; i < nrranks; i++) {
          if (rranks[i] == rank) continue;  // Ignore same-part global->local transfers
          num_remote_roots_total += roffset[i + 1] - roffset[i];
          num_ghost_interface_ranks++;
        }
        for (PetscInt i = 0; i < niranks; i++) {
          if (iranks[i] == rank) continue;
          num_remote_leaves_total += ioffset[i + 1] - ioffset[i];
          num_owned_interface_ranks++;
        }
      }
      PetscCallMPI(MPI_Gather(&num_remote_roots_total, 1, MPIU_INT, gather_buffer, 1, MPIU_INT, 0, comm));
      if (!rank) {
        PetscCall(PetscSortInt(comm_size, gather_buffer));
        part_boundary_dofs[0]           = gather_buffer[0];              // min
        part_boundary_dofs[1]           = gather_buffer[comm_size - 1];  // max
        part_boundary_dofs[2]           = gather_buffer[median_index];   // median
        PetscReal part_shared_dof_ratio = (PetscReal)part_boundary_dofs[1] / (PetscReal)part_boundary_dofs[2];
        PetscCall(PetscPrintf(
            comm, "  Ghost Interface %" PetscInt_FMT "-DoF nodes        : %" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ", %f\n",
            num_comp_q, part_boundary_dofs[0] / num_comp_q, part_boundary_dofs[1] / num_comp_q, part_boundary_dofs[2] / num_comp_q,
            part_shared_dof_ratio));
      }

      PetscCallMPI(MPI_Gather(&num_ghost_interface_ranks, 1, MPIU_INT, gather_buffer, 1, MPIU_INT, 0, comm));
      if (!rank) {
        PetscCall(PetscSortInt(comm_size, gather_buffer));
        part_neighbors[0]              = gather_buffer[0];              // min
        part_neighbors[1]              = gather_buffer[comm_size - 1];  // max
        part_neighbors[2]              = gather_buffer[median_index];   // median
        PetscReal part_neighbors_ratio = (PetscReal)part_neighbors[1] / (PetscReal)part_neighbors[2];
        PetscCall(PetscPrintf(comm, "  Ghost Interface Ranks              : %" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ", %f\n",
                              part_neighbors[0], part_neighbors[1], part_neighbors[2], part_neighbors_ratio));
      }

      PetscCallMPI(MPI_Gather(&num_remote_leaves_total, 1, MPIU_INT, gather_buffer, 1, MPIU_INT, 0, comm));
      if (!rank) {
        PetscCall(PetscSortInt(comm_size, gather_buffer));
        part_boundary_dofs[0]           = gather_buffer[0];              // min
        part_boundary_dofs[1]           = gather_buffer[comm_size - 1];  // max
        part_boundary_dofs[2]           = gather_buffer[median_index];   // median
        PetscReal part_shared_dof_ratio = (PetscReal)part_boundary_dofs[1] / (PetscReal)part_boundary_dofs[2];
        PetscCall(PetscPrintf(
            comm, "  Owned Interface %" PetscInt_FMT "-DoF nodes        : %" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ", %f\n",
            num_comp_q, part_boundary_dofs[0] / num_comp_q, part_boundary_dofs[1] / num_comp_q, part_boundary_dofs[2] / num_comp_q,
            part_shared_dof_ratio));
      }

      PetscCallMPI(MPI_Gather(&num_owned_interface_ranks, 1, MPIU_INT, gather_buffer, 1, MPIU_INT, 0, comm));
      if (!rank) {
        PetscCall(PetscSortInt(comm_size, gather_buffer));
        part_neighbors[0]              = gather_buffer[0];              // min
        part_neighbors[1]              = gather_buffer[comm_size - 1];  // max
        part_neighbors[2]              = gather_buffer[median_index];   // median
        PetscReal part_neighbors_ratio = (PetscReal)part_neighbors[1] / (PetscReal)part_neighbors[2];
        PetscCall(PetscPrintf(comm, "  Owned Interface Ranks              : %" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ", %f\n",
                              part_neighbors[0], part_neighbors[1], part_neighbors[2], part_neighbors_ratio));
      }
    }

    if (!rank) PetscCall(PetscFree(gather_buffer));
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
  if (user.create_mat) {
    Mat A;
    PetscCall(DMCreateMatrix(dm, &A));
    PetscCall(MatDestroy(&A));
  }
  PetscCall(DMPrintVecSizes(dm, &user));

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
