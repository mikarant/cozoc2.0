

#include "context.h"
#include "constants.h"
#include "daslice.h"
#include "defs.h"
#include "io.h"
#include "ops.h"
#include <math.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <stdbool.h>
#define PI 3.141592654

Context new_context (Options const options, Files const files) {
    Context   ctx;
    const int ncid = files.ncid_in;

    ctx.ncid = ncid;

    /* Read dimensions */

    ctx.mx = files.dimsize[DIM_X];
    ctx.my = files.dimsize[DIM_Y];
    ctx.mz = files.dimsize[DIM_Z];
    ctx.mt = files.dimsize[DIM_T];

    /* Set up the distributed 3D array layout for 1- and 2-component
     * fields */

    DMDACreate3d (
        PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE,
        DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, ctx.mx, ctx.my, ctx.mz,
        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, 0, 0, 0, &ctx.da);

    DMDAGetReducedDMDA (ctx.da, 2, &ctx.da2);

    /* Set up the distributed 2D array layout (xy-plane) */

    create_subdm_plane (DMDA_Z, ctx.da, &ctx.daxy);

    KSPCreate (PETSC_COMM_WORLD, &ctx.ksp);
    KSPSetDM (ctx.ksp, ctx.da);
    KSPSetFromOptions (ctx.ksp);

    /* Allocate context arrays*/

    PetscMalloc1 (ctx.mz, &ctx.Pressure);
    PetscMalloc1 (ctx.my, &ctx.Coriolis_parameter);
    PetscMalloc1 (ctx.my, &ctx.Latitude);
    PetscMalloc1 (ctx.mt, &ctx.Time_coordinate);

    DMCreateGlobalVector (ctx.daxy, &ctx.Surface_pressure);

    DMCreateGlobalVector (ctx.da, &ctx.Temperature);
    VecDuplicate (ctx.Temperature, &ctx.Sigma_parameter);
    VecDuplicate (ctx.Temperature, &ctx.Surface_attennuation);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity);
    VecDuplicate (ctx.Temperature, &ctx.Diabatic_heating);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_advection);
    VecDuplicate (ctx.Temperature, &ctx.Temperature_advection);
    VecDuplicate (ctx.Temperature, &ctx.Diabatic_heating_tendency);
    VecDuplicate (ctx.Temperature, &ctx.Friction_u_tendency);
    VecDuplicate (ctx.Temperature, &ctx.Friction_v_tendency);
    VecDuplicate (ctx.Temperature, &ctx.Friction_u);
    VecDuplicate (ctx.Temperature, &ctx.Friction_v);
    VecDuplicate (ctx.Temperature, &ctx.Temperature_tendency);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency_v);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency_t);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency_f);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency_q);
    VecDuplicate (ctx.Temperature, &ctx.Vorticity_tendency_a);
    VecDuplicate (ctx.Temperature, &ctx.Total_omega);
    for (size_t i = 0; i < NUM_GENERALIZED_OMEGA_COMPONENTS; i++) {
        VecDuplicate (ctx.Temperature, &ctx.omega[i]);
    }

    DMCreateGlobalVector (ctx.da2, &ctx.Horizontal_wind);
    VecDuplicate (ctx.Horizontal_wind, &ctx.Friction);

    /* These are read here because they are constants throughout the
     * calculation */

    /* Time coordinate */
    {
        size_t start[1] = {0};
        size_t count[1] = {ctx.mt};
        file_read_array_double (
            ncid, "time", start, count, ctx.Time_coordinate);

        for (int i = 0; i < (int)ctx.mt; i++)
            ctx.Time_coordinate[i] *= (double)3600;
    }

    /* Pressure levels (z-coordinate) */
    {
        size_t start[1] = {0};
        size_t count[1] = {ctx.mz};
        file_read_array_double (ncid, "lev", start, count, 
ctx.Pressure);
    }

    /* Latitude) */
    {
      size_t start[1] = {0};
      size_t count[1] = {ctx.my};
      file_read_array_double (ncid, "lat", start, count, ctx.Latitude);
    }

    /* Coriollis parameter is taken to be a function of latitude, only */

    for (int j = 0; j < ctx.my; j++) {
      ctx.Latitude[j] *= PI/180.0;
      ctx.Coriolis_parameter[j] = 2*7.292e-5*sin(ctx.Latitude[j]);
    }

    /* Grid spacings */
    ctx.hx = 2.0*PI/ctx.mx;
    ctx.hy = PI/ctx.my;
    info("Grid spacings are hard-wired");

    ctx.hz = ctx.Pressure[1] - ctx.Pressure[0]; /* hz is negative!!! */

    ctx.first = max_of_size_t (0, options.first);
    ctx.last  = min_of_size_t (options.last, ctx.mt - 1);

    return ctx;
}

Vec new_vec (Context *ctx) {
    Vec vec;
    DMCreateGlobalVector (ctx->da, &vec);
    return vec;
}

int temperature (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec T,
    Vec Ttend, Context *ctx) {

    static Vec tmpvec   = 0;
    static Vec Tnext    = 0;
    static Vec Tprev    = 0;

    if (!Tnext) {    // The first step
        VecDuplicate (ctx->Temperature, &Tnext);
        VecDuplicate (ctx->Temperature, &Tprev);
        VecDuplicate (ctx->Temperature, &tmpvec);
    }

    if (step == first) {
        if (first == 0) {
            file_read_3d (ncid, step, "T", T);
            file_read_3d (ncid, step + 1, "T", Tnext);

            VecCopy (T, Ttend);
            VecAXPY (Ttend, -1.0, Tnext);
            VecScale (Ttend, -1.0 / (double)(t[step + 1] - t[step]));
        } else {
            file_read_3d (ncid, step - 1, "T", Tprev);
            file_read_3d (ncid, step, "T", T);
            file_read_3d (ncid, step + 1, "T", Tnext);

            VecCopy (Tprev, Ttend);
            VecAXPY (Ttend, -1.0, Tnext);
            VecScale (Ttend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    } else {
        if (step == mt - 1) {
            VecCopy (T, Tprev);
            VecCopy (Tnext, T);

            VecCopy (Tprev, Ttend);
            VecAXPY (Ttend, -1.0, T);
            VecScale (Ttend, -1.0 / (double)(t[step] - t[step - 1]));
        } else {
            VecCopy (T, Tprev);
            VecCopy (Tnext, T);
            file_read_3d (ncid, step + 1, "T", Tnext);

            VecCopy (Tprev, Ttend);
            VecAXPY (Ttend, -1.0, Tnext);
            VecScale (Ttend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Tprev);
        VecDestroy (&Tnext);
    }
    return (0);
}

int geopotential_height (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec Z,
    Vec Ztend, Context *ctx) {

    static Vec tmpvec   = 0;
    static Vec Znext    = 0;
    static Vec Zprev    = 0;

    if (!Znext) {    // The first step
        VecDuplicate (ctx->Geopotential_height, &Znext);
        VecDuplicate (ctx->Geopotential_height, &Zprev);
        VecDuplicate (ctx->Geopotential_height, &tmpvec);
    }

    if (step == first) {
        if (first == 0) {
            file_read_3d (ncid, step, "Z", Z);
            file_read_3d (ncid, step + 1, "Z", Znext);

            VecCopy (Z, Ztend);
            VecAXPY (Ztend, -1.0, Znext);
            VecScale (Ztend, -1.0 / (double)(t[step + 1] - t[step]));
        } else {
            file_read_3d (ncid, step - 1, "Z", Zprev);
            file_read_3d (ncid, step, "Z", Z);
            file_read_3d (ncid, step + 1, "Z", Znext);

            VecCopy (Zprev, Ztend);
            VecAXPY (Ztend, -1.0, Znext);
            VecScale (Ztend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    } else {
        if (step == mt - 1) {
            VecCopy (Z, Zprev);
            VecCopy (Znext, Z);

            VecCopy (Zprev, Ztend);
            VecAXPY (Ztend, -1.0, Z);
            VecScale (Ztend, -1.0 / (double)(t[step] - t[step - 1]));
        } else {
            VecCopy (Z, Zprev);
            VecCopy (Znext, Z);
            file_read_3d (ncid, step + 1, "Z", Znext);

            VecCopy (Zprev, Ztend);
            VecAXPY (Ztend, -1.0, Znext);
            VecScale (Ztend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    }

    // Convert geopotential tendency to geopotential height tendency
    VecScale (Ztend, 1 / 9.81);

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Zprev);
        VecDestroy (&Znext);
    }
    return (0);
}

int sigma_parameter (
    DM da, PetscInt mz, PetscScalar *p, Vec Tvec, Vec sigmavec) {
    const double   R   = Specific_gas_constant_of_dry_air;
    const double   c_p = Specific_heat_of_dry_air;
    PetscInt       zs, ys, xs, zm, ym, xm;
    PetscScalar ***T;
    PetscScalar ***sigma;

    DMDAVecGetArrayRead (da, Tvec, &T);
    DMDAVecGetArray (da, sigmavec, &sigma);

    DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                sigma[k][j][i] =
                    log (T[k][j][i]) - (R / c_p) * log (p[k] / 100000.0);
            }
        }
    }

    DMDAVecRestoreArray (da, sigmavec, &sigma);
    fpder (da, mz, NULL, p, sigmavec);
    DMDAVecGetArray (da, sigmavec, &sigma);

    for (int k = zs; k < zs + zm; k++) {
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                sigma[k][j][i] *= -R * T[k][j][i] / p[k];
            }
        }
    }

    DMDAVecRestoreArrayRead (da, Tvec, &T);
    DMDAVecRestoreArray (da, sigmavec, &sigma);

    return (0);
}

static int horizontal_wind_and_vorticity (
    const int ncid, const int step, DM da, DM da2, size_t my, PetscScalar hx,
    PetscScalar hy, PetscScalar *latitude, Vec tmpvec, Vec V, Vec zeta) {
    for (int i = 0; i < 2; i++) {
        char name[2][3] = {"U", "V"};
        file_read_3d (ncid, step, name[i], tmpvec);
        VecStrideScatter (tmpvec, i, V, INSERT_VALUES);
    }

    horizontal_rotor (da, da2, my, hx, hy, latitude, V, zeta);
    return (0);
}

int horizontal_wind_and_vorticity_and_vorticity_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, DM da, DM da2,
    size_t my, PetscScalar hx, PetscScalar hy, PetscScalar *latitude, Vec V, Vec zeta, Vec zetatend,
    Context *ctx) {

    static Vec tmpvec   = 0;

    static Vec Vnext    = 0;
    static Vec Vprev    = 0;
    static Vec zetanext = 0;
    static Vec zetaprev = 0;

    if (!Vnext) {    // The first step
        VecDuplicate (ctx->Horizontal_wind, &Vnext);
        VecDuplicate (ctx->Horizontal_wind, &Vprev);
        VecDuplicate (ctx->Vorticity, &zetanext);
        VecDuplicate (ctx->Vorticity, &zetaprev);
        VecDuplicate (ctx->Vorticity, &tmpvec);
    }

    if (step == first) {
        if (first == 0) {
            horizontal_wind_and_vorticity (
            ncid, step, da, da2, my, hx, hy, latitude, tmpvec, V, zeta);
            horizontal_wind_and_vorticity (
            ncid, step + 1, da, da2, my, hx, hy, latitude, tmpvec, Vnext, zetanext);

            VecCopy (zeta, zetatend);
            VecAXPY (zetatend, -1.0, zetanext);
            VecScale (zetatend, -1.0 / (double)(t[step + 1] - t[step]));
    } else {
            horizontal_wind_and_vorticity (
            ncid, step - 1, da, da2, my, hx, hy, latitude, tmpvec, Vprev, zetaprev);
            horizontal_wind_and_vorticity (
            ncid, step, da, da2, my, hx, hy, latitude, tmpvec, V, zeta);
            horizontal_wind_and_vorticity (
            ncid, step + 1, da, da2, my, hx, hy, latitude, tmpvec, Vnext, zetanext);

            VecCopy (zetaprev, zetatend);
            VecAXPY (zetatend, -1.0, zetanext);
            VecScale (zetatend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    } else {
        if (step == mt - 1) {
            VecCopy (V, Vprev);
            VecCopy (Vnext, V);
            VecCopy (zeta, zetaprev);
            VecCopy (zetanext, zeta);

            VecCopy (zetaprev, zetatend);
            VecAXPY (zetatend, -1.0, zeta);
            VecScale (zetatend, -1.0 / (double)(t[step] - t[step - 1]));
        } else {
            VecCopy (V, Vprev);
            VecCopy (Vnext, V);
            VecCopy (zeta, zetaprev);
            VecCopy (zetanext, zeta);
            horizontal_wind_and_vorticity (
            ncid, step + 1, da, da2, my, hx, hy, latitude, tmpvec, Vnext, zetanext);

            VecCopy (zetaprev, zetatend);
            VecAXPY (zetatend, -1.0, zetanext);
            VecScale (zetatend, -1.0 / (double)(t[step + 1] - t[step - 1]));
        }
    }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Vprev);
        VecDestroy (&Vnext);
        VecDestroy (&zetaprev);
        VecDestroy (&zetanext);
    }

    return (0);
}


int diabatic_heating (Context *ctx, const int ncid, const int step) {

    Vec           Q          = ctx->Diabatic_heating;
    diabatic_heating_tendency (ncid, step, ctx->first, ctx->mt, ctx->Time_coordinate,
                                 Q, ctx->Diabatic_heating_tendency,ctx);
    return (0);
}

int diabatic_heating_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec Q,
    Vec Qtend, Context *ctx) {

    static Vec tmpvec   = 0;
    static Vec Qnext    = 0;
    static Vec Qprev    = 0;

    // This function has been tailored for the non-accumulated tendency fields

    if (!Qnext) {    // The first step
        VecDuplicate (ctx->Diabatic_heating_tendency, &Qnext);
        VecDuplicate (ctx->Diabatic_heating_tendency, &Qprev);
        VecDuplicate (ctx->Diabatic_heating_tendency, &tmpvec);
    }

    // Take the average of the current and the next step
    if (step!=mt-1){
            file_read_3d (ncid, step, "Q", Q);
            file_read_3d (ncid, step + 1, "Q", Qnext);

            VecCopy (Q, Qtend);
            VecAXPY (Qtend, 0.5, Qnext);
            VecScale (Qtend, 1.0 / (double)(t[step + 1] - t[step]));
    // For the last time step, take only the last one
    } else{
      file_read_3d (ncid, step, "Q", Q);
      VecCopy (Q, Qtend);
      VecScale (Qtend, 1.0 /  t[step]);
    }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Qprev);
        VecDestroy (&Qnext);
    }

    return (0);
}

int friction (Context *ctx, const int ncid, const int step) {

    Vec           F          = ctx->Friction;
    Vec           Fu          = ctx->Friction_u;
    Vec           Fv          = ctx->Friction_v;

    friction_u_tendency (ncid, step, ctx->first, ctx->mt, ctx->Time_coordinate, "FU", Fu,
                       ctx->Friction_u_tendency ,ctx);

    VecStrideScatter (ctx->Friction_u_tendency, 0, F, INSERT_VALUES);

    friction_v_tendency (ncid, step, ctx->first, ctx->mt, ctx->Time_coordinate, "FV", Fv,
                        ctx->Friction_v_tendency, ctx);
    VecStrideScatter (ctx->Friction_v_tendency, 1, F, INSERT_VALUES);

    return (0);
}

int friction_u_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, const char *varname, Vec F,
    Vec Ftend, Context *ctx) {

    static Vec tmpvec   = 0;
    static Vec Fnext    = 0;
    static Vec Fprev    = 0;

    if (!Fnext) {    // The first step
        VecDuplicate (ctx->Friction_u_tendency, &Fnext);
        VecDuplicate (ctx->Friction_u_tendency, &Fprev);
        VecDuplicate (ctx->Friction_u_tendency, &tmpvec);
    }

    if (step!=mt-1){
            file_read_3d (ncid, step, varname, F);
            file_read_3d (ncid, step + 1, varname, Fnext);

            VecCopy (F, Ftend);
            VecAXPY (Ftend, 0.5, Fnext);
            VecScale (Ftend, 1.0 / (double)(t[step + 1] - t[step]));
    } else{
      file_read_3d (ncid, step, varname, F);
      VecCopy (F, Ftend);
      VecScale (Ftend, 1.0 /  t[step]);
    }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Fprev);
        VecDestroy (&Fnext);
    }
    return (0);
}

int friction_v_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, const char *varname, Vec F,
    Vec Ftend, Context *ctx) {

    static Vec tmpvec   = 0;
    static Vec Fnext    = 0;
    static Vec Fprev    = 0;

    if (!Fnext) {    // The first step
        VecDuplicate (ctx->Friction_v_tendency, &Fnext);
        VecDuplicate (ctx->Friction_v_tendency, &Fprev);
        VecDuplicate (ctx->Friction_v_tendency, &tmpvec);
    }

    if (step!=mt-1){
            file_read_3d (ncid, step, varname, F);
            file_read_3d (ncid, step + 1, varname, Fnext);

            VecCopy (F, Ftend);
            VecAXPY (Ftend, 0.5, Fnext);
            VecScale (Ftend, 1.0 / (double)(t[step + 1] - t[step]));
    } else{
      file_read_3d (ncid, step, varname, F);
      VecCopy (F, Ftend);
      VecScale (Ftend, 1.0 /  t[step]);
    }

    if (step == ctx->last) {
        VecDestroy (&tmpvec);
        VecDestroy (&Fprev);
        VecDestroy (&Fnext);
    }
    return (0);
}

/*int zwack_okossi (Context *ctx, Vec bvec) {

  DM           da        = ctx->da;
  DM           daxy      = ctx->daxy;
  const double R         = Specific_gas_constant_of_dry_air;
  PetscScalar* f         = ctx->Coriolis_parameter;
  Vec          zetagtend = ctx->Geostrophic_vorticity_tendency;
  Vec          tadv      = ctx->Temperature_advection;
  Vec          omega     = ctx->Total_omega;
  Vec          Q         = ctx->Diabatic_heating_tendency;
  Vec          sigma     = ctx->Sigma_parameter;
  Vec          temptend,tmpvec, tmpvec2, inttemp, meantemp, meanzvec, inttot;
  PetscScalar ***a, ***b, **c,***z, **meanz, ***barray, ***itarray;
  PetscInt     i, j, k, zs, ys, xs, zm, ym, xm;

  // Initialize some temporary vectors
  DMGetGlobalVector (da, &tmpvec);
  DMGetGlobalVector (da, &tmpvec2);
  DMGetGlobalVector (da, &inttemp);
  DMGetGlobalVector (da, &temptend);
  DMGetGlobalVector (da, &inttot);
  DMGetGlobalVector (daxy, &meantemp);
  DMGetGlobalVector (daxy, &meanzvec);

  VecZeroEntries(inttot);
  VecZeroEntries(meanzvec);
  VecZeroEntries(temptend);
  VecZeroEntries(meantemp);

  VecAXPY(temptend, 1.0, tadv);

  VecCopy(sigma, tmpvec);
  calculate_static_stability(ctx, tmpvec);
  VecPointwiseMult (tmpvec2, omega, tmpvec);

  VecAXPY (temptend, 1.0, tmpvec2);

  VecAXPY (temptend, 1.0, Q);

  plaplace (temptend, ctx);

  DMDAVecGetArray (da, temptend, &a);
  VecZeroEntries(inttemp);

  DMDAVecGetArray (da, inttemp, &b);
  DMDAVecGetArray (daxy, meantemp, &c);
  DMDAVecGetArrayRead (da, zetagtend, &z);
  DMDAVecGetArray (daxy, meanzvec, &meanz);
  DMDAVecGetArray (da, bvec, &barray);
  DMDAVecGetArray (da, inttot, &itarray);

  DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

  // Vertical average of geostrophic vorticity tendency (zetagtend)
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        meanz[j][i] = meanz[j][i] + z[k][j][i]; } } }

  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      meanz[j][i] = meanz[j][i] / (zs + zm - 1); } }

  // Multiply temptend by R/f and divide by k
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        a[k][j][i] =  (R/f[j]) * a[k][j][i] / (zs + zm - k);
      } } }

  // integrate temperature tendency
  for (k = zs + 1; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        b[k][j][i] = b[k-1][j][i] + (a[k][j][i] + a[k-1][j][i]) / 2.0;
      } } }

  // take vertical average of integrated temperature tendency
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        c[j][i] = c[j][i] + b[k][j][i];
      } } }

  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      c[j][i] = c[j][i] / (zs + zm - 1);
      } }

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        itarray[k][j][i] = -c[j][i] + b[k][j][i]; } } }

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        barray[k][j][i] = (meanz[j][i] + itarray[k][j][i]); } } }


  DMDAVecRestoreArray (da, temptend, &a);
  DMDAVecRestoreArray (da, inttemp, &b);
  DMDAVecRestoreArray (daxy, meantemp, &c);
  DMDAVecRestoreArray (da, bvec, &barray);
  DMDAVecRestoreArray (da, inttot, &itarray);
  DMDAVecRestoreArrayRead (da, zetagtend, &z);
  DMDAVecRestoreArray (daxy, meanzvec, &meanz);

  VecCopy(inttot, bvec);

  DMRestoreGlobalVector (da, &tmpvec);
  DMRestoreGlobalVector (da, &tmpvec2);
  DMRestoreGlobalVector (da, &inttemp);
  DMRestoreGlobalVector (da, &temptend);
  DMRestoreGlobalVector (da, &inttot);
  DMRestoreGlobalVector (daxy, &meantemp);
  DMRestoreGlobalVector (daxy, &meanzvec);

  return(0);
  }*/

int calculate_static_stability (Context *ctx, Vec bvec) {

  DM           da    = ctx->da;
  PetscScalar* p     = ctx->Pressure;
  const double R     = Specific_gas_constant_of_dry_air;
  PetscInt     i, j, k, zs, ys, xs, zm, ym, xm;
  Vec          avec;
  PetscScalar ***a, ***b;

  DMGetLocalVector (da, &avec);
  DMGlobalToLocalBegin (da, bvec, INSERT_VALUES, avec);
  DMGlobalToLocalEnd (da, bvec, INSERT_VALUES, avec);

  DMDAVecGetArrayRead (da, avec, &a);
  DMDAVecGetArray (da, bvec, &b);

  DMDAGetCorners (da, &xs, &ys, &zs, &xm, &ym, &zm);

  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        b[k][j][i] = (p[k] / R) * a[k][j][i];
      } } }

  DMDAVecRestoreArrayRead (da, avec, &a);
  DMDAVecRestoreArray (da, bvec, &b);

  DMRestoreLocalVector (da, &avec);

  return(0);
}

int vorticity_tendency_v (Context *ctx) {

  DM           da        = ctx->da;
  size_t       mz        = ctx->mz;
  PetscScalar* p         = ctx->Pressure;
  PetscScalar* f         = ctx->Coriolis_parameter;
  Vec          vorttend  = ctx->Vorticity_tendency_v;
  Vec          vadv      = ctx->Vorticity_advection;
  Vec          omega     = ctx->omega[GENERALIZED_OMEGA_COMPONENT_V];
  Vec          vort      = ctx->Vorticity;
  Vec          tmpvec, tmpvec2, tmpvec3,bvec;

  // Initialize some temporary vectors
  DMGetGlobalVector (da, &tmpvec);
  DMGetGlobalVector (da, &tmpvec2);
  DMGetGlobalVector (da, &tmpvec3);
  DMGetGlobalVector (da, &bvec);

  // Initialiaze  vorticity tendency with zeroes
  VecZeroEntries(vorttend);

  // add vorticity advection term
  VecAXPY(vorttend,1.0,vadv);

  // Calculate the vertical advection term (nr 2)
  VecCopy (vort,tmpvec);
  fpder (da, mz, NULL, p, tmpvec);
  VecPointwiseMult (tmpvec2, omega, tmpvec);

  // add vertical advection term (with minus sign as in the equation)
  VecAXPY (vorttend, -1.0, tmpvec2);

  // Calculate the divergence term
  VecCopy (vort, tmpvec);
  field_array1d_add (tmpvec, f, DMDA_Y);
  VecCopy (omega,tmpvec2);
  fpder (da, mz, NULL, p, tmpvec2);
  VecPointwiseMult (tmpvec3, tmpvec, tmpvec2);

  // Add the divergence term
  VecAXPY (vorttend, 1.0, tmpvec3);

  // Calculate tilting term
  tilting (tmpvec3, ctx->Horizontal_wind, omega, ctx);

  // Add tilting term
  VecAXPY (vorttend, 1.0, tmpvec3);

  DMRestoreGlobalVector (da, &tmpvec);
  DMRestoreGlobalVector (da, &tmpvec2);
  DMRestoreGlobalVector (da, &tmpvec3);
  DMRestoreGlobalVector (da, &bvec);

  return(0);
}

int vorticity_tendency_t (Context *ctx) {

  DM           da        = ctx->da;
  size_t       mz        = ctx->mz;
  PetscScalar* p         = ctx->Pressure;
  PetscScalar* f         = ctx->Coriolis_parameter;
  Vec          vorttend  = ctx->Vorticity_tendency_t;
  Vec          omega     = ctx->omega[GENERALIZED_OMEGA_COMPONENT_T];
  Vec          vort      = ctx->Vorticity;
  Vec          tmpvec, tmpvec2, tmpvec3,bvec;

  // Initialize some temporary vectors
  DMGetGlobalVector (da, &tmpvec);
  DMGetGlobalVector (da, &tmpvec2);
  DMGetGlobalVector (da, &tmpvec3);
  DMGetGlobalVector (da, &bvec);

  // Initialiaze  vorticity tendency with zeroes
  VecZeroEntries(vorttend);

  // Calculate the vertical advection term (nr 2)
  VecCopy (vort,tmpvec);
  fpder (da, mz, NULL, p, tmpvec);
  VecPointwiseMult (tmpvec2, omega, tmpvec);

  // add vertical advection term (with minus sign as in the equation)
  VecAXPY (vorttend, -1.0, tmpvec2);

  // Calculate the divergence term
  VecCopy (vort, tmpvec);
  field_array1d_add (tmpvec, f, DMDA_Y);
  VecCopy (omega,tmpvec2);
  fpder (da, mz, NULL, p, tmpvec2);
  VecPointwiseMult (tmpvec3, tmpvec, tmpvec2);

  // Add the divergence term
  VecAXPY (vorttend, 1.0, tmpvec3);

  // Calculate tilting term
  tilting (tmpvec3, ctx->Horizontal_wind, omega, ctx);

  // Add tilting term
  VecAXPY (vorttend, 1.0, tmpvec3);

  DMRestoreGlobalVector (da, &tmpvec);
  DMRestoreGlobalVector (da, &tmpvec2);
  DMRestoreGlobalVector (da, &tmpvec3);
  DMRestoreGlobalVector (da, &bvec);

  return(0);
}

int vorticity_tendency_f (Context *ctx) {

  DM           da        = ctx->da;
  size_t       mz        = ctx->mz;
  PetscScalar* p         = ctx->Pressure;
  PetscScalar* f         = ctx->Coriolis_parameter;
  Vec          vorttend  = ctx->Vorticity_tendency_f;
  Vec          omega     = ctx->omega[GENERALIZED_OMEGA_COMPONENT_F];
  Vec          F         = ctx->Friction;
  Vec          vort      = ctx->Vorticity;
  Vec          tmpvec, tmpvec2, tmpvec3,bvec;

  // Initialize some temporary vectors
  DMGetGlobalVector (da, &tmpvec);
  DMGetGlobalVector (da, &tmpvec2);
  DMGetGlobalVector (da, &tmpvec3);
  DMGetGlobalVector (da, &bvec);

  // Initialiaze  vorticity tendency with zeroes
  VecZeroEntries(vorttend);

  // Calculate the vertical advection term (nr 2)
  VecCopy (vort,tmpvec);
  fpder (da, mz, NULL, p, tmpvec);
  VecPointwiseMult (tmpvec2, omega, tmpvec);

  // add vertical advection term (with minus sign as in the equation)
  VecAXPY (vorttend, -1.0, tmpvec2);

  // Calculate the divergence term
  VecCopy (vort, tmpvec);
  field_array1d_add (tmpvec, f, DMDA_Y);
  VecCopy (omega,tmpvec2);
  fpder (da, mz, NULL, p, tmpvec2);
  VecPointwiseMult (tmpvec3, tmpvec, tmpvec2);

  // Add the divergence term
  VecAXPY (vorttend, 1.0, tmpvec3);

  // Calculate tilting term
  tilting (tmpvec3, ctx->Horizontal_wind, omega, ctx);

  // Add tilting term
  VecAXPY (vorttend, 1.0, tmpvec3);

  // Calculate friction term
  horizontal_rotor (da, ctx->da2, ctx->my, ctx->hx, ctx->hy, ctx->Latitude, F, tmpvec);

  // Add friction term
  VecAXPY (vorttend, 1.0, tmpvec);

  DMRestoreGlobalVector (da, &tmpvec);
  DMRestoreGlobalVector (da, &tmpvec2);
  DMRestoreGlobalVector (da, &tmpvec3);
  DMRestoreGlobalVector (da, &bvec);

  return(0);
}

int vorticity_tendency_q (Context *ctx) {

  DM           da        = ctx->da;
  size_t       mz        = ctx->mz;
  PetscScalar* p         = ctx->Pressure;
  PetscScalar* f         = ctx->Coriolis_parameter;
  Vec          vorttend  = ctx->Vorticity_tendency_q;
  Vec          omega     = ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q];
  Vec          vort      = ctx->Vorticity;
  Vec          tmpvec, tmpvec2, tmpvec3,bvec;

  // Initialize some temporary vectors
  DMGetGlobalVector (da, &tmpvec);
  DMGetGlobalVector (da, &tmpvec2);
  DMGetGlobalVector (da, &tmpvec3);
  DMGetGlobalVector (da, &bvec);

  // Initialiaze geostrophic vorticity tendency with zeroes
  VecZeroEntries(vorttend);

  // Calculate the vertical advection term
  VecCopy (vort,tmpvec);
  fpder (da, mz, NULL, p, tmpvec);
  VecPointwiseMult (tmpvec2, omega, tmpvec);

  // add vertical advection term (with minus sign as in the equation)
  VecAXPY (vorttend, -1.0, tmpvec2);

  // Calculate the divergence term
  VecCopy (vort, tmpvec);
  field_array1d_add (tmpvec, f, DMDA_Y);
  VecCopy (omega,tmpvec2);
  fpder (da, mz, NULL, p, tmpvec2);
  VecPointwiseMult (tmpvec3, tmpvec, tmpvec2);

  // Add the divergence term
  VecAXPY (vorttend, 1.0, tmpvec3);

  // Calculate tilting term
  tilting (tmpvec3, ctx->Horizontal_wind, omega, ctx);

  // Add tilting term
  VecAXPY (vorttend, 1.0, tmpvec3);

  DMRestoreGlobalVector (da, &tmpvec);
  DMRestoreGlobalVector (da, &tmpvec2);
  DMRestoreGlobalVector (da, &tmpvec3);
  DMRestoreGlobalVector (da, &bvec);

  return(0);
}

int vorticity_tendency_a (Context *ctx) {

  DM           da        = ctx->da;
  size_t       mz        = ctx->mz;
  PetscScalar* p         = ctx->Pressure;
  PetscScalar* f         = ctx->Coriolis_parameter;
  Vec          vorttend  = ctx->Vorticity_tendency_a;
  Vec          omega     = ctx->omega[GENERALIZED_OMEGA_COMPONENT_A];
  Vec          vort      = ctx->Vorticity;
  Vec          tmpvec, tmpvec2, tmpvec3,bvec;

  // Initialize some temporary vectors
  DMGetGlobalVector (da, &tmpvec);
  DMGetGlobalVector (da, &tmpvec2);
  DMGetGlobalVector (da, &tmpvec3);
  DMGetGlobalVector (da, &bvec);

  // Initialiaze geostrophic vorticity tendency with zeroes
  VecZeroEntries(vorttend);

  // Calculate the vertical advection term
  VecCopy (vort,tmpvec);
  fpder (da, mz, NULL, p, tmpvec);
  VecPointwiseMult (tmpvec2, omega, tmpvec);

  // add vertical advection term (with minus sign as in the equation)
  VecAXPY (vorttend, -1.0, tmpvec2);

  // Calculate the divergence term
  VecCopy (vort, tmpvec);
  field_array1d_add (tmpvec, f, DMDA_Y);
  VecCopy (omega,tmpvec2);
  fpder (da, mz, NULL, p, tmpvec2);
  VecPointwiseMult (tmpvec3, tmpvec, tmpvec2);

  // Add the divergence term
  VecAXPY (vorttend, 1.0, tmpvec3);

  // Calculate tilting term
  tilting (tmpvec3, ctx->Horizontal_wind, omega, ctx);

  // Add tilting term
  VecAXPY (vorttend, 1.0, tmpvec3);

  DMRestoreGlobalVector (da, &tmpvec);
  DMRestoreGlobalVector (da, &tmpvec2);
  DMRestoreGlobalVector (da, &tmpvec3);
  DMRestoreGlobalVector (da, &bvec);

  return(0);
}

void update_context (size_t step, Files ncfile, Context *ctx) {

    DM           da       = ctx->da;
    DM           da2      = ctx->da2;
    size_t       my       = ctx->my;
    size_t       mz       = ctx->mz;
    size_t       mt       = ctx->mt;
    PetscScalar  hx       = ctx->hx;
    PetscScalar  hy       = ctx->hy;
    double *     time     = ctx->Time_coordinate;
    PetscScalar *p        = ctx->Pressure;
    Vec          psfc     = ctx->Surface_pressure;
    Vec          Z        = ctx->Geopotential_height;
    Vec         T        = ctx->Temperature;
    Vec         Ttend    = ctx->Temperature_tendency;
    Vec          sigma    = ctx->Sigma_parameter;
    Vec          V        = ctx->Horizontal_wind;
    Vec          zeta     = ctx->Vorticity;
    Vec          zetatend = ctx->Vorticity_tendency;

    const int ncid = ctx->ncid;

    read2D (ncid, step, "SP", psfc);
    file_read_3d (ncid, step, "Z", Z);
    temperature (ncid, step, ctx->first, mt, time, T, Ttend, ctx);
    sigma_parameter (da, mz, p, T, sigma);
    horizontal_wind_and_vorticity_and_vorticity_tendency (
    ncid, step, ctx->first, mt, time, da, da2, my, hx, hy, ctx->Latitude,V, zeta,
        zetatend, ctx);
    diabatic_heating (ctx, ncid, step);
    friction (ctx, ncid, step);
}

void free_context (Context *ctx) {
    PetscFree (ctx->Pressure);
    PetscFree (ctx->Coriolis_parameter);

    VecDestroy (&ctx->Surface_pressure);

    VecDestroy (&ctx->Temperature);
    VecDestroy (&ctx->Sigma_parameter);
    VecDestroy (&ctx->Vorticity);
    VecDestroy (&ctx->Geopotential_height);
    VecDestroy (&ctx->Diabatic_heating);
    VecDestroy (&ctx->Temperature_tendency);
    VecDestroy (&ctx->Vorticity_tendency);
    for (size_t i = 0; i < NUM_GENERALIZED_OMEGA_COMPONENTS; i++) {
        VecDestroy (&ctx->omega[i]);
    }

    VecDestroy (&ctx->Horizontal_wind);
    VecDestroy (&ctx->Friction);

    KSPDestroy (&ctx->ksp);
    DMDestroy (&ctx->da);
}
