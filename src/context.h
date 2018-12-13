#pragma once

#include "io.h"
#include "omega.h"
#include "options.h"
#include <petscksp.h>

typedef struct Context Context;
struct Context {
    int          ncid; // Input file
    size_t       current, first, last;
    DM           da, da2;
    DM           daxy;
    KSP          ksp;
    PetscInt     mx, my, mz, mt;    // Global grid sizes
    PetscScalar  hx, hy, hz;        // Grid spacings
    double *     Time_coordinate;    // In seconds
    PetscScalar *Pressure;
    PetscScalar *Coriolis_parameter;
    PetscScalar *Latitude;
    Vec          Diabatic_heating;
    Vec          Vorticity_advection;
    Vec          Temperature_advection;
    Vec          Diabatic_heating_tendency;
    Vec          Friction;
    Vec          Friction_u_tendency;
    Vec          Friction_u;
    Vec          Friction_v_tendency;
    Vec          Friction_v;
    Vec          Geopotential_height;
    Vec          Geopotential_height_tendency;
    Vec          Horizontal_wind;
    Vec          omega[NUM_GENERALIZED_OMEGA_COMPONENTS];
    Vec          Total_omega;
    Vec          Temperature;
    Vec          Temperature_tendency;
    Vec          Sigma_parameter;
    Vec          Surface_pressure;
    Vec          Surface_attennuation;
    Vec          Vorticity;
    Vec          Vorticity_tendency;
    Vec          Vorticity_tendency_v;
    Vec          Vorticity_tendency_t;
    Vec          Vorticity_tendency_f;
    Vec          Vorticity_tendency_q;
    Vec          Vorticity_tendency_a;
    Vec          Streamfunction;
    Vec          Ur;
    Vec          Vr;
    Vec          Rotational_wind;
    Vec          Vorticity_advection_by_vr;
};

Context new_context (Options, Files);

Vec new_vec (Context *);

void free_context (Context *ctx);

void update_context (size_t, Files, Context *);

int diabatic_heating (Context *, const int ncid, const int time);
int diabatic_heating_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec Q,
    Vec Qtend, Context *ctx);
int friction (Context *, const int ncid, const int time);
int friction_u_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, const char *varname, Vec F,
    Vec Ftend, Context *ctx);
int friction_v_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, const char *varname, Vec F,
    Vec Ftend, Context *ctx);
int horizontal_wind_and_vorticity_and_vorticity_tendency (
    int ncid, size_t step, size_t first, size_t mt, double *t, DM da, DM da2,
    size_t my, PetscScalar hx, PetscScalar hy, PetscScalar *latitude, Vec V, Vec zeta,
    Vec zetatend, Context *ctx);
int temperature (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec T,
    Vec Ttend, Context *);
int geopotential_height (
    int ncid, size_t step, size_t first, size_t mt, double *t, Vec Z,
    Vec Ztend, Context *);
int sigma_parameter (
    DM da, PetscInt mz, PetscScalar *p, Vec Tvec, Vec sigmavec);
int vorticity_tendency_v (Context *);
int vorticity_tendency_t (Context *);
int vorticity_tendency_f (Context *);
int vorticity_tendency_q (Context *);
int vorticity_tendency_a (Context *);
int calculate_static_stability (Context *, Vec bvec);
