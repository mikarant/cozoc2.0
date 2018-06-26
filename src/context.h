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
    Vec          Diabatic_heating_forcing;
    Vec          Vorticity_advection_forcing;
    Vec          Vorticity_advection;
    Vec          Temperature_advection_forcing;
    Vec          Temperature_advection;
    Vec          Diabatic_heating_tendency;
    Vec          Friction;
    Vec          Friction_u_tendency;
    Vec          Friction_u;
    Vec          Friction_v_tendency;
    Vec          Friction_v;
    Vec          Geopotential_height;
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
int sigma_parameter (
    DM da, PetscInt mz, PetscScalar *p, Vec Tvec, Vec sigmavec);
