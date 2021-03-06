#include "rules.h"
#include "constants.h"
#include "context.h"
#include "equation.h"
#include "fields.h"
#include "io.h"
#include "omega.h"
#include "operators.h"
#include "ops.h"
#include "targets.h"
#include <petscerror.h>
#include <petscksp.h>
#include <stdbool.h>

static void compute_diabatic_heating (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_v (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_t (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_f (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_q (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_a (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_vr (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_tendency_vd (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vorticity_advection (
    TARGET, Targets *, const Rules *, Context *);
static void compute_temperature_advection (
    TARGET, Targets *, const Rules *, Context *);
static void compute_friction (
    TARGET, Targets *, const Rules *, Context *);
static void compute_total_omega (
    TARGET, Targets *, const Rules *, Context *);
static void compute_horizontal_wind_etc (
    TARGET, Targets *, const Rules *, Context *);
static void compute_omega_operator (
    TARGET, Targets *, const Rules *, Context *);
static void compute_omega_component (
    TARGET, Targets *, const Rules *, Context *);
static void read_streamfunction (
    TARGET, Targets *, const Rules *, Context *); 
static void read_velocity_potential (
    TARGET, Targets *, const Rules *, Context *); 
static void compute_vadvr (
    TARGET, Targets *, const Rules *, Context *);
static void compute_vadvd (
    TARGET, Targets *, const Rules *, Context *);
static void compute_temperature_and_tendency (
    TARGET, Targets *, const Rules *, Context *);
static void compute_geopotential_height_and_tendency (
    TARGET, Targets *, const Rules *, Context *);
static void compute_sigma_parameter (
    TARGET, Targets *, const Rules *, Context *);
static void compute_surface_attennuation (
    TARGET, Targets *, const Rules *, Context *);
static void read_field_2d (TARGET, Targets *, const Rules *, Context *);
static void read_field_3d (TARGET, Targets *, const Rules *, Context *);

Rules new_rules (void) {
    Rules rules = {{
            [TARGET_FIELD_DIABATIC_HEATING] =
                (Rule){.prerequisites = 0,
                       .recipe        = 0},

            [TARGET_FIELD_VORTICITY_ADVECTION] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_HORIZONTAL_WIND),
                   .recipe = compute_vorticity_advection},

            [TARGET_FIELD_TEMPERATURE_ADVECTION] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_TEMPERATURE,
                    TARGET_FIELD_HORIZONTAL_WIND),
                   .recipe = compute_temperature_advection},

            [TARGET_FIELD_DIABATIC_HEATING_TENDENCY] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_DIABATIC_HEATING),
                   .recipe = compute_diabatic_heating},

            [TARGET_FIELD_FRICTION] =
                (Rule){.prerequisites = 0,
                       .recipe        = compute_friction},
/*
            [TARGET_FIELD_GEOPOTENTIAL_HEIGHT] =
                (Rule){.prerequisites = 0,
                       .recipe = 0},//compute_geopotential_height_and_tendency},

            [TARGET_FIELD_GEOPOTENTIAL_HEIGHT_TENDENCY] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_GEOPOTENTIAL_HEIGHT),
                   .recipe = 0},
*/
            [TARGET_FIELD_HORIZONTAL_WIND] =
                (Rule){.prerequisites = 0,
                       .recipe        = compute_horizontal_wind_etc},

            [TARGET_FIELD_OMEGA_OPERATOR] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_HORIZONTAL_WIND),
                   .recipe = compute_omega_operator},

            [TARGET_FIELD_OMEGA_V] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_VORTICITY_ADVECTION,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                   .recipe = compute_omega_component},

            [TARGET_FIELD_OMEGA_T] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_TEMPERATURE_ADVECTION,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                   .recipe = compute_omega_component},

            [TARGET_FIELD_OMEGA_Q] =
                (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_DIABATIC_HEATING_TENDENCY,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
		       .recipe = compute_omega_component},

            [TARGET_FIELD_OMEGA_F] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_FRICTION_U_TENDENCY,
                    TARGET_FIELD_FRICTION_V_TENDENCY,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                   .recipe = compute_omega_component},

            [TARGET_FIELD_OMEGA_A] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_TEMPERATURE,
                    TARGET_FIELD_VORTICITY_TENDENCY,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                   .recipe = compute_omega_component},

            [TARGET_FIELD_TOTAL_OMEGA] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_V,
                    TARGET_FIELD_OMEGA_T,
                    TARGET_FIELD_OMEGA_Q,
                    TARGET_FIELD_OMEGA_F,
                    TARGET_FIELD_OMEGA_A,),
                   .recipe = compute_total_omega},

            [TARGET_FIELD_TEMPERATURE] =
                (Rule){.prerequisites = 0,
                       .recipe        = compute_temperature_and_tendency},

            [TARGET_FIELD_TEMPERATURE_TENDENCY] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_TEMPERATURE),
                       .recipe = 0},

            [TARGET_FIELD_SIGMA_PARAMETER] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_TEMPERATURE),
                       .recipe = compute_sigma_parameter},

            [TARGET_FIELD_SURFACE_PRESSURE] =
                (Rule){.prerequisites = 0,
                       .recipe = read_field_2d},

            [TARGET_FIELD_SURFACE_ATTENNUATION] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_SURFACE_PRESSURE),
                       .recipe = compute_surface_attennuation},

            [TARGET_FIELD_VORTICITY] = (Rule){.prerequisites = new_target_list (
                                                  TARGET_FIELD_HORIZONTAL_WIND),
                                              .recipe = 0},

            [TARGET_FIELD_VORTICITY_TENDENCY] =
                (Rule){.prerequisites =
                           new_target_list (TARGET_FIELD_HORIZONTAL_WIND),
                       .recipe = 0},

            [TARGET_FIELD_VORTICITY_TENDENCY_V] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_VORTICITY_ADVECTION,
                    TARGET_FIELD_OMEGA_V),
                   .recipe = compute_vorticity_tendency_v},

            [TARGET_FIELD_VORTICITY_TENDENCY_T] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_OMEGA_T,
                    TARGET_FIELD_TEMPERATURE_ADVECTION),
                   .recipe = compute_vorticity_tendency_t},

            [TARGET_FIELD_VORTICITY_TENDENCY_F] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_OMEGA_F,
                    TARGET_FIELD_FRICTION),
                   .recipe = compute_vorticity_tendency_f},

            [TARGET_FIELD_VORTICITY_TENDENCY_Q] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_OMEGA_Q,
                    TARGET_FIELD_DIABATIC_HEATING_TENDENCY),
                   .recipe = compute_vorticity_tendency_q},

            [TARGET_FIELD_VORTICITY_TENDENCY_A] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_OMEGA_A),
                   .recipe = compute_vorticity_tendency_a},

            [TARGET_FIELD_FRICTION_U_TENDENCY] =
            (Rule){.prerequisites = new_target_list(TARGET_FIELD_FRICTION),
                   .recipe = 0},

            [TARGET_FIELD_FRICTION_V_TENDENCY] =
            (Rule){.prerequisites = new_target_list(TARGET_FIELD_FRICTION),
                   .recipe = 0},

            [TARGET_FIELD_STREAMFUNCTION] =
            (Rule){.prerequisites = new_target_list(TARGET_FIELD_VORTICITY),
                   .recipe = read_streamfunction},

            [TARGET_FIELD_VELOCITY_POTENTIAL] =
            (Rule){.prerequisites = new_target_list(TARGET_FIELD_VORTICITY),
                   .recipe = read_velocity_potential},

            [TARGET_FIELD_VORTICITY_ADVECTION_BY_VR] =
            (Rule){.prerequisites = new_target_list(
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_STREAMFUNCTION),
                   .recipe = compute_vadvr},

            [TARGET_FIELD_VORTICITY_ADVECTION_BY_VD] =
            (Rule){.prerequisites = new_target_list(
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_VELOCITY_POTENTIAL),
                   .recipe = compute_vadvd},

            [TARGET_FIELD_OMEGA_VR] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_VORTICITY_ADVECTION_BY_VR,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                   .recipe = compute_omega_component},

            [TARGET_FIELD_OMEGA_VD] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_OMEGA_OPERATOR,
                    TARGET_FIELD_VORTICITY_ADVECTION_BY_VD,
                    TARGET_FIELD_SURFACE_ATTENNUATION,
                    TARGET_FIELD_SIGMA_PARAMETER,
                    TARGET_FIELD_VORTICITY),
                   .recipe = compute_omega_component},

            [TARGET_FIELD_VORTICITY_TENDENCY_VR] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_VORTICITY_ADVECTION_BY_VR,
                    TARGET_FIELD_OMEGA_VR),
                   .recipe = compute_vorticity_tendency_vr},

            [TARGET_FIELD_VORTICITY_TENDENCY_VD] =
            (Rule){.prerequisites = new_target_list (
                    TARGET_FIELD_HORIZONTAL_WIND,
                    TARGET_FIELD_VORTICITY,
                    TARGET_FIELD_VORTICITY_ADVECTION_BY_VD,
                    TARGET_FIELD_OMEGA_VD),
                   .recipe = compute_vorticity_tendency_vd},
    }};

    return rules;
}

static char *target_name (TARGET id, const Targets *targets) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD:
        return (char *)&targets->target[id].field.name;
        break;
    case TARGET_TYPE_OPERATOR:
        return (char *)&targets->target[id].operator.name ;
        break;
    default:
        info ("ERROR: Not implemented in target_name()");
        return 0;
    }
}

void draw (const Rules *rules, const Targets *targets, const char *fname) {
    FILE *fd = fopen (fname, "w");
    PetscPrintf (
        PETSC_COMM_WORLD, "Writing target dependencies to file %s\n\n", fname);
    PetscFPrintf (PETSC_COMM_WORLD, fd, "digraph Rules {\n");
    for (size_t i = 0; i < NUM_TARGET; i++) {
        TARGETS *prereq = rules->rule[i].prerequisites;
        while (prereq) {
            PetscFPrintf (
                PETSC_COMM_WORLD, fd, "  %s -> %s\n",
                target_name (prereq->this, targets), target_name (i, targets));
            prereq = prereq->next;
        }
    }
    PetscFPrintf (PETSC_COMM_WORLD, fd, "}\n");
    fclose (fd);
}

/* Constructs a list of rules that can be executed in parallel.
 *
 * Rule is eligible for execution if
 *   (1) all prerequisites are already updated to next step, and that
 *   (2) none of the recipes of the children are pending on the current step,
 *       and
 *   (3) the rule is not at the last time step (finished).
 */

bool more_todo (
    const Rules *rules, Targets *targets, TARGETS *todo[], Context *ctx) {
    bool eligible[NUM_TARGET];
    for (size_t i = 0; i < NUM_TARGET; i++) {
        eligible[i] = true;
    }
    for (size_t i = 0; i < NUM_TARGET; i++) {
        TARGETS *prereq = rules->rule[i].prerequisites;
        while (prereq) {
            if (targets->target[i].time !=
                targets->target[prereq->this].time - 1) {    // (1)
                eligible[i] = false;
            }
            if (targets->target[i].time !=
                targets->target[prereq->this].time) {    // (2)
                eligible[prereq->this] = false;
            }
            prereq = prereq->next;
        }
        if (targets->target[i].time == ctx->last) {    // (3)
            eligible[i] = false;
        }
    }
    for (size_t i = 0; i < NUM_TARGET; i++) {
        if (eligible[i]) {
            *todo = push (i, *todo);
        }
    }
    return (*todo ? true : false);
}

void print_target_list (
    const char *title, TARGETS *head, const Targets *targets) {
    PetscPrintf (PETSC_COMM_WORLD, "%s: ", title);
    while (head) {
        info (
            "%s[%zu]%s", target_name (head->this, targets),
            targets->target[head->this].time + 1, (head->next ? ", " : ""));
        head = head->next;
    }
    info ("\n");
}

static void compute_vorticity_tendency_v (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_v (ctx, ctx->omega[GENERALIZED_OMEGA_COMPONENT_V], ctx->Horizontal_wind, 
                          ctx->Vorticity_advection, ctx->Vorticity_tendency_v);
}

static void compute_vorticity_tendency_t (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_t (ctx);
}

static void compute_vorticity_tendency_f (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_f (ctx);
}

static void compute_vorticity_tendency_q (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_q (ctx);
}

static void compute_vorticity_tendency_a (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_a (ctx);
}

static void compute_vorticity_tendency_vr (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_v (ctx, ctx->omega_vr, ctx->Rotational_wind, ctx->Vorticity_advection_by_vr, ctx->Vorticity_tendency_vr);
}

static void compute_vorticity_tendency_vd (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    vorticity_tendency_v (ctx, ctx->omega_vd, ctx->Divergent_wind, ctx->Vorticity_advection_by_vd, ctx->Vorticity_tendency_vd);
}

static void compute_diabatic_heating (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    diabatic_heating (ctx, ctx->ncid, targets->target[id].time);
}

static void compute_friction (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    friction (ctx, ctx->ncid, targets->target[id].time);
}

static void compute_total_omega (TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    VecZeroEntries(ctx->Total_omega);
    VecAXPY(ctx->Total_omega,1.0,ctx->omega[GENERALIZED_OMEGA_COMPONENT_V]);
    VecAXPY(ctx->Total_omega,1.0,ctx->omega[GENERALIZED_OMEGA_COMPONENT_T]);
    VecAXPY(ctx->Total_omega,1.0,ctx->omega[GENERALIZED_OMEGA_COMPONENT_F]);
    VecAXPY(ctx->Total_omega,1.0,ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]);
    VecAXPY(ctx->Total_omega,1.0,ctx->omega[GENERALIZED_OMEGA_COMPONENT_A]);
}

static void compute_omega_operator (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    KSPSetComputeOperators (ctx->ksp, omega_compute_operator, ctx);
}

static void compute_omega_component (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
//    KSPSetComputeOperators (ctx->ksp, omega_compute_operator, ctx);

    switch (id) {
    case TARGET_FIELD_OMEGA_V: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_V, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega[GENERALIZED_OMEGA_COMPONENT_V]);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    case TARGET_FIELD_OMEGA_T: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_T, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega[GENERALIZED_OMEGA_COMPONENT_T]);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    case TARGET_FIELD_OMEGA_Q: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_Q, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    case TARGET_FIELD_OMEGA_F: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_F, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega[GENERALIZED_OMEGA_COMPONENT_F]);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    case TARGET_FIELD_OMEGA_A: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_A, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega[GENERALIZED_OMEGA_COMPONENT_A]);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    case TARGET_FIELD_OMEGA_VR: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_Vr, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega_vr);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    case TARGET_FIELD_OMEGA_VD: {
        KSPSetComputeRHS (ctx->ksp, omega_compute_rhs_F_Vd, ctx);
        KSPSolve (ctx->ksp, 0, ctx->omega_vd);
        // KSPGetSolution (ctx->ksp, &x);
        break;
    }
    default: { info ("Not implemented in compute_omega_component.\n"); }
    }
}

static void read_streamfunction (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx)  {
    file_read_3d (ctx->ncid, targets->target[id].time, "STRF", ctx->Streamfunction);
}

static void read_velocity_potential (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx)  {
    file_read_3d (ctx->ncid, targets->target[id].time, "VPOT", ctx->Velocity_potential);
}


static void compute_vadvd (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    Vec          vpot = ctx->Velocity_potential;
    Vec          Vd = ctx->Divergent_wind;
    Vec          zeta = ctx->Vorticity;
    PetscScalar* f    = ctx->Coriolis_parameter;
    Vec          vadvd = ctx->Vorticity_advection_by_vd;
    Vec          b,c;
    const double r_inv = 1.0 / earth_radius;

    // Calculate u and v winds from the derivatives of velocity potential
    DMGetGlobalVector (ctx->da, &b);
    VecCopy (vpot, b);    
    xder (b, ctx);
    VecScale(b,r_inv);
    DMGetGlobalVector (ctx->da, &c);
    VecCopy (vpot, c);    
    yder (c, ctx);
    VecScale(c,r_inv);

    // Insert the wind components to total wind vector
    VecStrideScatter (b, 0, Vd, INSERT_VALUES);
    VecStrideScatter (c, 1, Vd, INSERT_VALUES);

    // calculate the advection of total vorticity
    DMGetGlobalVector (ctx->da, &b);
    VecCopy (zeta, b);

    field_array1d_add (b, f, DMDA_Y);

    horizontal_advection (b, Vd, ctx);

    // scale by -1 
    VecCopy (b, vadvd);
    VecScale(vadvd,-1.0);
    DMRestoreGlobalVector (ctx->da, &b);
    DMRestoreGlobalVector (ctx->da, &c);
}

static void compute_vadvr (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    Vec          strf = ctx->Streamfunction;
    Vec          Vr =ctx->Rotational_wind;
    Vec          zeta = ctx->Vorticity;
    PetscScalar* f    = ctx->Coriolis_parameter;
    Vec          vadvr = ctx->Vorticity_advection_by_vr;
    Vec          b,c;
    const double r_inv = 1.0 / earth_radius;

    // Calculate u and v winds from the derivatives of streamfunction
    DMGetGlobalVector (ctx->da, &b);
    VecCopy (strf, b);    
    yder (b, ctx);
    VecScale(b,-1*r_inv);
    DMGetGlobalVector (ctx->da, &c);
    VecCopy (strf, c);    
    xder (c, ctx);
    VecScale(c,r_inv);

    // Insert the wind components to total wind vector
    VecStrideScatter (b, 0, Vr, INSERT_VALUES);
    VecStrideScatter (c, 1, Vr, INSERT_VALUES);

    // calculate the advection of total vorticity
    DMGetGlobalVector (ctx->da, &b);
    VecCopy (zeta, b);

    field_array1d_add (b, f, DMDA_Y);

    horizontal_advection (b, Vr, ctx);

    // scale by -1 
    VecCopy (b, vadvr);
    VecScale(vadvr,-1.0);
    DMRestoreGlobalVector (ctx->da, &b);
    DMRestoreGlobalVector (ctx->da, &c);
}

static void compute_horizontal_wind_etc (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    horizontal_wind_and_vorticity_and_vorticity_tendency (
        ctx->ncid, targets->target[id].time, ctx->first, ctx->mt,
        ctx->Time_coordinate, ctx->da, ctx->da2, ctx->my, ctx->hx, ctx->hy,
        ctx->Latitude,ctx->Horizontal_wind, ctx->Vorticity, ctx->Vorticity_tendency, ctx);
}

static void compute_temperature_and_tendency (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    static Vec Tnext = 0;

    if (!Tnext) {    // The first step
        VecDuplicate (ctx->Temperature, &Tnext);
    }

    temperature (
        ctx->ncid, targets->target[id].time, ctx->first, ctx->mt,
        ctx->Time_coordinate, ctx->Temperature, ctx->Temperature_tendency, ctx);

    if (targets->target[id].time == ctx->last) {
        VecDestroy (&Tnext);
    }
}

static void compute_geopotential_height_and_tendency (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    static Vec Znext = 0;

    if (!Znext) {    // The first step
        VecDuplicate (ctx->Geopotential_height, &Znext);
    }

    geopotential_height (
        ctx->ncid, targets->target[id].time, ctx->first, ctx->mt,
        ctx->Time_coordinate, ctx->Geopotential_height, ctx->Geopotential_height_tendency, ctx);

    if (targets->target[id].time == ctx->last) {
        VecDestroy (&Znext);
    }
}

static void compute_sigma_parameter (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    sigma_parameter (
        ctx->da, ctx->mz, ctx->Pressure, ctx->Temperature,
        ctx->Sigma_parameter);
}

static void compute_surface_attennuation (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    mul_fact (ctx, ctx->Surface_attennuation);
}

static void compute_vorticity_advection (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    Vec          V    = ctx->Horizontal_wind;
    Vec          zeta = ctx->Vorticity;
    PetscScalar* f    = ctx->Coriolis_parameter;
    Vec          vadv = ctx->Vorticity_advection;
    Vec          b;

    DMGetGlobalVector (ctx->da, &b);
    VecCopy (zeta, b);

    field_array1d_add (b, f, DMDA_Y);

    horizontal_advection (b, V, ctx);
    VecCopy (b, vadv);
    VecScale(vadv,-1.0);
    DMRestoreGlobalVector (ctx->da, &b);
}

static void compute_temperature_advection (
    TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    Vec          V    = ctx->Horizontal_wind;
    Vec          T = ctx->Temperature;
    Vec          tadv = ctx->Temperature_advection;
    Vec          b;

    DMGetGlobalVector (ctx->da, &b);
    VecCopy (T, b);

    horizontal_advection (b, V, ctx);
    VecCopy (b, tadv);
    VecScale(tadv,-1.0);
    DMRestoreGlobalVector (ctx->da, &b);
}

static void
read_field_2d (TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    PetscScalar **psfcArray;
    Vec          psfc     = ctx->Surface_pressure;
    PetscInt      ys, xs, ym, xm;
    int            idd;

//    read2D (ctx->ncid, targets->target[id].time, "SP",psfc);

    DMDAGetCorners (ctx->daxy, &xs, &ys, 0, &xm, &ym, 0);
    DMDAVecGetArray (ctx->daxy, psfc, &psfcArray);
    size_t start[4] = {targets->target[id].time, 0 , ys, xs};
    size_t count[4] = {1, 1 , ym, xm};
    nc_inq_varid (ctx->ncid, targets->target[id].field.name, &idd);
    nc_get_vara_double (ctx->ncid, idd, start, count, &psfcArray[start[2]][start[3]]);
    DMDAVecRestoreArray (ctx->daxy, psfc, &psfcArray);

}

static void
read_field_3d (TARGET id, Targets *targets, const Rules *rules, Context *ctx) {
    file_read_3d (
        ctx->ncid, targets->target[id].time, targets->target[id].field.name,
        targets->target[id].field.vec);
}

/*static void
read_target (TARGET id, size_t time, Targets *targets, Context *ctx) {
    switch (targets->target[id].type) {
    case TARGET_TYPE_FIELD: {
        Field *f = &targets->target[id].field;
        file_read_3d (ctx->ncid, time, f->name, f->vec);
        break;
    }
    case TARGET_TYPE_OPERATOR: {
        info ("Function read() is not implemented for TARGET_TYPE_OPERATOR.\n");
        break;
    }
    default:
        info ("Internal error in rules.c, function read().\n");
    }
}


void recipe (RULE id, Rules *rules) {
    Rule *field = &rules->field[id];
    field->step++;
    if (field->ncid_in)
        read_field (id, rules);
    if (field->recipe)
        field->recipe (id, rules);
}
*/

void update (TARGET id, const Rules *rules, Targets *targets, Context *ctx) {
    Target *    target = &targets->target[id];
    const Rule *rule   = &rules->rule[id];
    target->time++;
    switch (target->type) {
    case TARGET_TYPE_FIELD: {
        if (rule->recipe) {
            info ("Computing %s[%zu]\n", target->field.name, target->time);
            rules->rule[id].recipe (id, targets, rules, ctx);
        }
        if (target->field.write) {
            write3D (
                ctx->ncid, target->time, target->field.name, target->field.vec);
        }
        break;
    }
    default:
        info ("Update of TARGET_TYPE_OPERATOR not implemented yet.\n");
    }
}

void run (const Rules *rules, Targets *targets, Context *ctx) {
    TARGETS *todo = 0;
    while ((more_todo (rules, targets, &todo, ctx))) {
        while (todo) {
            print_target_list ("todo", todo, targets);
            TARGETS *head = pop (&todo);
            update (head->this, rules, targets, ctx);
            free (head);
        }
    }
}
