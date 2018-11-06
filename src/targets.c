#include "targets.h"
#include "context.h"
#include "fields.h"
#include "netcdf.h"
#include "operators.h"
#include "options.h"
#include <stdbool.h>

Targets new_targets (Options options, Files files, Context *ctx) {
    Targets targets = {
        .target = {

                [TARGET_FIELD_DIABATIC_HEATING] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "Q",
                                              .description = "Diabatic heating",
                                              .units       = "K s**-1",
                                              .vec = ctx->Diabatic_heating},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_ADVECTION] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){.write = false,
                            .name  = "vadv",
                            .description =
                            "Vorticity advection",
                            .units = 0,
                            .vec   = ctx->Vorticity_advection},
                    .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE_ADVECTION] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){.write = false,
                            .name  = "tadv",
                            .description =
                            "Forcing due to temperature advection",
                            .units = 0,
                            .vec   = ctx->Temperature_advection},
                    .time = options.first - 1},

                [TARGET_FIELD_DIABATIC_HEATING_TENDENCY] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){.write = false,
                            .name  = "Qtend",
                            .description =
                            "Diabatic heating from accumulated field",
                            .units = 0,
                            .vec   = ctx->Diabatic_heating_tendency},
                    .time = options.first - 1},

                [TARGET_FIELD_FRICTION] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "F",
                                              .description = "Friction",
                                              .units       = 0,
                                              .vec         = ctx->Friction},
                             .time = options.first - 1},

                [TARGET_FIELD_HORIZONTAL_WIND] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "V",
                                              .description = "Horizontal wind",
                                              .units       = 0,
                                              .vec = ctx->Horizontal_wind},
                             .time = options.first - 1},

                [TARGET_FIELD_OMEGA_OPERATOR] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){
                        .write       = false,
                        .name        = "omega_operator",
                        .description = "",
                        .units       = "",
                        .vec = 0 },
                    .time = options.first - 1},

                [TARGET_FIELD_OMEGA_V] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){
                        .write       = true,
                        .name        = "cozoc_ome_v",
                        .description = "Omega due to vorticity advection",
                        .units       = "Pa s-1",
                        .vec =
                        ctx->omega[GENERALIZED_OMEGA_COMPONENT_V]},
                    .time = options.first - 1},

                [TARGET_FIELD_OMEGA_T] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){
                        .write       = true,
                        .name        = "cozoc_ome_t",
                        .description = "Omega due to temperature advection",
                        .units       = "Pa s-1",
                        .vec =
                        ctx->omega[GENERALIZED_OMEGA_COMPONENT_T]},
                    .time = options.first - 1},

                [TARGET_FIELD_OMEGA_Q] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){
                                .write       = true,
                                .name        = "cozoc_ome_q",
                                .description = "Omega due to diabatic heating",
                                .units       = "Pa s-1",
                                .vec =
                                    ctx->omega[GENERALIZED_OMEGA_COMPONENT_Q]},
                        .time = options.first - 1},

                [TARGET_FIELD_OMEGA_F] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){
                        .write       = true,
                        .name        = "cozoc_ome_f",
                        .description = "Omega due to friction",
                        .units       = "Pa s-1",
                        .vec =
                        ctx->omega[GENERALIZED_OMEGA_COMPONENT_F]},
                    .time = options.first - 1},

                [TARGET_FIELD_OMEGA_A] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){
                        .write       = true,
                        .name        = "cozoc_ome_a",
                        .description = "Omega due to imbalance term",
                        .units       = "Pa s-1",
                        .vec =
                        ctx->omega[GENERALIZED_OMEGA_COMPONENT_A]},
                    .time = options.first - 1},

                [TARGET_FIELD_TOTAL_OMEGA] =
                (Target){
                    .type = TARGET_TYPE_FIELD,
                    .field =
                    (Field){
                        .write       = false,
                        .name        = "cozoc_ome_tot",
                        .description = "Total calculated omega",
                        .units       = "Pa s-1",
                        .vec =
                        ctx->Total_omega},
                    .time = options.first - 1},

                [TARGET_FIELD_SURFACE_ATTENNUATION] =
                    (Target){
                        .type = TARGET_TYPE_FIELD,
                        .field =
                            (Field){.write = false,
                                    .name  = "Attenuation",
                                    .description =
                                        "Pressure level surface attennuation",
                                    .units = 0,
                                    .vec   = ctx->Surface_attennuation},
                        .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "TT",
                                              .description = "Temperature",
                                              .units       = "K",
                                              .vec         = 0},
                             .time = options.first - 1},

                [TARGET_FIELD_TEMPERATURE_TENDENCY] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.write       = false,
                                         .name        = "Ttend",
                                         .description = "Temperature tendency",
                                         .units       = "K -s",
                                         .vec         = ctx->Temperature_tendency},
                             .time = options.first - 1},

                [TARGET_FIELD_SIGMA_PARAMETER] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "sigma",
                                              .description = "Sigma parameter",
                                              .units       = "",
                                              .vec = ctx->Sigma_parameter},
                             .time = options.first - 1},

                [TARGET_FIELD_SURFACE_PRESSURE] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "SP",
                                              .description = "Surface pressure",
                                              .units       = "",
                                              .vec = ctx->Surface_pressure},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY] =
                    (Target){.type  = TARGET_TYPE_FIELD,
                             .field = (Field){.write       = false,
                                              .name        = "zeta",
                                              .description = "Vorticity",
                                              .units       = "",
                                              .vec         = ctx->Vorticity},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY] =
                    (Target){.type = TARGET_TYPE_FIELD,
                             .field =
                                 (Field){.write       = false,
                                         .name        = "zetatend",
                                         .description = "Vorticity tendency",
                                         .units       = "",
                                         .vec = ctx->Vorticity_tendency},
                             .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY_V] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = true,
                                 .name        = "vtend_v",
                                 .description = "Vorticity tendency due to vorticity advection",
                                 .units       = "",
                                 .vec = ctx->Vorticity_tendency_v},
                         .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY_T] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = true,
                                 .name        = "vtend_t",
                                 .description = "Vorticity tendency due to temperature advection",
                                 .units       = "",
                                 .vec = ctx->Vorticity_tendency_t},
                         .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY_F] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = true,
                                 .name        = "vtend_f",
                                 .description = "Vorticity tendency due to friction",
                                 .units       = "",
                                 .vec = ctx->Vorticity_tendency_f},
                         .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY_Q] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = true,
                                 .name        = "vtend_q",
                                 .description = "Vorticity tendency due to diabatic heating",
                                 .units       = "",
                                 .vec = ctx->Vorticity_tendency_q},
                         .time = options.first - 1},

                [TARGET_FIELD_VORTICITY_TENDENCY_A] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = true,
                                 .name        = "vtend_a",
                                 .description = "Vorticity tendency due to imbalance term",
                                 .units       = "",
                                 .vec = ctx->Vorticity_tendency_a},
                         .time = options.first - 1},

                [TARGET_FIELD_FRICTION_U_TENDENCY] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = false,
                                 .name        = "FUtend",
                                 .description = "Friction u tendency",
                                 .units       = "",
                                 .vec = ctx->Friction_u_tendency},
                         .time = options.first - 1},

                [TARGET_FIELD_FRICTION_V_TENDENCY] =
                (Target){.type = TARGET_TYPE_FIELD,
                         .field =
                         (Field){.write       = false,
                                 .name        = "FVtend",
                                 .description = "Friction v tendency",
                                 .units       = "",
                                 .vec = ctx->Friction_v_tendency},
                         .time = options.first - 1},

        }};

    nc_redef (files.ncid_out);
    for (size_t i = 0; i < NUM_TARGET; i++) {
        Target *t = &targets.target[i];
        if (t->type == TARGET_TYPE_FIELD && t->field.write) {
            file_def_var (files.ncid_out, t->field.name, &files);
        }
    }
    nc_enddef (files.ncid_out);

    return targets;
}

TARGETS *push (TARGET target, TARGETS *oldhead) {
    TARGETS *newhead = (TARGETS *)malloc (sizeof (TARGETS));
    newhead->this    = target;
    newhead->next    = oldhead;
    return newhead;
}

TARGETS *pop (TARGETS **head) {
    TARGETS *node = *head;
    if (node) {
        *head      = node->next;
        node->next = 0;
    }
    return node;
}

TARGETS *_new_target_list (const size_t n, const TARGET f[]) {
    TARGETS *p = 0;
    for (size_t i = 0; i < n; i++) {
        p = push (f[i], p);
    }
    return p;
}
