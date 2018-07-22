#ifndef SL_STRUCT_H
#define SL_STRUCT_H

#include "constants.h"

/* The data structures from SL */
/** @brief (actual) joint space state for each DOF */
struct SL_Jstate {
    double   th;   /*!< theta */
    double   thd;  /*!< theta-dot */
    double   thdd; /*!< theta-dot-dot */
    double   ufb;  /*!< feedback portion of command */
    double   u;    /*!< torque command */
    double   load; /*!< sensed torque */
};

/** @brief (desired) joint space state commands for each DOF */
struct SL_DJstate { /*!< desired values for controller */
    double   th;   /*!< theta */
    double   thd;  /*!< theta-dot */
    double   thdd; /*!< theta-dot-dot */
    double   uff;  /*!< feedforward torque command */
    double   uex;  /*!< externally imposed torque */
};

/** @brief (actual) Cartesian state */
struct SL_Cstate {
    double x[const_tt::NCART+1];    /*!< Position [x,y,z] */
    double xd[const_tt::NCART+1];   /*!< Velocity */
    double xdd[const_tt::NCART+1];  /*!< Acceleration */
};

/** @brief Vision blob info coming from SL (after calibration). */
struct blob_state {
    int status = 0; //!< was ball detected reliably in cameras
    double pos[const_tt::NCART] = {0.0,0.0,0.0}; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
};

#endif // SL_STRUCT_H
