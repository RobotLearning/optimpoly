/*
 * dmp_task.c
 *
 *  Created on: June 10, 2018
 *      Author: okoc
 */

/*
 * SERVE BALL WITH A DMP
 *  Created on: June 10, 2018
 *      Author: okoc
 */

#include "table_tennis_common.h"
#include "sl_interface.h"

void add_serve_with_dmp_task( void );
static int change_serve_with_dmp_task(void);
static int init_serve_with_dmp_task(void);
static int run_serve_with_dmp_task(void);
static void compute_torques();
static void check_safety();
static int goto_custom_posture(double custom_pose[]);

int init_dmp = FALSE;

/*
 * Adds the task to the task menu
 */
void add_serve_with_dmp_task( void ) {
	int i;
	char varname[30];

	addTask("SERVE WITH DMP task", init_serve_with_dmp_task, run_serve_with_dmp_task, change_serve_with_dmp_task);
}

/*
 * Changes task parameters
 */
static int change_serve_with_dmp_task(void) {
	int i,j;
	return TRUE;
}

/*
 * Initialization for task
 *
 */
static int init_serve_with_dmp_task(void) {

	int ready; // flags
	double custom_pose[N_DOFS] = {0.0};

	/* check whether any other task is running */
	if (strcmp(current_task_name,NO_TASK) != 0) {
		printf("Task can only be run if no other task is running!\n");
		return FALSE;
	}

	// init dmp
	init_dmp_serve(custom_pose,&init_dmp);

	/* go to a save posture */
	goto_custom_posture(custom_pose);
	//if (!goto_custom_posture(custom_pose))
	//	return FALSE;

	// for real setup
	setDefaultEndeffector();
	endeff[RIGHT_HAND].x[_Z_] = .3;

	/* ready to go */
	ready = 999;
	while (ready == 999) {
		if (!get_int("Enter 1 to start or anything else to abort ...",ready, &ready))
			return FALSE;
	}
	if (ready != 1)
		return FALSE;

	// turn on/off real time
	changeRealTime(TRUE);
	return TRUE;

}

/*
 * Runs the task from the task servo: REAL TIME requirements!
 *
 */
static int run_serve_with_dmp_task(void) {

    // serve ball with a feed-forward DMP policy
	serve_with_dmp(joint_state,joint_des_state,&init_dmp);

	// compute torques based on inverse dynamics
	compute_torques();

	return TRUE;
}

/*
 * Check the safety of the desired joint states explicitly
 * and then calculate the u_ff with inverse dynamics
 *
 * If friction compensation is turned on, then add some compensation on top of u_ff
 *
 */
static void compute_torques() {

	check_safety();

	// control the robot
	// calculate the feedforward commands with inverse dynamics
	SL_InvDyn(NULL, joint_des_state, endeff, &base_state, &base_orient);
	/*if (friction_comp) {
		addFrictionModel(joint_des_state);
	}*/

}

/*
 * Checks the safety of the calculated desired state
 * explicitly
 *
 */
static void check_safety() {

	if (!check_joint_limits(joint_state, SLACK)) {
		printf("Joint limits are exceeding limits! Freezing...\n");
		freeze();
	}
	if (!check_des_joint_limits(joint_des_state, SLACK)) {
		printf("Joint des limits are exceeding limits! Freezing...\n");
		freeze();
	}
	/*if(!check_range(joint_des_state)) {
		printf("Exceeding torque limits! Freezing...\n");
		freeze();
	}*/
	if (collision_detection(racket_state)) {
		printf("Collision with table detected!\n");
		freeze();
	}
}

/*
 * Robot goes to a custom posture
 */
static int goto_custom_posture(double custom_pose[N_DOFS]) {

    int i;
    int move = 0;
    SL_DJstate init_joint_state[N_DOFS+1];
    bzero((char *)&(init_joint_state[1]), N_DOFS * sizeof(init_joint_state[1]));
    printf("Min joint lim: [");
    for (i = 0; i < N_DOFS; i++) {
        printf("%.3f ",joint_range[i+1][MIN_THETA]);
    }
    printf("\nCustom pose: [");
    for (i = 0; i < N_DOFS; i++) {
        printf("%.3f ",custom_pose[i]);
        init_joint_state[i+1].th = custom_pose[i];
    }
    printf("\nMax joint lim: [");
    for (i = 0; i < N_DOFS; i++) {
        printf("%.3f ",joint_range[i+1][MAX_THETA]);
    }
    printf("\n");
    get_int("Do you want to move to json startup posture? 1 = YES, 0 = NO ...",move, &move);
    if (move) {
        go_target_wait(init_joint_state);
        if (!go_target_wait_ID(init_joint_state))
            return FALSE;
    }
    return TRUE;
}
