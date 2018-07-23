/*
 * ilc_task.c
 *
 *  Created on: Aug 23, 2017
 *      Author: okoc
 */

// SL general includes of system headers
#include "table_tennis_common.h"
#include "sl_interface.h"

void add_grav_comp_task( void );
static int change_grav_comp_task(void);
static int init_grav_comp_task(void);
static int run_grav_comp_task(void);
static int goto_posture(int posture);
static int goto_center_posture();
static int goto_right_posture();
static int goto_left_posture();
static void turn_off_fb();
static void turn_on_fb();

SL_DJstate init_joint_state[N_DOFS+1];

/*
 * Adds the task to the task menu
 */
void add_grav_comp_task( void ) {
	int i;
	char varname[30];

	addTask("GRAV COMP task", init_grav_comp_task, run_grav_comp_task, change_grav_comp_task);
}

/*
 * Changes task parameters
 */
static int change_grav_comp_task(void) {
	int i,j;
	return TRUE;
}

/*
 * Initialization for task
 *
 */
static int init_grav_comp_task(void) {

	int i;
	int posture, ready, turn_off_pd; // flags

	for (i = 1; i <= N_DOFS; i++) {
		joint_des_state[i].th = joint_state[i].th;
		joint_des_state[i].thd = 0.0;
		joint_des_state[i].thdd = 0.0;
		joint_des_state[i].uff = 0.0;
	}

	turn_on_fb();

	// for real setup
	setDefaultEndeffector();
	double value;
	get_double("Endeffector mass?.\n", 0, &value);
	endeff[RIGHT_HAND].m = value;
	get_double("Endeffector mcm[Z]?.\n", 0, &value);
	endeff[RIGHT_HAND].mcm[_Z_] = value;
	endeff[RIGHT_HAND].x[_Z_] = 0.3;

	/* check whether any other task is running */
	if (strcmp(current_task_name,NO_TASK) != 0) {
		printf("Task can only be run if no other task is running!\n");
		return FALSE;
	}

	/* go to a save posture */
	get_int("Which posture? 0 = CENTRE, 1 = RIGHT, 2 = LEFT.\n", 0, &posture);
	if (!goto_posture(posture))
		return FALSE;

	get_int("Turn off PD? 0 = NO, 1 = YES.\n", 0, &turn_off_pd);

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
	if (turn_off_pd)
		turn_off_fb();
	return TRUE;

}

/*
 * Runs the task from the task servo: REAL TIME requirements!
 *
 */
static int run_grav_comp_task(void) {


	int i = 0;
	for ( i = 1; i <= N_DOFS; i++ ) {
		joint_des_state[i].th = joint_state[i].th;//+joint_state[i].thd/task_servo_rate; //+joint_state[i].thdd/task_servo_rate/task_servo_rate; 
		joint_des_state[i].thd = 0.0;
		joint_des_state[i].thdd = 0.0;
		joint_des_state[i].uff = 0.0;
	}

	//SL_InvDyn( NULL, joint_des_state, endeff, &base_state, &base_orient );
	SL_InvDyn(joint_state, joint_des_state, endeff, &base_state, &base_orient );
	//check_range( joint_des_state );

	save_joint_data(joint_state,joint_des_state,0);
	save_ball_data("tcp://helbe:7660",1); //0 = NO DEBUG

	return TRUE;
}

/*
 * Switch between motor servo feedback and task servo feedback
 * MAYBE DELAYED?
 * Can be used before starting trajectory tracking
 */
static void turn_off_fb() {

	// to change feedback
	static int firsttime = TRUE;
	static Vector zerogain;

	if (firsttime) {
		firsttime = FALSE;
		zerogain = my_vector(1,7);
	}
	changePIDGains(zerogain, zerogain, zerogain);
}

/*
 * Switch between motor servo feedback and task servo feedback
 * MAYBE DELAYED?
 * Can be used before starting trajectory tracking
 */
static void turn_on_fb() {

	static int firsttime = TRUE;

	// to change feedback
	static Vector pgain, dgain, igain;

	if (firsttime) {
		firsttime = FALSE;
		pgain = my_vector(1,7);
		dgain = my_vector(1,7);
		igain = my_vector(1,7);
		read_gains(config_files[GAINS],pgain, dgain, igain);
	}


	changePIDGains(pgain, dgain, igain);
}

/*
 * Initialize robot posture at right, center and left sides
 */
static int goto_posture(int posture) {

	int i;
	if (posture == 1)
		goto_right_posture();
	else if (posture == 0)
		goto_center_posture();
	else if (posture == 2)
		goto_left_posture();
	else { // unrecognized input
		printf("Unrecognized posture. Quitting...\n");
		return FALSE;
	}

	for (i = 1; i <= N_DOFS; i++) {
		joint_des_state[i].th = init_joint_state[i].th;
		joint_des_state[i].thd = 0.0;
		joint_des_state[i].thdd = 0.0;
	}
	return TRUE;
}

/*
 * Robot goes to a center posture
 */
static int goto_center_posture() {

	int i;
	bzero((char *)&(init_joint_state[1]), N_DOFS * sizeof(init_joint_state[1]));
	init_joint_state[1].th = 0.0;
	init_joint_state[2].th = 0.0;
	init_joint_state[3].th = 0.0;
	init_joint_state[4].th = 1.5;
	init_joint_state[5].th = -1.75;
	init_joint_state[6].th = 0.0;
	init_joint_state[7].th = 0.0;

	go_target_wait(init_joint_state);
		return FALSE;

	return TRUE;
}

/*
 * Robot waits on the right-hand side
 */
static int goto_right_posture() {

	int i;
	bzero((char *)&(init_joint_state[1]), N_DOFS * sizeof(init_joint_state[1]));
	init_joint_state[1].th = 1.0;
	init_joint_state[2].th = -0.2;
	init_joint_state[3].th = -0.1;
	init_joint_state[4].th = 1.8;
	init_joint_state[5].th = -1.57;
	init_joint_state[6].th = 0.1;
	init_joint_state[7].th = 0.3;

	go_target_wait(init_joint_state);
	if (!go_target_wait_ID(init_joint_state))
		return FALSE;
	return TRUE;
}

/*
 * Robot waits on the left-hand side
 */
static int goto_left_posture() {

	int i;
	bzero((char *)&(init_joint_state[1]), N_DOFS * sizeof(init_joint_state[1]));
	init_joint_state[1].th = -1.0;
	init_joint_state[2].th = 0.0;
	init_joint_state[3].th = 0.0;
	init_joint_state[4].th = 1.5;
	init_joint_state[5].th = -1.57;
	init_joint_state[6].th = 0.1;
	init_joint_state[7].th = 0.3;

	go_target_wait(init_joint_state);
	if (!go_target_wait_ID(init_joint_state))
		return FALSE;
	return TRUE;
}
