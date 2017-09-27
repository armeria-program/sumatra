//
//  sca.h
//  Sumatra
//
//  Created by Gokhan Selamet on 19/06/17.
//  Copyright © 2017 Ferhat Ayaz. All rights reserved.
//

#ifndef sca_h
#define sca_h

#include <stdio.h>
#include "vec3.h"
#include "smtr.h"

struct simstep_interval {
	void *DATA;
	int interval;
};

// -------------------------------------------
// Event Constractors
// -------------------------------------------
struct event_data{
	int (*condition_function)(void *);
	void *condition_data;
	int (*action_function)(void *);
	void *action_data;
};

typedef struct event_data eventData;

eventData* sca_init_eventData ();
int sca_event (eventData *E);
int sca_condition_simstep_interval (simtime *intervaltime);
int sca_condition_simstep_limit (simtime *limit);
int sca_condition_is_folded (int *limit);

int sca_event_action_break (void *dummy);
// -------------------------------------------
// Event Constractors - End
// -------------------------------------------

void sca_shiftParticles();
int sca_shiftParticles_Event (void *data);

int sca_recordState_Event (void *Precord);
int sca_recordState(char *filepath);
int sca_loadState(char *filepath);

int run_loop (simtime timelimit);



struct simstate {
	/*!< recording simulation state
	 */
	simtime currentTimeStep;
	int seed;
	vec3 *particles;
	vec3 *velocities;
	vec3 *forces;
};
typedef struct simstate sca_simstate;

sca_simstate sca_simstate_init ();
void sca_simstate_free (struct simstate *C);
void sca_simstate_record (struct simstate *C);
void sca_simstate_load (struct simstate *C);
void sca_simstate_setnull (struct simstate *C);
void sca_simstate_writeTofile (struct simstate *C, char *filepath);
void sca_simstate_loadfromFile (sca_simstate *stateOut, char *filepath);

struct cache_state {
	/* for multiple states */
	int size;
	int current_index;
	int counter;
	struct simstate *data;
};

typedef struct cache_state cacheStates;
cacheStates* sca_cache_init (int size);
void sca_cache_setnull_cacheStates (cacheStates *C);
int sca_cache_record_event (simtime time_interval, int cache_state_size, cacheStates **C);
void sca_cache_load_lastState (cacheStates *C);
void sca_cache_load_oldestState (cacheStates *C);
void sca_cache_recordCurrentState (cacheStates *C);
sca_simstate* sca_cache_get_element (cacheStates *C, int elementIndex);
void sca_cache_recordState (cacheStates *C, sca_simstate S);
vec3 sum_vec3List_vectorel (vec3 *vec3list, int N);
float sum_vec3List_scalar (vec3 *vec3list, int N);

#endif /* sca_h */



