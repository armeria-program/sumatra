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
int sca_condition_simstep_interval (stime *intervaltime);
int sca_condition_simstep_limit (stime *limit);
// -------------------------------------------
// Event Constractors - End
// -------------------------------------------

void sca_shiftParticles();
int sca_shiftParticles_Event (void *data);

int sca_recordState_Event (void *Precord);
int sca_recordState(char *filepath);
int sca_loadState(char *filepath);

int run_loop (stime timelimit);

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
int sca_cache_record_event (stime time_interval, int cache_state_size, cacheStates **C);
void sca_cache_load_lastState (cacheStates *C);
void sca_cache_load_oldestState (cacheStates *C);
#endif /* sca_h */



