//
//  sca.c
//  Sumatra
//
//  Created by Gokhan Selamet on 19/06/17.
//  Copyright © 2017 Ferhat Ayaz. All rights reserved.
//

/*
 Sumatra Common Actions File
 */

#include "sca.h"
#include <string.h>
#include <stdlib.h>

vec3 sum_vec3List_vectorel (vec3 *vec3list, int N);
float sum_vec3List_scalar (vec3 *vec3list, int N);
int sca_recordState(char *filepath);

int condition_simstep_interval (void *data);



// -------------------------------------------
// Cache State Action
// -------------------------------------------

int record_cache_state_event();

struct simstate sca_simstate_init () {
	struct simstate O;
	O.particles = (vec3*) calloc(smtr_ctx->particleCount, sizeof(vec3));
	O.velocities = (vec3*) calloc(smtr_ctx->particleCount, sizeof(vec3));
	O.forces = (vec3*) calloc(smtr_ctx->particleCount, sizeof(vec3));
	return O;
}

cacheStates* sca_cache_init (int size) {
	int i=0;
	cacheStates *O;
	O = (cacheStates*) calloc (1, sizeof(struct cache_state));
	O->data = (struct simstate*) calloc(size, sizeof(struct simstate));
	
	O->size = size;
	O->current_index = -1;
	O->counter = 0;
	
	for (i=0; i<size; i++) {
		O->data[i] = sca_simstate_init();
	}
	return O;
}

void sca_simstate_free (struct simstate *C) {
	free(C->particles);
	free(C->velocities);
	free(C->forces);
}

void sca_simstate_record (struct simstate *C) {
	/*!
	 Record current state from smtr_ctx to given `C`
	 */
	C->currentTimeStep = smtr_ctx->currentTimeStep;
	C->seed = smtr_ctx->SEED;
	int N = smtr_ctx->particleCount;
	memcpy(C->particles, smtr_ctx->particles, N * sizeof(vec3));
	memcpy(C->velocities, smtr_ctx->velocities, N * sizeof(vec3));
	memcpy(C->forces, smtr_ctx->forces, N * sizeof(vec3));
}

void sca_simstate_load (struct simstate *C) {
	smtr_ctx->currentTimeStep = C->currentTimeStep;
	smtr_ctx->SEED = C->seed;
	srand(smtr_ctx->SEED);
	int N = smtr_ctx->particleCount;
	memcpy(smtr_ctx->particles, C->particles, N * sizeof(vec3));
	memcpy(smtr_ctx->velocities, C->velocities, N * sizeof(vec3));
	memcpy(smtr_ctx->forces, C->forces, N * sizeof(vec3));
	smtr_update_distances();
}

void sca_simstate_writeTofile (struct simstate *C, char *filepath) {
	/* Record sim state to to file */
	//if (smtr_ctx->currentTimeStep != 500) return 0;
	FILE *out;
	out = fopen(filepath, "wb");
	int N;
	// TODO: move write file here using cache_states actions
	
	N = smtr_ctx->particleCount;
	fwrite(&C->currentTimeStep, sizeof(simtime), 1, out);
	fwrite(&C->seed, sizeof(int), 1, out );
	fwrite(C->particles, sizeof(vec3), N, out);
	fwrite(C->velocities, sizeof(vec3), N, out);
	fwrite(C->forces, sizeof(vec3), N, out);
	fclose(out);
}

void free_cacheStates (cacheStates *C) {
	int i;
	struct simstate *data;
	for (i=0; i<C->size; i++) {
		data = &C->data[i];
		sca_simstate_free(data);
	}
	free(C->data);
	free(C);
}

void sca_simstate_setnull (struct simstate *C) {
	C->currentTimeStep = -1;
	C->seed = 0;
	memset(C->particles, 0, smtr_ctx->particleCount * sizeof(vec3));
	memset(C->velocities, 0, smtr_ctx->particleCount * sizeof(vec3));
	memset(C->forces, 0, smtr_ctx->particleCount * sizeof(vec3));
}

void sca_cache_recordCurrentState (cacheStates *C) {
	/*!
	 Record current state to cacheStates
	 If caches states is not full uses new block, else removes oldest cacheState data and records its place
	 */
	C->current_index = ( C->current_index + 1 ) % C->size;
	C->counter++;
	sca_simstate_record(&C->data[C->current_index]);
}

void sca_cache_setnull_cacheStates (cacheStates *C) {
	int i=0;
	C->current_index=-1;
	C->counter = 0;
	for (i=0; i<C->size; i++) {
		sca_simstate_setnull(&C->data[i]);
	}
}

sca_simstate* sca_cache_get_element (cacheStates *C, int elementIndex) {
	if (elementIndex >= C->size ) { return NULL; }
	if (C->counter < C->size) { return &C->data[elementIndex]; }
	else { return &C->data[ (C->current_index+1 % C->size) + elementIndex ]; }
}

void sca_cache_load_lastState (cacheStates *C) {
	int i = (C->current_index + C->size - 1 ) % C->size ;
	sca_simstate_load(&C->data[i]);
}

void sca_cache_load_oldestState (cacheStates *C) {
	int i;
	if (C->counter >= C->size)
		i = (C->current_index + 1) % C->size;
	else
		i = 0;
	sca_simstate_load(&C->data[i]);
}



int sca_cache_record_event (simtime time_interval, int cache_state_size, cacheStates **C) {
	struct event_data *E;
	simtime *interval;
	E = malloc(sizeof(struct event_data));
	interval = init_stime();
	*interval = time_interval;
	
	*C = sca_cache_init(cache_state_size);
	
	E->condition_function = (int(*)(void*)) sca_condition_simstep_interval;
	E->condition_data=interval;
	
	E->action_function= (int(*)(void*)) sca_cache_recordCurrentState;
	E->action_data=*C;
	
	smtr_subscribe_event( (SmtrCallbackFunc) sca_event, E);
	return 0;
}

// -------------------------------------------
// Cache State Action END
// -------------------------------------------

// -------------------------------------------
// Record and Load State Action
// -------------------------------------------
int sca_recordState_Event (void *Precord) {
	struct simstep_interval *P = Precord;
	if (smtr_ctx->currentTimeStep % P->interval != 0) return 0;
	char name[255];
	char *base = P->DATA;
	sprintf(name, "%s_%llu.bin", base, smtr_ctx->currentTimeStep);
	sca_recordState(name);
	return 0;
}

int sca_recordState(char *filepath) {
	/* Record sim state to to file */
	//if (smtr_ctx->currentTimeStep != 500) return 0;
	FILE *out;
	out = fopen(filepath, "wb");
	
	struct simstate C;
	C = sca_simstate_init();
	sca_simstate_record(&C);
	
	// TODO: move write file here using cache_states actions
	
	sca_simstate_free(&C);
	
	//_write_cache_state_tofile(&C, filepath);
	
	fwrite(&smtr_ctx->currentTimeStep, sizeof(simtime), 1, out);
	fwrite(&smtr_ctx->SEED, sizeof(int), 1, out );
	fwrite(smtr_ctx->particles, sizeof(vec3), smtr_ctx->particleCount, out);
	fwrite(smtr_ctx->velocities, sizeof(vec3), smtr_ctx->particleCount, out);
	fwrite(smtr_ctx->forces, sizeof(vec3), smtr_ctx->particleCount, out);
	
	/*
	 int i;
	 printf("F %d %llu\n", smtr_ctx->SEED, smtr_ctx->currentTimeStep);
	 for (i=0; i<smtr_ctx->particleCount; i++) {
	 printf ("%d %f;%f;%f %f;%f;%f %f;%f;%f\n", i, smtr_ctx->particles[i].x,smtr_ctx->particles[i].y,smtr_ctx->particles[i].z,
	 smtr_ctx->velocities[i].x,smtr_ctx->particles[i].y,smtr_ctx->particles[i].z,
	 smtr_ctx->forces[i].x,smtr_ctx->particles[i].y,smtr_ctx->particles[i].z);
	 }*/
	fclose(out);
	return 0;
}

int sca_loadState(char *filepath) {
	/* Record sim state to continue later */
	FILE *out;
	out = fopen(filepath, "rb");
	fread(&smtr_ctx->currentTimeStep, sizeof(simtime), 1, out);
	fread(&smtr_ctx->SEED, sizeof(int), 1, out );
	srand(smtr_ctx->SEED);
	fread(smtr_ctx->particles, sizeof(vec3), smtr_ctx->particleCount, out);
	//updateDistance here when moved in sumatra
	fread(smtr_ctx->velocities, sizeof(vec3), smtr_ctx->particleCount, out);
	fread(smtr_ctx->forces, sizeof(vec3), smtr_ctx->particleCount, out);
	
	/*
	 int i;
	 printf("F %d %llu\n", smtr_ctx->SEED, smtr_ctx->currentTimeStep);
	 for (i=0; i<smtr_ctx->particleCount; i++) {
	 printf ("%d %f;%f;%f %f;%f;%f %f;%f;%f\n", i, smtr_ctx->particles[i].x,smtr_ctx->particles[i].y,smtr_ctx->particles[i].z,
	 smtr_ctx->velocities[i].x,smtr_ctx->particles[i].y,smtr_ctx->particles[i].z,
	 smtr_ctx->forces[i].x,smtr_ctx->particles[i].y,smtr_ctx->particles[i].z);
	 }*/
	fclose(out);
	return 0;
}
// -------------------------------------------
// Record and Load State Action - END
// -------------------------------------------

eventData* sca_init_eventData () {
	return malloc(1 * sizeof(eventData));
}

struct _time_limit {
	simtime timelimit;
};

int sca_event (eventData *E) {
	//struct event_data *E = D;
	if ( E->condition_function (E->condition_data) )
		E->action_function (E->action_data);
		return 0;
	}

int sca_condition_simstep_interval (simtime *intervaltime) {
	/*!<
	 For every `intervaltime` steps, return True
	 Usefull such as printing trajectory files each 400 timesteps etc.
	 */
	if (smtr_ctx->currentTimeStep % *intervaltime == 0)
		return 1;
	return 0;
}

int sca_condition_simstep_limit (simtime *limit) {
	/*!
	 After given `limit` time step expired return True
	 Usefull such as ending simulation in givin time step.
	 */
	if (smtr_ctx->currentTimeStep > *limit)
		return 1;
	return 0;
}

int sca_condition_is_folded (int *limit) {
	/*!
	*/
	return 0;
}

int sca_event_action_break (void *dummy) {
	/*!<
	 This function casted original `smtr_break function into from
	 that sca_event can use as an action
	 usage: pass null `dummy` is not using
	 */
	smtr_break();
	return 0;
}

int run_loop (simtime timelimit) {
	simtime *T;
	T = init_stime();
	*T = timelimit;
	
	struct event_data *E;
	E = malloc(1 * sizeof(struct event_data));
	
	E->condition_function = sca_condition_simstep_limit;
	E->condition_data = T;
	
	E->action_function = sca_event_action_break;
	E->action_data = NULL;
	
	smtr_subscribe_event(sca_event, E);
	smtr_run();
	
	// TODO: free E and T
	return 0;
}

// -------------------------------------------
// Shift Particles Action
// -------------------------------------------
int sca_shiftParticles_Event (void *data) {
	struct simstep_interval *P = data;
	if (smtr_ctx->currentTimeStep % P->interval != 0 ) return 0;
	sca_shiftParticles();
	return 0;
}

void sca_shiftParticles() {
	/* 
	 Shift all particles into center point
	 */
	vec3 center;
	int i, N = smtr_ctx->particleCount;
	vec3 *P = smtr_ctx->particles;
	center = sum_vec3List_vectorel(P, N);
	center = vec3_div(center, N);
	for (i=0; i<N; i++) {
		P[i] = vec3_sub(P[i], center);
	}
}
// -------------------------------------------
// Record and Load State - END
// -------------------------------------------

/* 
 Events in the below was used in during testing. 
 Some of them should be usefull for later.
 */
int printForce (void *data) {
	int i, x;
	vec3 *F, *f;
	
	for (i=0; i<smtr_ctx->forceCount; i++){
		F = smtr_ctx->forceClasses[i].forces;
		for (x=0; x<smtr_ctx->particleCount; x++) {
			f = & F[x];
			printf ("%d - %d %f %f %f\n", i, x, f->x, f->y, f->z );
		}
	}
	return 0;
}

int printVelocity (void *data) {
	int i = 0 , n=smtr_ctx->particleCount;
	vec3 *V;
	float sum = 0.0f;
	
	for (i=0; i<n; i++) {
		V = & smtr_ctx->velocities[i];
		//printf("V - %d - %f %f %f\n", i, V->x, V->y, V->z);
		sum += fabsf(V->x) + fabsf(V->y) + fabsf(V->z);
	}
	printf ("Velocity = %f\n", sum);
	
	return 0;
}

int printDistanceDisplacement (void *data) {
	struct simstep_interval *C = data;
	int i,N = smtr_ctx->particleCount;
	float sum,dist;
	vec3 *coor = smtr_ctx->particles;
	static vec3 *prevcoor = NULL;
	if (!prevcoor) prevcoor = malloc( N * sizeof(vec3));
	
	if (smtr_ctx->currentTimeStep % (C->interval) == C->interval-1 )
		memcpy(prevcoor, coor, N * sizeof(vec3));
	else if (smtr_ctx->currentTimeStep % (C->interval) == 0 ) {
		sum = 0.0;
		fprintf(C->DATA,"%llu ", smtr_ctx->currentTimeStep);
		for (i=0; i<N; i++) {
			dist = vec3_dist(prevcoor[i], coor[i]);
			fprintf (C->DATA,"%.2f " , dist);
			sum += dist;
		}
		fprintf(C->DATA,"%.3f\n", sum / N);
	}
	return 0;
}

int forceSumTest (void *data) {
	int i, x;
	vec3 *F, *f;
	vec3 sum={0.0, 0.0, 0.0};
	
	for (i=1; i<smtr_ctx->forceCount; i++){
		F = smtr_ctx->forceClasses[i].forces;
		sum = (vec3){0.0, 0.0, 0.0};
		for (x=0; x<smtr_ctx->particleCount; x++) {
			f = & F[x];
			sum = vec3_add(sum, *f);
		}
		printf ("Force   Sum x:%f y:%f z:%f\n", sum.x, sum.y, sum.z);
	}
	sum = (vec3){0.0, 0.0, 0.0};
	for (i=0; i<smtr_ctx->particleCount; i++) {
		f = & smtr_ctx->forces[i];
		sum = vec3_add(sum, *f);
	}
	printf ("Force T Sum x:%f y:%f z:%f\n", sum.x, sum.y, sum.z);
	
	return 0;
}

int printDistances(void *data) {
	int i, j, N;
	N = smtr_ctx->particleCount;
	
	printf ("%llu %6.3f \n", smtr_ctx->currentTimeStep,svec3_dist(9, 14));
	/*for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
	 printf("%2.1f ", svec3_dist(i, j) );
		}
		printf("\n");
	 }
	 */
	//printf("------------------------------------------");
	return 0;
}

int printDistanceDisplacement_Event (void *data) {
	struct simstep_interval *C = data;
	int i,N = smtr_ctx->particleCount;
	float sum,dist;
	vec3 *coor = smtr_ctx->particles;
	static vec3 *prevcoor = NULL;
	if (!prevcoor) prevcoor = malloc( N * sizeof(vec3));
	
	if (smtr_ctx->currentTimeStep % (C->interval) == C->interval-1 )
		memcpy(prevcoor, coor, N * sizeof(vec3));
	else if (smtr_ctx->currentTimeStep % (C->interval) == 0 ) {
		sum = 0.0;
		fprintf(C->DATA,"%llu ", smtr_ctx->currentTimeStep);
		for (i=0; i<N; i++) {
			dist = vec3_dist(prevcoor[i], coor[i]);
			fprintf (C->DATA,"%.2f " , dist);
			sum += dist;
		}
		fprintf(C->DATA,"%.3f\n", sum / N);
	}
	return 0;
}

// -------------------------------------------
// Helpfull Common Functions
// -------------------------------------------
vec3 sum_vec3List_vectorel (vec3 *vec3list, int N) {
	/*
	 Returns vectorel sum of given vec3list
	 */
	int i;
	vec3 sum = {0,0,0};
	for (i=0; i<N; i++) {
		sum = vec3_add(sum, vec3list[i]);
	}
	return sum;
}

float sum_vec3List_scalar (vec3 *vec3list, int N) {
	/*
	 Returns scalar sum of given vec3list
	 */
	int i;
	float sum = 0;
	for (i=0; i<N; i++) {
		sum += vec3_mag(vec3list[i]);
	}
	return sum;
}
/*
int run_kinetic (int desiredSampleSize, cacheStates startState) {
    int i;
	eventData *ed;
	ed = sca_init_eventData();
	
	ed->action_function = sca_event_action_break;
	ed->action_data = NULL;
	
	ed->condition_data = NULL;
	ed->condition_function = NULL;
	
    for (i=0; i<desiredSampleSize; i++) {
        smtr_run();
		// Action for init particles
    }
    return 0;
}
*/
