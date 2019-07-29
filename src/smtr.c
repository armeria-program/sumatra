//
//  sim.c
//  Sumatra
//
//  Created by Ferhat Ayaz on 02/06/16.
//  Copyright © 2016 Ferhat Ayaz. All rights reserved.
//

#include "smtr.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> // temp

SmtrContext _smtr;
SmtrContext *smtr_ctx = &_smtr;


// gamma: friction coefficient / mass
static float kGamma;

// delta: integration time step
float kIntegrationStep;

// Time unit
float kTimeUnit;

// velocity based verlet algorithm coeeficients
// used in equation of motion
#define c0  (kIntegrationStep * kGamma * 0.5)
float kVelocityScaleFactor;

// eta: random number generator
static float eta(void);

// Local functions
void calculate_forces(void);
static void update_coordinates(void);
//static void update_distances();
static void update_random_forces(void *userData, vec3 *results);
static int notify_subscribers(void);
static void add_random_force(float temperature);

// Extension to vec3
static vec3 vec3_random(float scalar);

simtime* init_stime() {
	return malloc(1 * sizeof(simtime));
}

int smtr_init(vec3 *partcls, float *mass, int size, float temperature)
{
  int i;
  
  kTimeUnit = sqrt(kMass/kEpsilon)*kAlpha;
  kGamma = 0.05/kTimeUnit;
  kIntegrationStep = 0.005*kTimeUnit;
  kVelocityScaleFactor = (1 - c0) * (1 - c0 + c0 * c0);

  _smtr.particleCount = size;
  _smtr.forceCount = 0;
  
  _smtr.distances = calloc(size * size, sizeof(float));
  _smtr.neighbours = calloc(size, sizeof(int) * MAX_NEIGHBOUR_COUNT);


  _smtr.forces = calloc(size, sizeof(vec3));
  _smtr.velocities = calloc(size, sizeof(vec3));
  _smtr.particles = calloc(size, sizeof(vec3));
  _smtr.mass = calloc(size, sizeof(float));
  _smtr.vforceScalars = calloc(size, sizeof(float));

  memcpy(_smtr.particles, partcls, size * sizeof(vec3));
  memcpy(_smtr.mass, mass, size * sizeof(float));

	struct vec3 norm = {0,0,0};
	for (i = 0; i < size; i++) {
    _smtr.velocities[i] = vec3_random( sqrt(kEpsilon*temperature/mass[i]) );
		norm = vec3_add(norm, _smtr.velocities[i]);
	}
	norm = vec3_div(norm, _smtr.particleCount);
// Normalize initial random velocitys
	for (i = 0; i < size; i++) {
		_smtr.velocities[i] = vec3_sub(_smtr.velocities[i], norm);
	}

  for (i = 0; i < size; i++)
    _smtr.vforceScalars[i] = (1 - c0 + c0 * c0)*( kIntegrationStep/(2*mass[i]) );

  _smtr.currentTimeStep = 0;
  _smtr.subscriberCount = 0;
  
  smtr_update_distances();
 //   printf ("ktimeUnit : %f\n kGamma:%f \n kIntegrationStep: %f\n kVelocityScaleFactor: %f\n C0: %f\n vforceScalar: %f\n", kTimeUnit, kGamma, kIntegrationStep, kVelocityScaleFactor , c0, _smtr.vforceScalars[0]);
	

  if (temperature > 0) add_random_force(temperature);
	 
  return 0;
}

int smtr_add_force(SmtrUpdateForceFunc func, void *userData)
{
  assert(_smtr.forceCount < MAX_FORCE_COUNT);
  
  SmtrForceClass *fc = _smtr.forceClasses + _smtr.forceCount++;
  fc->forces = calloc(_smtr.particleCount, sizeof(vec3));
  fc->userData = userData;
  fc->update = func;
	
	return _smtr.forceCount - 1;
}

void smtr_add_energy(SmtrCalculateEnergyFunc func, void *userData )
{
	assert(_smtr.forceCount < MAX_ENERGY_COUNT);
	
	SmtrEnergyClass *fc = _smtr.energyClasses + _smtr.energyCount++;
	fc->func = func;
	fc->userData = userData;
}

float smtr_calcPotentialEnergy()
{
	int i;
	float total=0.0f;
	SmtrEnergyClass *fc;
	for (i=0; i<_smtr.energyCount; i++) {
		fc = & _smtr.energyClasses[i];
		fc->energy = fc->func(fc->userData);
		total += fc->energy;
	}
	return total;
}


static
void update_random_forces(void *userData, vec3 *results)
{
  int i;
  float *m = userData;
  
  for (i = 0; i < smtr_ctx->particleCount; i++)
    results[i] = vec3_random(m[i]);
}

static void add_random_force(float temperature)
{
  float *forceConsts;
  float multiplier;
  int i;
  
  // Calculate constants
  forceConsts = calloc(smtr_ctx->particleCount, sizeof(float));
  // multiplier = sqrt(20*kEpsilon*temperature)/kTimeUnit;
    multiplier = 2*kEpsilon*temperature*kGamma/kIntegrationStep;
    
    for (i = 0; i < smtr_ctx->particleCount; i++) {
    forceConsts[i] = sqrt(multiplier * smtr_ctx->mass[i]);
       // printf("Rastgele Kuvvvet Sabiti %d %f \n", i, forceConsts[i]);
    }
  // Add force
  smtr_add_force(update_random_forces, forceConsts);
}


void smtr_subscribe_event(SmtrCallbackFunc func, void *userData)
{
  assert(_smtr.subscriberCount < MAX_SUBSCRIBER_COUNT);
  
  SmtrSubscriber *sbc = _smtr.subscribers + _smtr.subscriberCount++;
  sbc->userData = userData;
  sbc->callback = func;
}


static
void update_velocities()
{
  int i;
  vec3 tmp1, tmp2;
	float sumx=0, sumy=0, sumz=0, sum=0.0;
	
  for (i = 0; i < _smtr.particleCount; i++)
  {
    tmp1 = vec3_mul(_smtr.velocities[i], kVelocityScaleFactor);
    tmp2 = vec3_mul(_smtr.forces[i], _smtr.vforceScalars[i]);
    _smtr.velocities[i] = vec3_add(tmp1, tmp2);
	}
	

  calculate_forces();
  
  for (i = 0; i < _smtr.particleCount; i++)
  {
    tmp1 = vec3_mul(_smtr.forces[i], _smtr.vforceScalars[i]);
    _smtr.velocities[i] = vec3_add(_smtr.velocities[i], tmp1);
		//sum += fabsf (_smtr.velocities[i].x) + fabsf (_smtr.velocities[i].y) + fabsf (_smtr.velocities[i].z);
	}
	//printf ("Velocity Sum: %f\n", sum );
}


#if 0
void calculate_forces()
{
  int i;
  SmtrForceClass *fc;
  
  memset(_smtr.forces, 0, sizeof(vec3) * _smtr.particleCount);
  
  for (i = 0; i < _smtr.forceCount; i++)
  {
    fc = _smtr.forceClasses + i;
    fc->update(fc->userData, _smtr.forces);
  }
}
#endif

#if 1
//static
void calculate_forces()
{
  int i, j;
  SmtrForceClass *fc;
  
  memset(_smtr.forces, 0, sizeof(vec3) * _smtr.particleCount);
  
  for (i = 0; i < _smtr.forceCount; i++)
  {
    fc = _smtr.forceClasses + i;
    memset(fc->forces, 0, smtr_ctx->particleCount * sizeof(vec3));
    fc->update(fc->userData, fc->forces);
    for (j = 0; j < _smtr.particleCount; j++)
      _smtr.forces[j] = vec3_add(_smtr.forces[j], fc->forces[j]);
  }
}
#endif

static
void update_coordinates()
{
  int i;
  vec3 *v = _smtr.velocities;
  vec3 *f = _smtr.forces;
  vec3 *p = _smtr.particles;
  float *m = _smtr.mass;
  const float h = kIntegrationStep;
  const float gm = kGamma;
  for (i = 0; i < _smtr.particleCount; i++)
  {
    p[i].x += h * (v[i].x + h*(f[i].x - m[i]*gm*v[i].x) * 0.5/m[i]); //TODO: calculate one time constants
    p[i].y += h * (v[i].y + h*(f[i].y - m[i]*gm*v[i].y) * 0.5/m[i]);
    p[i].z += h * (v[i].z + h*(f[i].z - m[i]*gm*v[i].z) * 0.5/m[i]);
  }
}


void smtr_update_distances()
{
  int i, j, k;
  float dist;
  int *n;
  
  memset(_smtr.neighbours, 0, 
    sizeof(int) * MAX_NEIGHBOUR_COUNT * _smtr.particleCount);
  
  for (i = 0; i < _smtr.particleCount; i++)
  {
    n = _smtr.neighbours + i * MAX_NEIGHBOUR_COUNT;
    for (j = i + 1, k = 0; j < _smtr.particleCount; j++)
    {
      dist = vec3_dist(_smtr.particles[i], _smtr.particles[j]);
      if (dist < MAX_INTERACTING_DIST && k < MAX_NEIGHBOUR_COUNT - 2)
        n[k++] = j;
      _smtr.distances[i*_smtr.particleCount + j] = dist;
    }
  }
}

static
void update_neighbours()
{
  int i;
  int *n;
  float _tmp;
  
  for (i = 0; i < _smtr.particleCount; i++)
  {
    n = _smtr.neighbours + i * MAX_NEIGHBOUR_COUNT;
    for (; *n ; n++)
    {
      assert(i < *n);
      _smtr.distances[i*_smtr.particleCount + *n] =
        vec3_dist(_smtr.particles[i], _smtr.particles[*n]);
      _tmp = vec3_dist(_smtr.particles[i], _smtr.particles[*n]);
    }
  }
}


static
int notify_subscribers()
{
  int i, result;
  SmtrSubscriber *scb;
  
  for (i = 0; i < _smtr.subscriberCount; i++)
  {
    scb = _smtr.subscribers + i;
    result = scb->callback(scb->userData);
    if (result != 0)
      return result;
  }
  return 0;
}

void smtr_run_loop(unsigned long long steps)
{
  //_smtr.currentTimeStep = 0;
  smtr_update_distances(); // TODO: remove it from here move it to loadState when it is in sumatra project
  while (_smtr.currentTimeStep < steps)
  {
		if (notify_subscribers()) break;

    update_velocities();
    update_coordinates();
    if (_smtr.currentTimeStep % DIST_CALC_INTERVAL == 0)
      smtr_update_distances();
    else
      update_neighbours();

    _smtr.currentTimeStep++;
  }
}

void smtr_run_loop2(unsigned long long steps)	{
	
	
}
//
// Allen, M. P. & Tildesley, D. J. (1987).
// Computer Simulation of Liquids, Oxford University Press, Oxford, UK. p. 347
//
// Generates normal distributed random numbers in the range (-1, 1)
// with mean 0, standard deviation 1
//
static float eta()
{
  static const float
  aa1 = 3.949846138,
  aa3 = 0.252408784,
  aa5 = 0.076542912,
  aa7 = 0.008355968,
  aa9 = 0.029899776;
  
  int i;
  float r = 0, r2, a;
  
  for (i = 0; i < 12; i++)
  {
    smtr_ctx->SEED=rand();
    a = ( smtr_ctx->SEED / (float) (RAND_MAX));
    r += a;
  }
  
  r = (r - 6) * 0.25;
  r2 = r * r;
  
  return ((((aa9 * r2 + aa7) * r2 + aa5) * r2 + aa3) * r2 + aa1) * r;
}

static vec3 vec3_random(float scalar)
{
  return (vec3) {
    scalar * eta(),
    scalar * eta(),
    scalar * eta()
  };
}

float smtr_calcKineticEnergy () {
	int i;
	float sumKineticEnergy = 0.0f;
	float mass;
	vec3 V;
	
	for (i=0; i<_smtr.particleCount; i++) {
		mass = _smtr.mass[i];
		V = _smtr.velocities[i];
		sumKineticEnergy += mass * (V.x*V.x + V.y*V.y + V.z*V.z);
	}
	return sumKineticEnergy * 0.5;
}

// -------------------------------------------
// Running Action
// -------------------------------------------
void smtr_run() {
	_smtr._running = 1;
	smtr_update_distances(); // TODO: remove it from here move it to loadState when it is in sumatra project
	while (_smtr._running) {
		if (notify_subscribers()) break;
		
		update_velocities();
		update_coordinates();
		if (_smtr.currentTimeStep % DIST_CALC_INTERVAL == 0)
			smtr_update_distances();
		else
			update_neighbours();
		
		_smtr.currentTimeStep++;
	}
}

void smtr_break () {
	_smtr._running = 0;
}

FILE* fopenSmtr (char *path, char *type) {
  FILE* F;
  extern int errno;
  F = fopen(path, type);
  if (F == NULL) {
    fprintf(stderr, "Error on file %s with %s : %s", path, type, strerror(errno));
    exit (2);
  }
  return F;
}
