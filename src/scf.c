//
//  scf.c
//  Sumatra
//
//  Created by Gokhan SELAMET on 15/06/2017.
//  Copyright © 2017 Ferhat Ayaz. All rights reserved.
//

/*
 Sumatra Common Forces File
 */

#include "scf.h"
#include "smtr.h"

// -------------------------------------------
// Stretching Force
// -------------------------------------------
struct stretching_force_context {
    float ck_stretching;
    int size;
    float *nativeDistance;
};

void update_stretching_forces(void *userData, vec3 *results) {
    struct stretching_force_context *sfg;
    vec3 force, normv, sub;
    int i;
    
    float dEnergy;
    float dist;
    
    sfg = userData;
    vec3 *p = smtr_ctx->particles;
	
	for (i = 0; i < sfg->size; i++) {
        dist = svec3_dist(i, i + 1);
        dEnergy = -2 * sfg->ck_stretching * (dist - sfg->nativeDistance[i]);
        sub = vec3_sub(p[i], p[i + 1]);
        normv = vec3_mul(sub , 1/dist); // TODO: try div vs mul for speed. vec3_div(sub, dist);
        force = vec3_mul(normv, dEnergy);
        
        results[i] = vec3_add(results[i], force);
        results[i + 1] = vec3_sub(results[i + 1], force);
    }
}

float calculate_stretching_energy(void *userData) {
    float rij,nativeDist,sumEnergy=0.0f, energy = 0.0f;
    struct stretching_force_context *sfc=userData;
    int i;
    float ck_s = sfc->ck_stretching;
    
    for (i=0; i<sfc->size; i++) {
        nativeDist = sfc->nativeDistance[i];
        rij = svec3_dist(i, i+1);
        energy = ck_s * powf( (rij-nativeDist), 2.0f);
        sumEnergy += energy;
    }
    return sumEnergy;
}

int scf_add_streching_force(float ck_stretching) {
    int forceIndex;
    struct stretching_force_context *sfg;
    int i;
	
	sfg = calloc(1, sizeof(struct stretching_force_context));
    sfg->size = smtr_ctx->particleCount - 1;
    sfg->ck_stretching = ck_stretching;
    sfg->nativeDistance = calloc(sfg->size, sizeof(float));
	
	for (i = 0; i < sfg->size; i++) {
        sfg->nativeDistance[i] = vec3_dist(smtr_ctx->particles[i],
                                           smtr_ctx->particles[i+1]);
    }
	
    forceIndex = smtr_add_force(update_stretching_forces, sfg);
    smtr_add_energy(calculate_stretching_energy, sfg);
    return forceIndex;
}
// -------------------------------------------
// Stretching Force END
// -------------------------------------------

// -------------------------------------------
// Bending Force
// -------------------------------------------
struct bending_force_context {
    float *nativeTetha;
    float *currentTetha;
    int size;
    float ck_theta;
};
/* 
 Angle Calculation Defined as a macro:
 To have a consistent calculation between native and current angles.
 Force calculation function needs also u , v , and cosTheta values for directional derivative
 so that this values has calculated from macro as below. 
 This is not a pretty solution, but it is the one that i have found,
 solves both performance and code handling together so far.
*/
#define EPSTHETA 1.0 - 1.0e-6
#define bending_angle(i,j,k,u,v,cosTheta) { \
												u = svec3_unit(i, j); \
												v = svec3_unit(k, j); \
												cosTheta = vec3_dot(u, v); \
												if (fabs(cosTheta) > epsTheta) \
												{ \
													cosTheta = cosTheta >= 0 ? fabs(epsTheta) : -fabs(epsTheta); \
												} \
}

float calculate_bending_energy (void *userData) {
	/* Bending Energy
	 	ck_theta * ( theta - theta_0 )^2
	*/
	float sumEnergy = 0.0f, energy, delta;
    struct bending_force_context *bfc = userData;
    int i;
    
    for (i=0; i<bfc->size; i++) {
		delta = bfc->currentTetha[i] - bfc->nativeTetha[i];
		energy = bfc->ck_theta * delta * delta;
        sumEnergy += energy;
    }
    return sumEnergy;
}

static
void update_bending_force(int i, int j, int k, vec3 *results, void *userData) {
    static const float epsTheta = 1.0 - 1.0e-6;
    struct bending_force_context *bf = userData;
    float r_ij = svec3_dist(i, j);
    float r_kj = svec3_dist(j, k);
	vec3 u, v;
    vec3 dr_i, dr_k, tmp;
    float theta, sinTheta, cosTheta, theta0;
    float dEnergy;
	
	bending_angle(i, j, k, u, v, cosTheta);
	theta = acosf(cosTheta);
	bf->currentTetha[i] = theta;
	
	sinTheta = sinf(theta);
    theta0 = bf->nativeTetha[i];
    dEnergy = -2*bf->ck_theta*(theta - theta0);
    
    dr_i = vec3_div(vec3_sub(vec3_mul(u, cosTheta), v), r_ij*sinTheta);
    dr_k = vec3_div(vec3_sub(vec3_mul(v, cosTheta), u), r_kj*sinTheta);
    tmp = vec3_mul(vec3_add(dr_i, dr_k), -1);
    //printf ("%2d dE: %8.5f drI: %s \tdrK: %s \n", i, dEnergy, vec3_repr(dr_i), vec3_repr(dr_k));
    results[i] = vec3_add(results[i], vec3_mul(dr_i, dEnergy));
    results[j] = vec3_add(results[j], vec3_mul(tmp, dEnergy));
    results[k] = vec3_add(results[k], vec3_mul(dr_k, dEnergy));
    
    if (vec3_isnan(results[i]) || vec3_isnan(results[j]) || vec3_isnan(results[k])) {
        printf("#Yoooo Bending\n");
        fprintf(stderr, "BND %llu %d %.3f %.3f %.3f %.3f\n", smtr_ctx->currentTimeStep, i, cosTheta, theta, theta0, dEnergy);
        // Not good man. Not good.
    }
}

void update_bending_forces(void *userData, vec3 *results)
{
    int i;
    for(i = 0; i < smtr_ctx->particleCount - 2; i++)
    	update_bending_force(i, i + 1, i + 2, results, userData);
}

int scf_add_bending_force(float ck_theta) {
    int i, size, forceIndex;
    struct bending_force_context *bf;
	vec3 _u, _v;
	float _cosTheta, epsTheta = EPSTHETA;
    
    bf = calloc(1, sizeof(struct bending_force_context));
    size = smtr_ctx->particleCount - 2;
    bf->size = size;
    bf->nativeTetha = calloc(size, sizeof(float));
    bf->currentTetha = calloc(size, sizeof(float));
	
	for (i = 0; i < size; i++) {
		bending_angle(i, i+1, i+2, _u, _v, _cosTheta);
		bf->nativeTetha[i] = acosf(_cosTheta);
		bf->currentTetha[i] = bf->nativeTetha[i];
	}
    
    bf->ck_theta = kEpsilon * ck_theta;
    forceIndex = smtr_add_force(update_bending_forces, bf);
    smtr_add_energy(calculate_bending_energy, bf);
    return forceIndex;
}
// -------------------------------------------
// Bending Force END
// -------------------------------------------

// -------------------------------------------
// Torsion Force
// -------------------------------------------

struct torsion_force_context {
	float *nativePhi, *currentPhi;
	float ck_phi1;
	float ck_phi3;
	int size;
};

float torsion_phi_angle (int i, int j, int k, int l) {
	/*
	 TODO: make sth like in the bending such as;
	 - Angle using in the Forces must have calculate in One place.
	 	to have a consistant value btw native, current, forces, energy calculations
	 - But at force calculation some values such as v_ij , v_kj, v_mj ... 
	 	needed for the direction also. Make them some how calculated only once during calculating angle
 	 */
	static const float eps = 1.0e-3;
	
	float s;
	float r_kj, r_nk,r_mj;
	
	float phi, cosphi;
	
	vec3 *p = smtr_ctx->particles;
	
	vec3 v_ij, v_kj, v_kl, v_mj, v_nk, v_il;
	float comp;
	
	
	
	
	v_ij = vec3_sub(p[i], p[j]);
	v_kj = vec3_sub(p[k], p[j]);
	v_kl = vec3_sub(p[k], p[l]);
	
	v_mj = vec3_cross(v_ij, v_kj);
	v_nk = vec3_cross(v_kj, v_kl);
	v_il = vec3_cross(v_mj, v_nk);
	
	r_kj = svec3_dist(j, k);
	r_nk = vec3_mag(v_nk);
	r_mj = vec3_mag(v_mj);
	
	if(((r_nk * r_nk) < eps) || ((r_mj * r_mj) < eps)) {
		// continue;
		fprintf(stderr, "Force execption reached at native\n");
	}
	
	
	cosphi = vec3_dot(v_nk, v_mj)/(r_nk*r_mj);
	
	if (cosphi < -1.0) cosphi = -1;
	if (cosphi >  1.0) cosphi = 1;
	
	comp = vec3_dot(v_kj, v_il);
	s = (fabs(comp) <= 0.0000001 || comp > 0) ? 1.0 : -1.0;
	//phi = (fabs(cosphi) <= 0.000001) ? s * M_PI_2 : s * acos(cosphi);
	phi = (fabs(cosphi) == 0.0e0) ? s * M_PI_2 : s * acos(cosphi);
	return phi;
}

void update_torsion_forces(void *userData, vec3 *results)
{
	struct torsion_force_context *tf = userData;
	int i, j, k, l;
	static const float eps = 1.0e-3;
	
	float s;
	float r_kj, r_nk,r_mj;
	float rijrkj,rklrkj;
	
	float phi, phi_0, diff_phi,
	cosphi, sinphi, cosdphi, sindphi,
	cos3dphi, sinp3p0, cos3p0, sin3p0;
	
	float derivative;
	float derivative1;
	float derivative2;
	
	vec3 *p = smtr_ctx->particles;
	
	vec3 v_ij, v_kj, v_kl, v_mj, v_nk, v_il;
	vec3 dr_i, dr_l, dr_j, dr_k;
	float comp;
	
	
	for(i = 0; i < tf->size; i++) {
		j = i + 1;
		k = i + 2;
		l = i + 3;
		
		v_ij = vec3_sub(p[i], p[j]);
		v_kj = vec3_sub(p[k], p[j]);
		v_kl = vec3_sub(p[k], p[l]);
		
		v_mj = vec3_cross(v_ij, v_kj);
		v_nk = vec3_cross(v_kj, v_kl);
		v_il = vec3_cross(v_mj, v_nk);
		
		r_kj = svec3_dist(j, k);
		r_nk = vec3_mag(v_nk);
		r_mj = vec3_mag(v_mj);
		
		if(((r_nk * r_nk) < eps) || ((r_mj * r_mj) < eps))
			continue;
		
		cosphi = vec3_dot(v_nk, v_mj)/(r_nk*r_mj);
		
		if (cosphi < -1.0) cosphi = -1;
		if (cosphi >  1.0) cosphi = 1;
		
		comp = vec3_dot(v_kj, v_il);
		s = (fabs(comp) <= 0.0000001 || comp > 0) ? 1.0 : -1.0;
		//phi = (fabs(cosphi) <= 0.000001) ? s * M_PI_2 : s * acos(cosphi);
		phi = (fabs(cosphi) == 0.0e0) ? s * M_PI_2 : s * acos(cosphi);
		tf->currentPhi[i] = phi;
		
		// phi = vec3_pos_dih(p[i], p[j], p[k], p[l] );
		
		phi_0 = tf->nativePhi[i];
		
		diff_phi = phi - phi_0;
		sinphi = sin(phi);
		cosdphi = cos(diff_phi);
		sindphi = sin(diff_phi);
		cos3dphi = cos(3*diff_phi);
		sinp3p0 = sin(phi-3*phi_0);
		cos3p0 = cos(3*phi_0);
		sin3p0 = sin(3*phi_0);
		
		derivative1 = tf->ck_phi1*sindphi;
		derivative2 = 3*tf->ck_phi3*(4*cosphi*cosphi*sinp3p0+3*cosphi*sin3p0-sinphi*cos3p0);
		derivative = -(derivative1+derivative2);
		
		rijrkj = vec3_dot(v_ij, v_kj)/(r_kj*r_kj);
		rklrkj = vec3_dot(v_kl, v_kj)/(r_kj*r_kj);
		
		dr_i = vec3_mul(v_mj, r_kj/(r_mj*r_mj));
		dr_l = vec3_mul(v_nk, -r_kj/(r_nk*r_nk));
		
		dr_j = vec3_sub(vec3_mul(dr_i, rijrkj - 1), vec3_mul(dr_l, rklrkj));
		dr_k = vec3_sub(vec3_mul(dr_l, rklrkj - 1), vec3_mul(dr_i, rijrkj));
		
		//printf ("%2d %9.5f dE: %8.5f drI: %s \tdrL: %s \tdrJ: %s \tdrK: %s\n", i, phi_0, derivative,
		//		vec3_repr(dr_i), vec3_repr(dr_l), vec3_repr(dr_j), vec3_repr(dr_k) );
		
		results[i] = vec3_add(results[i], vec3_mul(dr_i, derivative));
		results[j] = vec3_add(results[j], vec3_mul(dr_j, derivative));
		results[k] = vec3_add(results[k], vec3_mul(dr_k, derivative));
		results[l] = vec3_add(results[l], vec3_mul(dr_l, derivative));
		
		if (vec3_isnan(results[i]) || vec3_isnan(results[j]) ||
			vec3_isnan(results[k]) || vec3_isnan(results[l])) {
			fprintf(stderr, "TRSN %llu %d %.3f %.3f %.3f %.3f %.3f\n", smtr_ctx->currentTimeStep, i, cosphi, phi, phi_0, derivative1, derivative2);
			// Breakpoint before armageddon
		}
		
	}
}

float calculate_torsion_energy (void *userData) {
	struct torsion_force_context *tfc = userData;
	float energy, sumEnergy = 0.0f;
	float cosdphi, diff, cos3dphi;
	int i;
	
	for (i=0; i<tfc->size; i++) {
		diff = tfc->currentPhi[i] - tfc->nativePhi[i];
		cosdphi = cosf( diff );
		cos3dphi = cosf( 3*diff );
		energy = tfc->ck_phi1 * (1 - cosdphi) + tfc->ck_phi3 * (1 - cos3dphi);
		sumEnergy += energy;
	}
	return sumEnergy;
}


int scf_add_torsion_force(float ck_phi1, float ck_phi3) {
	int i, forceIndex;
	//vec3 *p = smtr_ctx->particles;
	struct torsion_force_context *tf;
	
	tf = calloc(1, sizeof(struct torsion_force_context));
	tf->size = smtr_ctx->particleCount - 3;
	tf->nativePhi = calloc(tf->size, sizeof(float));
	tf->currentPhi = calloc(tf->size, sizeof(float));
	
	for (i = 0; i < tf->size; i++) {
		tf->nativePhi[i] = torsion_phi_angle(i, i + 1, i + 2, i + 3);
		tf->currentPhi[i] = tf->nativePhi[i];
	}
	tf->ck_phi1 = kEpsilon * ck_phi1;
	tf->ck_phi3 = kEpsilon * ck_phi3;
	
	forceIndex = smtr_add_force(update_torsion_forces, tf);
	smtr_add_energy(calculate_torsion_energy, tf);
	return forceIndex;
}

// -------------------------------------------
// Torsion Force END
// -------------------------------------------
