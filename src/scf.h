//
//  scf.h
//  Sumatra
//
//  Created by Gokhan SELAMET on 15/06/2017.
//  Copyright © 2017 Ferhat Ayaz. All rights reserved.
//

// Sumatra Common Forces File

#ifndef scf_h
#define scf_h

#include <stdio.h>

int scf_add_streching_force(float ck_stretching);
int scf_add_bending_force(float ck_theta);
int scf_add_torsion_force(float ck_phi1, float ck_phi3);

#endif /* scf_h */
