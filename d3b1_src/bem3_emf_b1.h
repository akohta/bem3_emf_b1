/*
 * bem3_emf_b1.h
 *
 *  Created on: Dec 14, 2018
 *      Author: ohta
 */
#ifndef BEM3_EMF_B1_H_
#define BEM3_EMF_B1_H_

#include "d3b1_const.h"
#include "d3b1_elem.h"

// ------- setup.c ----------------------------------------------------
void read_domd(int argc,char **argv,DOMD *md);    // read data file. this is for the solver
void print_domd(DOMD *md);                        // print readed data
void print_domd_mksa(DOMD *md);                   // print readed data in MKSA system of units
void initialize_domd(DOMD *md);                   // initialize the data. this is for the solver
void finalize_domd(DOMD *md);                     // free allocated memory.
int domain_id(double *rt,DOMD *md);               // return domain id of point rt, return the main domain id if rt is on boundary
void dat_read_domd(char *filename,DOMD *md);      // read datafile outputed by dat_write_domd()
void dat_write_domd(char *filename,DOMD *md);     // write all data in binary file with specified filename
void output_node_particles(char *fname,DOMD *md); // outputs the nodes as point cloud data ( .particles file ) 


// ------- solve_bieq.c -----------------------------------------------
void solve_bieq(DOMD *md);                     // solve boundary integral equations. this is for the solver


// ------- d3b1_field.c -----------------------------------------------
// modified electromagnetic potential by using boundary integral equations
int mEMP_s(double complex *U,double *rt,int type,DOMD *md); // scattered or internal field
int mEMP_t(double complex *U,double *rt,int type,DOMD *md); // total field
int mEMP_i(double complex *U,double *rt,int type,DOMD *md); // incident field
// outputs
// U[0]=U_x, U[1]=U_y, U[2]=U_z ( vector potential ), U[4]=phi ( scalar potential ), return domain id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of DOMD object.
// type=0:4-point GL, type=1:9-pont or 7-point GL, type=2:GLN-point GL, type=3:GHN-point GL, type=4:DE integration.
// the higher the type numbers are, the lower error is, the slower calculation speed is. it is usually setted to 0 or 1. 
// the type must be set higher number when analyse near the boundary ( in the order of a element ).

// electromagnetic field by using boundary integral equations
int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md); // scattered or internal field
int EH_t(double complex *E,double complex *H,double *rt,int type,DOMD *md); // total field
int EH_i(double complex *E,double complex *H,double *rt,int type,DOMD *md); // incident field
// outputs
// E[0]=E_x, E[1]=E_y, E[2]=E_z, H[0]=H_x, H[1]=H_y, H[2]=H_z, return domain id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of DOMD object.
// type is the same as mEMP_*() functions.

// electromagnetic field by using derivative boundary integral equations 
int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DOMD *md); // scattered or internal field
int EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,DOMD *md); // total field
int EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,DOMD *md); // incident field
// outputs and inputs are the same as EH_*() functions.
// these fields are calculated by modified electromagnetic potential using derivative boundary integral equations.
// these are slower than EH_*() functions. the error is smaller than EH_*() functions in far-field.

// electromagnetic field on the boundary
void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md);
void EH_t_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md);
void EH_i_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md);
// outputs
// E[0]=E_x, E[1]=E_y, E[2]=E_z, H[0]=H_x, H[1]=H_y, H[2]=H_z.
// inputs
// did:domain id, t:element id defined in each domain, zeta_t,eta_t: point parameter on element surface (main domain side coordinate), type, pointer of DOMD object.

// tangential directional derivatives of mEMPs at node
int dmEMP_dt_node(double complex *U,double complex *dUdz,double complex *dUde,double *vtz,double *vte,int did,int eid,int nid,int type,DOMD *md);
// outputs
// U[0]=U_x, U[1]=U_y, U[2]=U_z, dUdz[0]=dU_x/dzeta, dUdz[1]=dU_y/dzeta, dUdz[2]=dU_z/dzeta, dUde[0]=dU_x/deta, dUde[1]=dU_y/deta, dUde[2]=dU_z/deta,
// unit vector of zeta directon=(vtz[0],vtz[1],vtz[2]), unit vector of eta direction=(vte[0],vte[1],vte[2]),
// return 0:normal termination, return -1:abnormal termination.
// inputs 
// did:domain id, eid:element id defined in each domain, nid:node id, type, pointer of DOMD object


// ------- d3b1_force.c -----------------------------------------------
// calculation of net radiation force and torque
void force_FN(double *F,double *N,double *rc,int type,DOMD *md);
// outputs
// F[0]=F_x, F[1]=F_y, F[2]=F_z ( radiation force ), N[0]=N_x, N[1]=N_y, N[2]=N_z ( radiation torque ).
// inputs
// rc[0]=rc_x, rc[1]=rc_y, rc[2]=rc_z ( center of rotation ), type, pointer of DOMD object.
// type=0:4-point GL, type!=0:9-point or 7-point GL.
void force_FN2(double *F,double *N,double *rc,int did,int type,DOMD *md);
// integrate only the boundary shared with the open region and the domain "did".
// did=0 : the same as force_FN().
// others are the same as force_FN().


// ------- d3b1_absorb.c -----------------------------------------------
// calculation of net absorbed energy
int absorb_P(double *P,int type,DOMD *md);
// outputs
// P : absorbed energy 
// inputs
// type=0:4-point GL, type!=0:9-point or 7-point GL.
// return 
//  0 : normal termination
// -1 : loss of significance (catastrophic cancellation) occurred during the surface integral of Poynting vector.
//      The absorbed energy is even smaller than the returned value.
int absorb_P2(double *P,int did,int type,DOMD *md);
// integrate only the boundary shared with the open region and the domain "did".
// did=0 : the same as absorb_P().
// others are the same as absorb_P().


#endif /* BEM3_EMF_B1_H_ */
