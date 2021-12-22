/*
 * d3b1_force.c
 *
 *  Created on: Jan 26, 2019
 *      Author: ohta
 */

#include "bem3_emf_b1.h"

void bil_force_4p(double *F,double *N,double *rc,int s,DOMD *md);
void bil_force_9p(double *F,double *N,double *rc,int s,DOMD *md);
void lit_force_4p(double *F,double *N,double *rc,int s,DOMD *md);
void lit_force_7p(double *F,double *N,double *rc,int s,DOMD *md);


void force_FN(double *F,double *N,double *rc,int type,DOMD *md)
{
  double tf[3],tn[3],Fx,Fy,Fz,Nx,Ny,Nz;
  int t,td;

  Fx=0.0;    Fy=0.0;    Fz=0.0;
  Nx=0.0;    Ny=0.0;    Nz=0.0;

  #pragma omp parallel for schedule(dynamic) reduction(+:Fx,Fy,Fz,Nx,Ny,Nz) private(td,tf,tn)
  for(t=1;t<=md->bd.sb[0].Ne;t++){
    td=md->bd.sb[0].sid[t];
    if( ELT4==check_element_type(td,&(md->bd)) ){
      if(type==0) bil_force_4p(tf,tn,rc,t,md);
      else bil_force_9p(tf,tn,rc,t,md);
    }
    else {
      if(type==0) lit_force_4p(tf,tn,rc,t,md);
      else lit_force_7p(tf,tn,rc,t,md);
    }
    Fx+=tf[0];
    Fy+=tf[1];
    Fz+=tf[2];
    Nx+=tn[0];
    Ny+=tn[1];
    Nz+=tn[2];
  }

  F[0]=Fx;
  F[1]=Fy;
  F[2]=Fz;
  N[0]=Nx;
  N[1]=Ny;
  N[2]=Nz;
}

void force_FN2(double *F,double *N,double *rc,int did,int type,DOMD *md)
{
  double tf[3],tn[3],Fx,Fy,Fz,Nx,Ny,Nz;
  int t,td;

  Fx=0.0;    Fy=0.0;    Fz=0.0;
  Nx=0.0;    Ny=0.0;    Nz=0.0;

  #pragma omp parallel for schedule(dynamic) reduction(+:Fx,Fy,Fz,Nx,Ny,Nz) private(td,tf,tn)
  for(t=1;t<=md->bd.sb[0].Ne;t++){
    td=md->bd.sb[0].sid[t];
    if(md->bd.sd[abs(td)]!=did && did!=0) continue;
    if( ELT4==check_element_type(td,&(md->bd)) ){
      if(type==0) bil_force_4p(tf,tn,rc,t,md);
      else bil_force_9p(tf,tn,rc,t,md);
    }
    else {
      if(type==0) lit_force_4p(tf,tn,rc,t,md);
      else lit_force_7p(tf,tn,rc,t,md);
    }
    Fx+=tf[0];
    Fy+=tf[1];
    Fz+=tf[2];
    Nx+=tn[0];
    Ny+=tn[1];
    Nz+=tn[2];
  }

  F[0]=Fx;
  F[1]=Fy;
  F[2]=Fz;
  N[0]=Nx;
  N[1]=Ny;
  N[2]=Nz; 
}

/////////////////////////////////////////////////////////////////////////////
void bil_force_4p(double *F,double *N,double *rc,int s,DOMD *md)
{
  double complex te[3],th[3],e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3];
  int i,j,k,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  epd=creal(md->n[0])*creal(md->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      e[j]=md->bd.sb[0].E[s][i][j];
      h[j]=md->bd.sb[0].H[s][i][j];
      r[j]=md->bd.ren[asd][i][j];
      w[j]=md->bd.wen[asd][i][j];
    }
    calc_mfb_EH(te,th,r,&(md->mw));
    for(j=0;j<3;j++){
      e[j]+=te[j];
      h[j]+=th[j];
    }

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*md->bd.wt_44[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*md->bd.wt_44[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-1.0;
    N[j]*=-1.0;
  }
}

void bil_force_9p(double *F,double *N,double *rc,int s,DOMD *md)
{
  double complex e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,k,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  bil_copy_elem_const_rw(cr,cw,asd,&(md->bd));

  epd=creal(md->n[0])*creal(md->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,md->bd.zt_49[i],md->bd.et_49[i],cr,cw);
    EH_t_bd(e,h,0,s,md->bd.zt_49[i],md->bd.et_49[i],2,md);

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*md->bd.wt_49[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*md->bd.wt_49[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-1.0;
    N[j]*=-1.0;
  }
}

void lit_force_4p(double *F,double *N,double *rc,int s,DOMD *md)
{
  double complex te[3],th[3],e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,k,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(md->bd));

  epd=creal(md->n[0])*creal(md->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    if(i<3){
      for(j=0;j<3;j++){
        e[j]=md->bd.sb[0].E[s][i][j];
        h[j]=md->bd.sb[0].H[s][i][j];
        r[j]=md->bd.ren[asd][i][j];
        w[j]=md->bd.wen[asd][i][j];
      }
      calc_mfb_EH(te,th,r,&(md->mw));
      for(j=0;j<3;j++){
        e[j]+=te[j];
        h[j]+=th[j];
      }
    }
    else {
      lit_rw_zeta_eta(r,w,md->bd.zt_34[i],md->bd.et_34[i],cr,cw);
      EH_t_bd(e,h,0,s,md->bd.zt_34[i],md->bd.et_34[i],2,md);
    }

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*md->bd.wt_34[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*md->bd.wt_34[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-0.5;
    N[j]*=-0.5;
  }
}

void lit_force_7p(double *F,double *N,double *rc,int s,DOMD *md)
{
  double complex e[3],h[3];
  double epd,aE2,aH2,Tc,T[9],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,k,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(md->bd));

  epd=creal(md->n[0])*creal(md->n[0]);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,md->bd.zt_37[i],md->bd.et_37[i],cr,cw);
    EH_t_bd(e,h,0,s,md->bd.zt_37[i],md->bd.et_37[i],2,md);

    aE2=creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2]));
    aH2=creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2]));
    Tc=0.25*(epd*aE2+aH2);

    for(j=0;j<3;j++) for(k=0;k<3;k++) T[3*j+k]=0.5*(epd*creal(e[j]*conj(e[k]))+creal(h[j]*conj(h[k])));
    T[0]-=Tc;
    T[4]-=Tc;
    T[8]-=Tc;
    for(j=0;j<3;j++){
      F[j]+=( T[3*j+0]*w[0]+T[3*j+1]*w[1]+T[3*j+2]*w[2] )*md->bd.wt_37[i];
      N[j]+=( ((r[1]-rc[1])*T[3*j+2]-(r[2]-rc[2])*T[3*j+1])*w[0]
             +((r[2]-rc[2])*T[3*j+0]-(r[0]-rc[0])*T[3*j+2])*w[1]
             +((r[0]-rc[0])*T[3*j+1]-(r[1]-rc[1])*T[3*j+0])*w[2] )*md->bd.wt_37[i];
    }
  }

  for(j=0;j<3;j++){
    F[j]*=-0.5;
    N[j]*=-0.5;
  }
}
