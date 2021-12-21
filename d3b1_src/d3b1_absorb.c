/*
 * d3b1_absorb.c
 *
 *  Created on: Dec 21, 2021
 *      Author: ohta
 */

#include "bem3_emf_b1.h"

int absorb_P(double *P,int type,DOMD *md)
{
  void bil_absorb_4p(double *P,int s,DOMD *md);
  void bil_absorb_9p(double *P,int s,DOMD *md);
  void lit_absorb_4p(double *P,int s,DOMD *md);
  void lit_absorb_7p(double *P,int s,DOMD *md);
  
  double tp,p;
  int t,td;

  p=0.0;
  
  #pragma omp parallel for schedule(dynamic) reduction(+:p) private(td,tp)
  for(t=1;t<=md->bd.sb[0].Ne;t++){
    td=md->bd.sb[0].sid[t];
    if( ELT4==check_element_type(td,&(md->bd)) ){
      if(type==0) bil_absorb_4p(&tp,t,md);
      else bil_absorb_9p(&tp,t,md);
    }
    else {
      if(type==0) lit_absorb_4p(&tp,t,md);
      else lit_absorb_7p(&tp,t,md);
    }
    p+=tp;
  }
  *P=-p;
  
  if(*P<0.0){
    *P*=-1.0;
    return -1;
  }
  else return 0;
}

////////////////////////////////////////////////////////////////////////
void bil_absorb_4p(double *P,int s,DOMD *md)
{
  double complex te[3],th[3],e[3],h[3];
  double vp[3],r[3],w[3];
  int i,j,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  *P=0.0;

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      e[j]=md->bd.sb[0].E[s][i][j];
      h[j]=md->bd.sb[0].H[s][i][j];
      r[j]=md->bd.ren[asd][i][j];
      w[j]=md->bd.wen[asd][i][j];
    }
    calc_mfb_EH(te,th,r,&(md->mw)); // incident field
    for(j=0;j<3;j++){
      e[j]+=te[j];
      h[j]+=th[j];
    }
    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*md->bd.wt_44[i];
  }
  *P*=-0.5;
}

void bil_absorb_9p(double *P,int s,DOMD *md)
{
  double complex e[3],h[3];
  double vp[3],r[3],w[3],cr[3][4],cw[3][3];
  int i,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  bil_copy_elem_const_rw(cr,cw,asd,&(md->bd));

  *P=0.0;
  
  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,md->bd.zt_49[i],md->bd.et_49[i],cr,cw);
    EH_t_bd(e,h,0,s,md->bd.zt_49[i],md->bd.et_49[i],2,md); // total field on boundary

    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*md->bd.wt_49[i];
  }
  *P*=-0.5;
}

void lit_absorb_4p(double *P,int s,DOMD *md)
{
  double complex te[3],th[3],e[3],h[3];
  double vp[3],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(md->bd));

  *P=0.0;

  for(i=0;i<4;i++){
    if(i<3){
      for(j=0;j<3;j++){
        e[j]=md->bd.sb[0].E[s][i][j];
        h[j]=md->bd.sb[0].H[s][i][j];
        r[j]=md->bd.ren[asd][i][j];
        w[j]=md->bd.wen[asd][i][j];
      }
      calc_mfb_EH(te,th,r,&(md->mw)); // incident field
      for(j=0;j<3;j++){
        e[j]+=te[j];
        h[j]+=th[j];
      }
    }
    else {
      lit_rw_zeta_eta(r,w,md->bd.zt_34[i],md->bd.et_34[i],cr,cw);
      EH_t_bd(e,h,0,s,md->bd.zt_34[i],md->bd.et_34[i],2,md); // total field on boundary
    }
    
    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*md->bd.wt_34[i];
  }
  *P*=-0.25;
}

void lit_absorb_7p(double *P,int s,DOMD *md)
{
  double complex e[3],h[3];
  double vp[3],r[3],w[3],cr[3][4],cw[3][3];
  int i,asd;

  asd=abs(md->bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(md->bd));

  *P=0.0;

  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,md->bd.zt_37[i],md->bd.et_37[i],cr,cw);
    EH_t_bd(e,h,0,s,md->bd.zt_37[i],md->bd.et_37[i],2,md); // total field on boundary

    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*md->bd.wt_37[i];
  }
  *P*=-0.25;
}

