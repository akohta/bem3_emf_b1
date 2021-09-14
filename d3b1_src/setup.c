/*
 * setup.c
 *
 *  Created on: Dec 15, 2018
 *      Author: ohta
 */

#include "bem3_emf_b1.h"

void read_domd(int argc,char **argv,DOMD *md)
{
  void filename_chk(int argc,char **argv);
  void read_medium_data(char *med_fn,DOMD *md);
  void read_mesh_data(char *msh_fn,DOMD *md);
  
  filename_chk(argc,argv);
  // multi_wave
  init_mfb(&(md->mw));
  read_data_mfb(&(md->mw));
  read_medium_data(argv[1],md);
  read_mesh_data(argv[2],md);
  
  if(argc==11){
    md->rv[0]=atof(argv[ 4]);
    md->rv[1]=atof(argv[ 5]);
    md->rv[2]=atof(argv[ 6]);
    md->th   =atof(argv[ 7]);
    md->tv[0]=atof(argv[ 8]);
    md->tv[1]=atof(argv[ 9]);
    md->tv[2]=atof(argv[10]);
  }
  else {
    md->rv[0]=1.0;
    md->rv[1]=0.0;
    md->rv[2]=0.0;
    md->th=0.0;
    md->tv[0]=0.0;
    md->tv[1]=0.0;
    md->tv[2]=0.0;
  }
}

void print_domd(DOMD *md)
{
  void print_medium_data(DOMD *md);
  void print_mesh_data(DOMD *md);
  
  printf("-- multi_fbeam data --\n");
  print_data_mfb(&(md->mw));

  print_medium_data(md);
  printf("\n");

  print_mesh_data(md);
  printf("\n");
  
  if(md->th!=0.0 || vabs_d(md->tv)!=0.0){
    printf("-- rotation and translation settings --\n");
    if(md->th!=0.0){
      printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",md->rv[0],md->rv[1],md->rv[2]);
      printf("rotation angle           [rad]: %8.7g\n",md->th);
    }
    if(vabs_d(md->tv)!=0.0){
      printf("translation vector            :(%8.7g,%8.7g,%8.7g)\n",md->tv[0],md->tv[1],md->tv[2]);
    }
  }
  printf("\n");
}

void print_domd_mksa(DOMD *md)
{
  void print_medium_data(DOMD *md);
  void print_mesh_data(DOMD *md);
  
  printf("-- multi_fbeam data --\n");
  print_data_mfb_mksa(&(md->mw));

  print_medium_data(md);
  printf("\n");

  print_mesh_data(md);
  printf("\n");
  
  if(md->th!=0.0 || vabs_d(md->tv)!=0.0){
    printf("-- rotation and translation settings --\n");
    if(md->th!=0.0){
      printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",md->rv[0],md->rv[1],md->rv[2]);
      printf("rotation angle           [rad]: %8.7g\n",md->th);
    }
    if(vabs_d(md->tv)!=0.0){
      printf("translation vector         [m]:(%8.7g,%8.7g,%8.7g)\n",OSUtoMKSA_length(md->tv[0]),OSUtoMKSA_length(md->tv[1]),OSUtoMKSA_length(md->tv[2]));
    }
  }
  printf("\n");  
}

void initialize_domd(DOMD *md)
{
  void rotation_translation_obj(double *rv,double th,double *tv,DOMD *md);
  void init_elem_const(BOUD *bd);
  void malloc_sub_domain(DOMD *md);
  void init_sub_domain(DOMD *md);
  void init_boundary_data(DOMD *md);
  
  int i;

  // multi_wave
  setup_mfb(&(md->mw));
  // medium
  md->n[0]=md->mw.n_0;  
  for(i=0;i<=md->MN;i++) md->kn[i]=2.0*M_PI/md->mw.lambda_0*md->n[i];
  // rotation and translation object
  if(md->th!=0.0 || vabs_d(md->tv)!=0.0) rotation_translation_obj(md->rv,md->th,md->tv,md);  
  // element constant
  init_elem_const(&(md->bd));
  // sub domain
  malloc_sub_domain(md);
  init_sub_domain(md);
  // boundary data
  init_boundary_data(md);
}

void finalize_domd(DOMD *md)
{
  void mfree_node(BOUD *bd);
  void mfree_elem(BOUD *bd);
  void mfree_sub_domain(DOMD *md);
  
  mfree_sub_domain(md);

  mfree_elem(&(md->bd));
  mfree_node(&(md->bd));
  free(md->n);
  free(md->kn);

  free_mfb(&(md->mw));

}

int domain_id(double *rt,DOMD *md)
{
  double fid_calc_solid_angle(int type,double r[4][3],int *flg);
  
  double rv[4][3],omega;
  double *Og=(double *)m_alloc2(md->MN+1,sizeof(double),"b_utils.c, domain_id()");
  int i,j,k,d,flg;

  for(d=0;d<md->MN+1;d++){
    for(i=1;i<=md->bd.sb[d].Ne;i++){
      if(md->bd.sb[d].sid[i]>0){
        // read node data
        for(j=0;j<4;j++)
          for(k=0;k<3;k++) rv[j][k]=md->bd.rn[md->bd.ed[md->bd.sb[d].sid[i]][j]][k]-rt[k];
        omega=fid_calc_solid_angle(check_element_type(md->bd.sb[d].sid[i],&(md->bd)),rv,&flg);
        if(flg<0){ // on boundary
          free(Og);
          return d;
        }
        Og[d]+=omega;
        Og[md->bd.sd[md->bd.sb[d].sid[i]]]-=omega;

      } // end if
    }
    if(d==0 && fabs(Og[d])<2.0*M_PI){ // opened region
      free(Og);
      return d;
    }
    else if(Og[d]>2.0*M_PI){ // closed region
      free(Og);
      return d;
    }
  }

  free(Og);
  return -1; // error
}

void dat_read_domd(char *fname,DOMD *md)
{
  void malloc_node(BOUD *bd);
  void malloc_elem(BOUD *bd);
  void init_elem_const(BOUD *bd);
  void malloc_sub_domain(DOMD *md);
  
  FILE *fp;
  int i,j,d,tmp;

  if((fp=fopen(fname,"rb"))==NULL){    printf("dat_read(), Failed to open the %s file.\n",fname);    exit(1);  }

  // fname
  fread(md->med_fn,sizeof(char),128,fp);
  fread(md->msh_fn,sizeof(char),128,fp);
  // rotation and translation data
  fread(md->rv,sizeof(double),3,fp);
  fread(&(md->th),sizeof(double),1,fp);
  fread(md->tv,sizeof(double),3,fp);
  // material def
  fread(&(md->MN),sizeof(int),1,fp);
  md->n =(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(), md->n"); // malloc
  md->kn=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(),md->kn"); // malloc
  fread(md->n,sizeof(double complex),md->MN+1,fp);
  fread(md->kn,sizeof(double complex),md->MN+1,fp);
  // multi_wave_nt
  fread(&(md->mw),sizeof(Bobj),1,fp);
  md->mw.bd.ipw=(Ipw *)m_alloc2(md->mw.n_ipw,sizeof(Ipw),"setup.c,dat_read(),md->mw.bd.ipw"); // malloc
  md->mw.bd.fpw=(Fpw *)m_alloc2(md->mw.n_fpw,sizeof(Fpw),"setup.c,dat_read(),md->mw.bd.fpw"); // malloc
  md->mw.bd.lgb=(LGb *)m_alloc2(md->mw.n_lgb,sizeof(LGb),"setup.c,dat_read(),md->mw.bd.lgb"); // malloc
  md->mw.bd.bsb=(Bsb *)m_alloc2(md->mw.n_bsb,sizeof(Bsb),"setup.c,dat_read(),md->mw.bd.bsb"); // malloc
  md->mw.bd.blg=(BsLGb *)m_alloc2(md->mw.n_blg,sizeof(BsLGb),"setup.c,dat_read(),md->mw.bd.blg"); // malloc
  md->mw.bd.rab=(RAb *)m_alloc2(md->mw.n_rab,sizeof(RAb),"setup.c,dat_read(),md->mw.bd.rab"); // malloc
  fread(md->mw.bd.ipw,sizeof(Ipw),md->mw.n_ipw,fp);
  fread(md->mw.bd.fpw,sizeof(Fpw),md->mw.n_fpw,fp);
  fread(md->mw.bd.lgb,sizeof(LGb),md->mw.n_lgb,fp);
  fread(md->mw.bd.bsb,sizeof(Bsb),md->mw.n_bsb,fp);
  fread(md->mw.bd.blg,sizeof(BsLGb),md->mw.n_blg,fp);
  fread(md->mw.bd.rab,sizeof(RAb),md->mw.n_rab,fp);
  setup_mfb(&(md->mw)); // setup multi_wave
  // BOUD
  fread(&(md->bd.Nn),sizeof(int),1,fp);
  fread(&(md->bd.Ne),sizeof(int),1,fp);
  malloc_node(&(md->bd)); // malloc
  malloc_elem(&(md->bd)); // malloc
  for(i=0;i<=md->bd.Nn;i++) fread(md->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fread(md->bd.ed[i],sizeof(int),4,fp);
  fread(&(md->bd.NN),sizeof(int),1,fp);
  for(i=0;i<=md->bd.Ne;i++) fread(md->bd.eni[i],sizeof(int),4,fp);
  fread(md->bd.md,sizeof(int),md->bd.Ne+1,fp);
  fread(md->bd.sd,sizeof(int),md->bd.Ne+1,fp);
  fread(md->bd.gd,sizeof(int),md->bd.Ne+1,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.wen[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd. Ui[i][j],sizeof(double complex),4,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.dUi[i][j],sizeof(double complex),4,fp);
  init_elem_const(&(md->bd)); // setup
  // sub domain data
  malloc_sub_domain(md); // malloc
  for(d=0;d<=md->MN;d++){
    fread(&tmp,sizeof(int),1,fp);
    fread(md->bd.sb[d].sid,sizeof(int),md->bd.sb[d].Ne+1,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(md->bd.sb[d]. U[i][j],sizeof(double complex),4,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(md->bd.sb[d].dU[i][j],sizeof(double complex),4,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(md->bd.sb[d]. E[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(md->bd.sb[d].dE[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(md->bd.sb[d]. H[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(md->bd.sb[d].dH[i][j],sizeof(double complex),3,fp);
  }

  fclose(fp);
}

void dat_write_domd(char *fname,DOMD *md)
{
  FILE *fp;
  int i,j,d;

  if((fp=fopen(fname,"wb"))==NULL){    printf("dat_write(), Failed to create the %s file.\n",fname);    exit(1);  }

  // fname
  fwrite(md->med_fn,sizeof(char),128,fp);
  fwrite(md->msh_fn,sizeof(char),128,fp);
  // rotation and translation data
  fwrite(md->rv,sizeof(double),3,fp);
  fwrite(&(md->th),sizeof(double),1,fp);
  fwrite(md->tv,sizeof(double),3,fp);
  // material def
  fwrite(&(md->MN),sizeof(int),1,fp);
  fwrite(md->n,sizeof(double complex),md->MN+1,fp);
  fwrite(md->kn,sizeof(double complex),md->MN+1,fp);
  // multi_fbeam
  fwrite(&(md->mw),sizeof(Bobj),1,fp);
  fwrite(md->mw.bd.ipw,sizeof(Ipw),md->mw.n_ipw,fp);
  fwrite(md->mw.bd.fpw,sizeof(Fpw),md->mw.n_fpw,fp);
  fwrite(md->mw.bd.lgb,sizeof(LGb),md->mw.n_lgb,fp);
  fwrite(md->mw.bd.bsb,sizeof(Bsb),md->mw.n_bsb,fp);
  fwrite(md->mw.bd.blg,sizeof(BsLGb),md->mw.n_blg,fp);
  fwrite(md->mw.bd.rab,sizeof(RAb),md->mw.n_rab,fp);
  // BOUD
  fwrite(&(md->bd.Nn),sizeof(int),1,fp);
  fwrite(&(md->bd.Ne),sizeof(int),1,fp);
  for(i=0;i<=md->bd.Nn;i++) fwrite(md->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.ed[i],sizeof(int),4,fp);
  fwrite(&(md->bd.NN),sizeof(int),1,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.eni[i],sizeof(int),4,fp);
  fwrite(md->bd.md,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.sd,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.gd,sizeof(int),md->bd.Ne+1,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.wen[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd. Ui[i][j],sizeof(double complex),4,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.dUi[i][j],sizeof(double complex),4,fp);
  // sub domain data
  for(d=0;d<=md->MN;d++){
    fwrite(&(md->bd.sb[d].Ne),sizeof(int),1,fp);
    fwrite(md->bd.sb[d].sid,sizeof(int),md->bd.sb[d].Ne+1,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.sb[d]. U[i][j],sizeof(double complex),4,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.sb[d].dU[i][j],sizeof(double complex),4,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.sb[d]. E[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.sb[d].dE[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.sb[d]. H[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=md->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.sb[d].dH[i][j],sizeof(double complex),3,fp);
  }

  fclose(fp);
}

void output_node_particles(char *fname,DOMD *md)
{
  FILE *fp;
  int s1,s2,i,j;
  char *sd,fo[128]="";

  sd=strrchr(fname,'.');
  if(sd==NULL){ // no file extension
    sprintf(fo,"%s.particles",fname);
  }
  else {
    s1=strlen(fname);
    s2=strlen(sd);
    strncpy(fo,fname,s1-s2);
    sprintf(fo,"%s.particles",fo);
  }
  
  if((fp=fopen(fo,"wt"))==NULL){    printf("Can not open the %s file.\n",fo);    exit(1);  }
  fprintf(fp,"# x y z object_id\n");
  
  for(i=1;i<=md->bd.Ne;i++){
    for(j=0;j<4;j++){
      fprintf(fp,"%15.14e %15.14e %15.14e %d\n",md->bd.ren[i][j][0],md->bd.ren[i][j][1],md->bd.ren[i][j][2],0);
    }
  }

  fclose(fp);
}

/////////////////////////////////////////////////////////////////////////////
void rotation_translation_obj(double *rv,double th,double *tv,DOMD *md)
{
  double ct,st,r[3],M[9],nv[3];
  size_t s,i;

  nv[0]=rv[0];
  nv[1]=rv[1];
  nv[2]=rv[2];
  vuni_d(nv);

  // rotation matrix
  st=sin(th);
  ct=cos(th);
  M[0]=ct+nv[0]*nv[0]*(1.0-ct);
  M[1]=nv[0]*nv[1]*(1.0-ct)-nv[2]*st;
  M[2]=nv[2]*nv[0]*(1.0-ct)+nv[1]*st;
  M[3]=nv[0]*nv[1]*(1.0-ct)+nv[2]*st;
  M[4]=ct+nv[1]*nv[1]*(1.0-ct);
  M[5]=nv[1]*nv[2]*(1.0-ct)-nv[0]*st;
  M[6]=nv[2]*nv[0]*(1.0-ct)-nv[1]*st;
  M[7]=nv[1]*nv[2]*(1.0-ct)+nv[0]*st;
  M[8]=ct+nv[2]*nv[2]*(1.0-ct);

  for(s=1;s<=md->bd.Nn;s++){
    for(i=0;i<3;i++) r[i]=M[3*i+0]*md->bd.rn[s][0]+M[3*i+1]*md->bd.rn[s][1]+M[3*i+2]*md->bd.rn[s][2]+tv[i];
    for(i=0;i<3;i++) md->bd.rn[s][i]=r[i];
  }
}

void filename_chk(int argc,char **argv)
{
  if(argc!=4 && argc!=11){
    printf("This program needs command line arguments as follows.\n");
    printf("%s medium_datafile_name mesh_datafile_name output_datafile_name [rv_x rv_y rv_z theta tr_x tr_y tr_z ](optional)\n",argv[0]);
    printf("rv : vector defining rotation axis, theta : rotation angle ( using Rodrigues' rotation formula ), tr : translation vector\n"); 
    printf("Exit...\n");
    exit(0);
  }
}

void read_medium_data(char *med_fn,DOMD *md)
{
  FILE *fp;
  double td,td2;
  char buf[256]="";
  int i,ti;

  if((fp=fopen(med_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",med_fn);    exit(1);  }
  strcpy(md->med_fn,med_fn);
  fgets(buf,256,fp);
  fgets(buf,256,fp);
  fscanf(fp,"%d\n",&ti);  md->MN=ti;
  md->n=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(), md->n");
  md->kn=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(),md->kn");
  fgets(buf,256,fp);
  for(i=1;i<=md->MN;i++){
    fscanf(fp,"%lf",&td);
    fscanf(fp,"%lf",&td2); md->n[i]=td+td2*I;
  }
  fclose(fp);
}

void print_medium_data(DOMD *md)
{
  int i;
  printf("-- medium data --\n");
  printf("medium data file name                  : %s\n",md->med_fn);
  for(i=1;i<=md->MN;i++){
    printf("medium (domain) id %2d refractive index :%8.7g + %8.7gI\n",i,creal(md->n[i]),cimag(md->n[i]));
  }
}

void read_mesh_data(char *msh_fn,DOMD *md)
{
  void malloc_node(BOUD *bd);
  void malloc_elem(BOUD *bd);
  
  FILE *fp;
  char buf[256]="";
  double td;
  int ti,i,j,ti2,etype,tmpi,nc;

  if((fp=fopen(msh_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",msh_fn);    exit(1);  }
  strcpy(md->msh_fn,msh_fn);
  fgets(buf,256,fp);
  fscanf(fp,"%lf",&td);
  // check file version
  if(td<MSHVER){
    printf("This program supports mesh file version %g later. Reading file version is %g. Exit...\n",MSHVER,td);
    fclose(fp);
    exit(1);
  }
  fscanf(fp,"%d",&ti);
  //check data format
  if(ti!=MSHASCI){
    printf("This program supports 'ASCII' data format mesh file. Exit...\n");
    fclose(fp);
    exit(1);
  }
  fscanf(fp,"%d\n",&ti);
  //check data precision
  if(ti!=MSHPREC){
    printf("This program supports double precision mesh data. Exit...\n");
    fclose(fp);
    exit(1);
  }
  fgets(buf,256,fp);
  fgets(buf,256,fp);

  fscanf(fp,"%d\n",&ti); 
  md->bd.Nn=ti;
  malloc_node(&(md->bd));
  for(i=1;i<=md->bd.Nn;i++){
    fscanf(fp,"%d",&ti);    if(ti!=i)       printf("bad id %d\n",ti);
    fscanf(fp,"%lf",&td);    md->bd.rn[i][0]=td;
    fscanf(fp,"%lf",&td);    md->bd.rn[i][1]=td;
    fscanf(fp,"%lf\n",&td); md->bd.rn[i][2]=td;
  }
  fgets(buf,256,fp);
  fgets(buf,256,fp);

  fscanf(fp,"%d",&ti);
  md->bd.Ne=ti/2; 
  malloc_elem(&(md->bd));

  nc=0;
  for(i=1;i<=md->bd.Ne;i++){
    // element id
    fscanf(fp,"%d",&ti);
    if(ti!=i*2-1){
      printf("bad id :%d. Exit...\n",ti);
      exit(1);
    }
    // element type
    fscanf(fp,"%d",&ti);
    etype=ti;
    if(ti!=ELT3 && ti!=ELT4){
      printf("bad element type. element type must be %d or %d. Exit...\n",ELT3,ELT4);
      exit(1);
    }
    // number of tags
    fscanf(fp,"%d",&ti);
    for(j=0;j<ti;j++){
      fscanf(fp,"%d",&ti2);
      if(j==0){ // domain id ( gmsh physical entity)
        if(ti2==OPENDID) ti2=0;
        if(md->MN>=ti2) md->bd.md[i]=ti2;
        else {
          printf("domain id %d is not defined medium data. check domain and medium data. exit..\n",ti2);
          exit(1);
        }
      }
      else if(j==1){ // group id ( elementary geometrical entity )
        md->bd.gd[i]=ti2;
      }
    }
    // node id
    fscanf(fp,"%d",&ti);  md->bd.ed[i][0]=ti;
    fscanf(fp,"%d",&ti);  md->bd.ed[i][1]=ti;
    fscanf(fp,"%d",&ti);  md->bd.ed[i][2]=ti;
    if(etype==ELT3) md->bd.ed[i][3]=0;
    else {
      fscanf(fp,"%d",&ti);  md->bd.ed[i][3]=ti;
    }
    // element node id
    if(etype==ELT3){
      md->bd.eni[i][0]=nc++;
      md->bd.eni[i][1]=nc++;
      md->bd.eni[i][2]=nc++;
      md->bd.eni[i][3]=-1;
    }
    else {
      md->bd.eni[i][0]=nc++;
      md->bd.eni[i][1]=nc++;
      md->bd.eni[i][2]=nc++;
      md->bd.eni[i][3]=nc++;
    }

    // element id
    fscanf(fp,"%d",&ti);
    if(ti!=i*2){
      printf("bad id :%d. Exit...\n",ti);
      exit(1);
    }
    // element type
    fscanf(fp,"%d",&ti);
    etype=ti;
    if(ti!=ELT3 && ti!=ELT4){
      printf("bad element type. element type must be %d or %d. Exit...\n",ELT3,ELT4);
      exit(1);
    }
    // number of tags
    fscanf(fp,"%d",&ti);
    for(j=0;j<ti;j++){
      fscanf(fp,"%d",&ti2);
      if(j==0){ // domain id
        if(ti2==OPENDID) ti2=0;
        if(md->MN>=ti2) md->bd.sd[i]=ti2;
        else {
          printf("domain id %d is not defined medium data. check domain and medium data! exit..\n",ti2);
          exit(1);
        }
      }
    }
    // check node id
    if(etype==ELT3){
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][0]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][2]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][1]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
    }
    else {
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][0]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][3]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][2]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(md->bd.ed[i][1]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
    }
    // exchange open region domain to main domain
    if(md->bd.sd[i]==0){
      md->bd.sd[i]=md->bd.md[i];
      md->bd.md[i]=0;
      if(etype==ELT3){
        tmpi=md->bd.ed[i][1];
        md->bd.ed[i][1]=md->bd.ed[i][2];
        md->bd.ed[i][2]=tmpi;
      }
      else {
        tmpi=md->bd.ed[i][3];
        md->bd.ed[i][3]=md->bd.ed[i][1];
        md->bd.ed[i][1]=tmpi;
      }
    }
  }
  fclose(fp);
  md->bd.NN=nc;
 
}

void print_mesh_data(DOMD *md)
{
  printf("-- mesh data --\n");
  printf("mesh data file name    : %s\n",md->msh_fn);
  printf("node number            : %8d\n",md->bd.Nn);
  printf("defined element number : %8d\n",md->bd.Ne*2);
}

void malloc_node(BOUD *bd)
{
  int i,N=bd->Nn;

  bd->rn=(double **)m_alloc2(N+1,sizeof(double*),"setup.c, malloc_node(), bd->rn");
  for(i=0;i<=N;i++){
    bd->rn[i]=(double *)m_alloc2(3,sizeof(double),"setup.c, malloc_node(), bd->rn[i]");
  }
}

void mfree_node(BOUD *bd)
{
  int i,N=bd->Nn;

  for(i=0;i<=N;i++) free(bd->rn[i]);
  free(bd->rn);
  bd->Nn=0;
}

void malloc_elem(BOUD *bd)
{
  int i,j,Ne=bd->Ne;

  bd->ed =(int **)m_alloc2(Ne+1,sizeof(int *),"setup.c, malloc_elem(), bd->ed");
  bd->eni=(int **)m_alloc2(Ne+1,sizeof(int *),"setup.c, malloc_elem(), bd->eni");
  for(i=0;i<=Ne;i++){
    bd->ed [i]=(int *)m_alloc2(4,sizeof(int ),"setup.c, malloc_elem(), bd->ed[i]");
    bd->eni[i]=(int *)m_alloc2(4,sizeof(int ),"setup.c, malloc_elem(), bd->eni[i]");
  }

  bd->md=(int *)m_alloc2(Ne+1,sizeof(int),"setup.c, malloc_elem(), bd->md");
  bd->sd=(int *)m_alloc2(Ne+1,sizeof(int),"setup.c, malloc_elem(), bd->sd");
  bd->gd=(int *)m_alloc2(Ne+1,sizeof(int),"setup.c, malloc_elem(), bd->gd");

  // element constant
  bd->cr =(double ***)m_alloc2(Ne+1,sizeof(double **),"setup.c, malloc_elem(), bd->cr");
  bd->cw =(double ***)m_alloc2(Ne+1,sizeof(double **),"setup.c, malloc_elem(), bd->cw");
  bd->ren=(double ***)m_alloc2(Ne+1,sizeof(double **),"setup.c, malloc_elem(), bd->ren");
  bd->wen=(double ***)m_alloc2(Ne+1,sizeof(double **),"setup.c, malloc_elem(), bd->wen");
  bd-> Ui=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup.c, malloc_elem(), bd->Ui");
  bd->dUi=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup.c, malloc_elem(), bd->dUi");
  for(i=0;i<=Ne;i++){
    bd->cr[i]=(double **)m_alloc2(3,sizeof(double *),"setup.c, malloc_elem(), bd->cr[i]");
    bd->cw[i]=(double **)m_alloc2(3,sizeof(double *),"setup.c, malloc_elem(), bd->cw[i]");
    for(j=0;j<3;j++){
      bd->cr[i][j]=(double *)m_alloc2(4,sizeof(double),"setup.c, malloc_elem(), bd->cr[i][j]");
      bd->cw[i][j]=(double *)m_alloc2(3,sizeof(double),"setup.c, malloc_elem(), bd->cw[i][j]");
    }

    bd->ren[i]=(double **)m_alloc2(4,sizeof(double *),"setup.c, malloc_elem(), bd->ren[i]");
    bd->wen[i]=(double **)m_alloc2(4,sizeof(double *),"setup.c, malloc_elem(), bd->wen[i]");
    bd-> Ui[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_elem(), bd->Ui[i]");
    bd->dUi[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_elem(), bd->dUi[i]");
    for(j=0;j<4;j++){
      bd->ren[i][j]=(double *)m_alloc2(3,sizeof(double),"setup.c, malloc_elem(), bd->ren[i][j]");
      bd->wen[i][j]=(double *)m_alloc2(3,sizeof(double),"setup.c, malloc_elem(), bd->wen[i][j]");
      bd-> Ui[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup.c, malloc_elem(), bd->Ui[i][j]");
      bd->dUi[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup.c, malloc_elem(), bd->dUi[i][j]");
    }
  }
}

void mfree_elem(BOUD *bd)
{
  int i,j,Ne=bd->Ne;

  for(i=0;i<=Ne;i++){
    free(bd->ed[i]);
    free(bd->eni[i]);
  }
  free(bd->ed);
  free(bd->eni);

  free(bd->md);
  free(bd->sd);
  free(bd->gd);

  for(i=0;i<=Ne;i++){
    for(j=0;j<3;j++){
      free(bd->cr[i][j]);
      free(bd->cw[i][j]);
    }
    free(bd->cr[i]);
    free(bd->cw[i]);

    for(j=0;j<4;j++){
      free(bd->ren[i][j]);
      free(bd->wen[i][j]);
      free(bd->Ui[i][j]);
      free(bd->dUi[i][j]);
    }
    free(bd->ren[i]);
    free(bd->wen[i]);
    free(bd->Ui[i]);
    free(bd->dUi[i]);
  }
  free(bd->cr);
  free(bd->cw);
  free(bd->ren);
  free(bd->wen);
  free(bd->Ui);
  free(bd->dUi);

  bd->Ne=0;
}

void init_elem_const(BOUD *bd)
{
  int i,j,d,Ne,a,b;
  double rc[3][4];

  Ne=bd->Ne;

  // geometric constant
  for(i=1;i<=Ne;i++){
    for(d=0;d<3;d++)
      for(j=0;j<4;j++) rc[d][j]=bd->rn[bd->ed[i][j]][d];

    if(bd->ed[i][3]!=0){ // bi-linear element
      for(d=0;d<3;d++){
        bd->cr[i][d][0]=0.25*( rc[d][0]+rc[d][1]+rc[d][2]+rc[d][3]);
        bd->cr[i][d][1]=0.25*(-rc[d][0]+rc[d][1]+rc[d][2]-rc[d][3]);
        bd->cr[i][d][2]=0.25*(-rc[d][0]-rc[d][1]+rc[d][2]+rc[d][3]);
        bd->cr[i][d][3]=0.25*( rc[d][0]-rc[d][1]+rc[d][2]-rc[d][3]);
      }

      a=1;      b=2;
      bd->cw[i][0][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][0][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][0][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
      a=2;      b=0;
      bd->cw[i][1][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][1][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][1][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
      a=0;      b=1;
      bd->cw[i][2][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][2][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][2][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
    }
    else { // linear triangular element
      for(d=0;d<3;d++){
        bd->cr[i][d][0]=1.0/3.0*( rc[d][0]+rc[d][1]+rc[d][2]);
        bd->cr[i][d][1]=1.0/3.0*(-rc[d][0]+2.0*rc[d][1]-rc[d][2]);
        bd->cr[i][d][2]=1.0/sqrt(3.0)*( rc[d][2]-rc[d][0]);
        bd->cr[i][d][3]=0.0;
      }

      bd->cw[i][0][0]=( (rc[1][1]-rc[1][0])*(rc[2][2]-rc[2][0]) - (rc[2][1]-rc[2][0])*(rc[1][2]-rc[1][0]) );
      bd->cw[i][0][1]=0.0;
      bd->cw[i][0][2]=0.0;
      bd->cw[i][1][0]=( (rc[2][1]-rc[2][0])*(rc[0][2]-rc[0][0]) - (rc[0][1]-rc[0][0])*(rc[2][2]-rc[2][0]) );
      bd->cw[i][1][1]=0.0;
      bd->cw[i][1][2]=0.0;
      bd->cw[i][2][0]=( (rc[0][1]-rc[0][0])*(rc[1][2]-rc[1][0]) - (rc[1][1]-rc[1][0])*(rc[0][2]-rc[0][0]) );
      bd->cw[i][2][1]=0.0;
      bd->cw[i][2][2]=0.0;
    }
  }

  // element constant
  // gaussian quadrature node and weight
  bd->zt_44[0]=-P44_N;  bd->zt_44[1]= P44_N;  bd->zt_44[2]= P44_N;  bd->zt_44[3]=-P44_N;
  bd->et_44[0]=-P44_N;  bd->et_44[1]=-P44_N;  bd->et_44[2]= P44_N;  bd->et_44[3]= P44_N;
  bd->wt_44[0]= P44_W;  bd->wt_44[1]= P44_W;  bd->wt_44[2]= P44_W;  bd->wt_44[3]= P44_W;

  bd->zt_49[0]=-P49_N;   bd->zt_49[1]= P49_N;   bd->zt_49[2]= P49_N;
  bd->zt_49[3]=-P49_N;   bd->zt_49[4]= 0.0;     bd->zt_49[5]= P49_N;
  bd->zt_49[6]= 0.0;     bd->zt_49[7]=-P49_N;   bd->zt_49[8]= 0.0;
  bd->et_49[0]=-P49_N;   bd->et_49[1]=-P49_N;   bd->et_49[2]= P49_N;
  bd->et_49[3]= P49_N;   bd->et_49[4]=-P49_N;   bd->et_49[5]= 0.0;
  bd->et_49[6]= P49_N;   bd->et_49[7]= 0.0;     bd->et_49[8]= 0.0;
  bd->wt_49[0]= P49_W0;  bd->wt_49[1]= P49_W0;  bd->wt_49[2]= P49_W0;
  bd->wt_49[3]= P49_W0;  bd->wt_49[4]= P49_W1;  bd->wt_49[5]= P49_W1;
  bd->wt_49[6]= P49_W1;  bd->wt_49[7]= P49_W1;  bd->wt_49[8]= P49_W2;

  bd->zt_34[0]=-P34_N0;  bd->zt_34[1]= 2.0*P34_N0;   bd->zt_34[2]=-P34_N0;  bd->zt_34[3]= 0.0;
  bd->et_34[0]=-P34_N1;  bd->et_34[1]= 0.0;          bd->et_34[2]= P34_N1;  bd->et_34[3]= 0.0;
  bd->wt_34[0]= P34_W0;  bd->wt_34[1]= P34_W0;       bd->wt_34[2]= P34_W0;  bd->wt_34[3]=-P34_W1;

  bd->zt_37[0]=-P37_N0;  bd->zt_37[1]= 2.0*P37_N0;  bd->zt_37[2]=-P37_N0;
  bd->zt_37[3]= P37_N1;  bd->zt_37[4]=-2.0*P37_N1;  bd->zt_37[5]= P37_N1;  bd->zt_37[6]= 0.0;
  bd->et_37[0]=-P37_N2;  bd->et_37[1]= 0.0;         bd->et_37[2]= P37_N2;
  bd->et_37[3]= P37_N3;  bd->et_37[4]= 0.0;         bd->et_37[5]=-P37_N3;  bd->et_37[6]= 0.0;
  bd->wt_37[0]= P37_W0;     bd->wt_37[1]= P37_W0;  bd->wt_37[2]= P37_W0;
  bd->wt_37[3]= P37_W1;     bd->wt_37[4]= P37_W1;  bd->wt_37[5]= P37_W1;     bd->wt_37[6]= P37_W2;

  // gauss-legendre GLN point rule
  gauleg(-1.0,1.0,bd->xli,bd->wli,GLN);
  gauleg(-1.0,1.0,bd->xhi,bd->whi,GHN);
}

void malloc_sub_domain(DOMD *md)
{
  int *Nc,i,j,k;
  Nc=(int *)m_alloc2(md->MN+1,sizeof(int),"setup.c, malloc_sub_domain(), Nc");
  for(i=0;i<=md->MN;i++) Nc[i]=0;

  for(i=1;i<=md->bd.Ne;i++){
    Nc[md->bd.md[i]]++;
    Nc[md->bd.sd[i]]++;
  }

  md->bd.sb=(SUBD *)m_alloc2(md->MN+1,sizeof(SUBD),"setup.c, malloc_sub_domain(), md->bd.sb");
  for(i=0;i<=md->MN;i++){
    md->bd.sb[i].Ne=Nc[i];
    md->bd.sb[i].sid=(int *)m_alloc2(Nc[i]+1,sizeof(int),"setup.c, malloc_sub_domain(), md->bd.sb[i],sid");
    md->bd.sb[i].U =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup.c, malloc_sub_domain(), md->bd.sb[i].U");
    md->bd.sb[i].dU=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup.c, malloc_sub_domain(), md->bd.sb[i].dU");
    md->bd.sb[i].E =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup.c, malloc_sub_domain(), md->bd.sb[i].E");
    md->bd.sb[i].dE=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup.c, malloc_sub_domain(), md->bd.sb[i].dE");
    md->bd.sb[i].H =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup.c, malloc_sub_domain(), md->bd.sb[i].H");
    md->bd.sb[i].dH=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup.c, malloc_sub_domain(), md->bd.sb[i].dH");
    for(j=0;j<=Nc[i];j++){
      md->bd.sb[i].U [j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_sub_domain(), md->bd.sb[i].U[j]");
      md->bd.sb[i].dU[j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_sub_domain(), md->bd.sb[i].dU[j]");
      md->bd.sb[i].E [j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_sub_domain(), md->bd.sb[i].E[j]");
      md->bd.sb[i].dE[j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_sub_domain(), md->bd.sb[i].dE[j]");
      md->bd.sb[i].H [j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_sub_domain(), md->bd.sb[i].H[j]");
      md->bd.sb[i].dH[j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup.c, malloc_sub_domain(), md->bd.sb[i].dH[j]");
      for(k=0;k<4;k++){
        md->bd.sb[i].U [j][k]=(double complex *)m_alloc2(4,sizeof(double complex),"setup.c, malloc_sub_domain(), md->bd.sb[i].U[j][k]");
        md->bd.sb[i].dU[j][k]=(double complex *)m_alloc2(4,sizeof(double complex),"setup.c, malloc_sub_domain(), md->bd.sb[i].dU[j][k]");
        md->bd.sb[i].E [j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup.c, malloc_sub_domain(), md->bd.sb[i].E[j][k]");
        md->bd.sb[i].dE[j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup.c, malloc_sub_domain(), md->bd.sb[i].dE[j][k]");
        md->bd.sb[i].H [j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup.c, malloc_sub_domain(), md->bd.sb[i].H[j][k]");
        md->bd.sb[i].dH[j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup.c, malloc_sub_domain(), md->bd.sb[i].dH[j][k]");
      }
    }
  }

  free(Nc);
}

void mfree_sub_domain(DOMD *md)
{
  int i,j,k;

  for(i=0;i<=md->MN;i++){
    for(j=0;j<=md->bd.sb[i].Ne;j++){
      for(k=0;k<4;k++){
        free(md->bd.sb[i].U[j][k]); free(md->bd.sb[i].dU[j][k]);
        free(md->bd.sb[i].E[j][k]); free(md->bd.sb[i].dE[j][k]);
        free(md->bd.sb[i].H[j][k]); free(md->bd.sb[i].dH[j][k]);
      }
      free(md->bd.sb[i].U[j]);      free(md->bd.sb[i].dU[j]);
      free(md->bd.sb[i].E[j]);      free(md->bd.sb[i].dE[j]);
      free(md->bd.sb[i].H[j]);      free(md->bd.sb[i].dH[j]);
    }
    free(md->bd.sb[i].sid);
    free(md->bd.sb[i].U);    free(md->bd.sb[i].dU);
    free(md->bd.sb[i].E);    free(md->bd.sb[i].dE);
    free(md->bd.sb[i].H);    free(md->bd.sb[i].dH);
  }

  free(md->bd.sb);
}

void init_sub_domain(DOMD *md)
{
  int d,i,c;

  for(d=0;d<=md->MN;d++){
    c=1;
    for(i=1;i<=md->bd.Ne;i++){
      if(md->bd.md[i]==d){
        md->bd.sb[d].sid[c]=i;
        c++;
      }
      else if(md->bd.sd[i]==d){
        md->bd.sb[d].sid[c]=-i;
        c++;
      }
    }

  }
}

void init_boundary_data(DOMD *md)
{
  double complex h[3],dh[3];
  double cr[3][4],cw[3][3],r[3],w[3];
  int i,j,l,m;

  // element node and incident field data
  for(i=1;i<=md->bd.Ne;i++){
    for(l=0;l<3;l++){
      for(m=0;m<4;m++) cr[l][m]=md->bd.cr[i][l][m];
      for(m=0;m<3;m++) cw[l][m]=md->bd.cw[i][l][m];
    }

    if(md->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        lit_rw_zeta_eta(r,w,md->bd.zt_34[j],md->bd.et_34[j],cr,cw);
        for(l=0;l<3;l++){
          md->bd.ren[i][j][l]=r[l];
          md->bd.wen[i][j][l]=w[l];
        }
        vuni_d(w);
        calc_mfb_EH_dv(md->bd.Ui[i][j],h,md->bd.dUi[i][j],dh,r,w,&(md->mw));
        md->bd.Ui[i][j][3]=0.0;
        md->bd.dUi[i][j][3]=0.0;
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        bil_rw_zeta_eta(r,w,md->bd.zt_44[j],md->bd.et_44[j],cr,cw);
        for(l=0;l<3;l++){
          md->bd.ren[i][j][l]=r[l];
          md->bd.wen[i][j][l]=w[l];
        }
        vuni_d(w);
        calc_mfb_EH_dv(md->bd.Ui[i][j],h,md->bd.dUi[i][j],dh,r,w,&(md->mw));
        md->bd.Ui[i][j][3]=0.0;
        md->bd.dUi[i][j][3]=0.0;
      }
    }
  }

}

double fid_calc_solid_angle(int type,double r[4][3],int *flg)
{
  double n01[3],n12[3],n20[3],nt[3],sa,ca;
  double Omega,a0,a1,a2,D;
  int st,i;

  Omega=0.0;
  *flg=0;

  st=0;
  if(type==ELT3) for(i=0;i<3;i++)  st+=vuni_d(r[i]);
  else for(i=0;i<4;i++) st+=vuni_d(r[i]);
  if(st<0){
    *flg=-1;    return 0; // on boundary
  }

  vcrs_d(n01,r[0],r[1]);
  vcrs_d(n12,r[1],r[2]);
  vcrs_d(n20,r[2],r[0]);

  vcrs_d(nt,n01,n20);
  sa=vabs_d(nt);
  ca=-vdot_d(n01,n20);
  a0=atan2(sa,ca);

  vcrs_d(nt,n12,n01);
  sa=vabs_d(nt);
  ca=-vdot_d(n12,n01);
  a1=atan2(sa,ca);

  vcrs_d(nt,n20,n12);
  sa=vabs_d(nt);
  ca=-vdot_d(n20,n12);
  a2=atan2(sa,ca);

  D=vdot_d(r[0],n12);
  if(D>0.0) Omega+=a0+a1+a2-M_PI;
  else if(D<0.0) Omega-=a0+a1+a2-M_PI;
  else { // on boundary
    *flg=-1;    return 0;
  }

  if(ELT4==type){
    n01[0]=-n20[0];    n01[1]=-n20[1];    n01[2]=-n20[2];

    vcrs_d(n12,r[2],r[3]);
    vcrs_d(n20,r[3],r[0]);

    vcrs_d(nt,n01,n20);
    sa=vabs_d(nt);
    ca=-vdot_d(n01,n20);
    a0=atan2(sa,ca);

    vcrs_d(nt,n12,n01);
    sa=vabs_d(nt);
    ca=-vdot_d(n12,n01);
    a1=atan2(sa,ca);

    vcrs_d(nt,n20,n12);
    sa=vabs_d(nt);
    ca=-vdot_d(n20,n12);
    a2=atan2(sa,ca);

    D=vdot_d(r[0],n12);
    if(D>0.0) Omega+=a0+a1+a2-M_PI;
    else if(D<0.0) Omega-=a0+a1+a2-M_PI;
    else { // on boundary
      *flg=-1;    return 0;
    }
  }
  return Omega;
}


