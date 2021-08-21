#include "bem3_emf_b1.h"

int main(int argc,char *argv[])
{
  DOMD md;
  double complex E[3],H[3];
  double F[3],N[3],r[3],rc[3];
  int type;

  dat_read_domd(argv[1],&md); // read datafile outputed by 'd3b1_bv_solver'
  print_domd(&md);            // print data 
  //print_domd_mksa(&md);     // for print in MKSA system of units

  r[0]=0.3;                   // set x-coordinate
  r[1]=0.0;                   // set y-coordinate
  r[2]=1.1;                   // set z-coordinate 
  type=0;                     // selsect 4 point Gauss-Legendre 
  EH_t(E,H,r,type,&md);       // calclation of total field ( add incident field to scattered field )
  //EH_mEMP_t(E,H,r,type,&md);// for far-field ( error is smaller than EH_t() in far-field )
  printf("Total electromagnetic field at r=( % g,% g,% g )\n",r[0],r[1],r[2]);
  printf("type=%d setting ( 4 point Gauss-Legendre )\n",type);
  printf("Ex = % 15.14e %+15.14e I \n",creal(E[0]),cimag(E[0]));
  printf("Ey = % 15.14e %+15.14e I \n",creal(E[1]),cimag(E[1])); 
  printf("Ez = % 15.14e %+15.14e I \n",creal(E[2]),cimag(E[2]));
  printf("Hx = % 15.14e %+15.14e I \n",creal(H[0]),cimag(H[0]));
  printf("Hy = % 15.14e %+15.14e I \n",creal(H[1]),cimag(H[1])); 
  printf("Hz = % 15.14e %+15.14e I \n",creal(H[2]),cimag(H[2]));  
  type=1;                     // selsect 9 point(quadrangular element) or 7 point(triangular element) Gauss-Legendre 
  EH_t(E,H,r,type,&md);       // calclation of total field ( add incident field to scattered field )
  //EH_mEMP_t(E,H,r,type,&md);// for far-field. error is smaller than EH_t() in far-field
  printf("type=%d setting ( 9 or 7 point Gauss-Legendre )\n",type);
  printf("Ex = % 15.14e %+15.14e I \n",creal(E[0]),cimag(E[0]));
  printf("Ey = % 15.14e %+15.14e I \n",creal(E[1]),cimag(E[1])); 
  printf("Ez = % 15.14e %+15.14e I \n",creal(E[2]),cimag(E[2]));
  printf("Hx = % 15.14e %+15.14e I \n",creal(H[0]),cimag(H[0]));
  printf("Hy = % 15.14e %+15.14e I \n",creal(H[1]),cimag(H[1])); 
  printf("Hz = % 15.14e %+15.14e I \n",creal(H[2]),cimag(H[2]));  
   
  rc[0]=0.0;                  // set x coordinate of rotation center 
  rc[1]=0.0;                  // set y coordinate of rotation center 
  rc[2]=0.0;                  // set z coordinate of rotation center 
  printf("\nRadiation force and torque\n");
  printf("type=0 setting ( 4 point Gauss-Legendre ) \n"); 
  type=0;                     // select 4 point Gauss-Legendre
  force_FN(F,N,rc,type,&md);
  printf("F = (% 15.14e,% 15.14e,% 15.14e)\n",F[0],F[1],F[2]);
  printf("N = (% 15.14e,% 15.14e,% 15.14e), center of rotation (% g,% g,% g)\n",N[0],N[1],N[2],rc[0],rc[1],rc[2]); 
  printf("type=1 setting ( 9 or 7 point Gauss-Legendre ) \n");
  type=1;                     // select 9 or 7 point Gauss-Legendre
  force_FN(F,N,rc,type,&md);
  printf("F = (% 15.14e,% 15.14e,% 15.14e)\n",F[0],F[1],F[2]);
  printf("N = (% 15.14e,% 15.14e,% 15.14e), center of rotation (% g,% g,% g)\n",N[0],N[1],N[2],rc[0],rc[1],rc[2]); 
  
  finalize_domd(&md);
  return 0;
}
