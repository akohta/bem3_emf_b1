#include "bem3_emf_b1.h"

int main(int argc,char **argv)
{
  DOMD md;

  read_domd(argc,argv,&md);
  print_domd(&md);
  //print_domd_mksa(&md);
  initialize_domd(&md);
  output_node_particles(argv[3],&md);

  solve_bieq(&md);
  dat_write_domd(argv[3],&md);

  finalize_domd(&md);
  return 0;  
}
