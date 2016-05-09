#include "logger.h"

void log(const char* where, const char* what, bool terminate, int importance){
  FILE *f = fopen("gneural_network.log", "a");
  if(f == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  for(int i=0; i<importance; i++){
    fprintf(f,"*");
  }
  fprintf(f," ",where,": ",what);
  if(terminate){ fprintf(f,"\n"); }

  fclose(f);
  
}
