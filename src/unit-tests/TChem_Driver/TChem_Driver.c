
#include <stdio.h>
#include <TChem_Driver.hpp>

int main(int argc, char* argv[]) {
  printf("TChem driver initialized successfully.\n");

   int len;
   int nBatch = 1;
   initialize(argv[1], argv[2], argv[3], nBatch);

   len = TChem_getLengthOfStateVector();
   double *state = (double *)calloc(len, sizeof(double));

   len = TChem_getNumberConcentrationVectorSize();
   double *num_conc = (double *)calloc(len, sizeof(double));

   state[0] = 0.0;
   state[1] = 81060;
   state[2] = 270.5;

   FILE *file = fopen("gas_species_names.txt", "w");

   char SpecName[100]; 
   for (int i = 0; i < TChem_getNumberOfSpecies(); i++){
    TChem_getSpeciesName(&i, SpecName, 100);
    fprintf(file, "%s\n", SpecName);
   }

   state[50] = 0.06;
   state[51] = 1e-06;
   state[36] = 0.001;
   state[81] = 1e8;

   int i_spec;
   i_spec = 3 + TChem_getNumberOfSpecies();
   for (int i = 0; i < len; i++){
     num_conc[i] = 10000.0;
     state[i_spec] = 1e-08;
     i_spec = i_spec + 5;
   }
 
   TChem_setStateVector(state, 0);
   TChem_setNumberConcentrationVector(num_conc, 0);

   FILE *aero_prop_file = fopen("aero_species_props.txt", "w");

   double density, mw, kappa;
   for (int i = 0; i< TChem_getNumberOfAeroSpecies(); i++){
      density = TChem_getAerosolSpeciesDensity(&i);
      mw = TChem_getAerosolSpeciesMW(&i);
      kappa = TChem_getAerosolSpeciesKappa(&i);
      TChem_getAerosolSpeciesName(&i, SpecName, 100);
      fprintf(aero_prop_file, "%s %f %f %f\n", SpecName, density, mw, kappa);
   }

   double del_t = 30.0;

   TChem_doTimestep(del_t);
  
   FILE *output_file = fopen("output.txt", "w");
 
   TChem_getStateVector(state, 0);
   for (int i = 0; i < TChem_getLengthOfStateVector(); i++){
      fprintf(output_file, "%e \n", state[i]);
   }

   finalize();
   printf("TChem driver finalized.\n");
   return 0;
}
