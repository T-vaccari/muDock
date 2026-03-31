#include <mudock/chem/vinardo_type.hpp>
#include <mudock/chem/vinardo_layer.hpp>

// Used to compile the new files: clang++ -std=c++20 -I mudock/include -c smoketest.cpp
// Simple smoke test to compile the new files.
int main() { 
   //Instantiate a molecule
   mudock::molecule<mudock::static_containers> mol;

   // Alloca spazio per 3 atomi e 0 legami
   mol.resize(3, 0);

   // Riempi la molecola assegnando i tipi di autodock (usando OA per Ossigeno)
   mol.autodock_type(0) = mudock::autodock_ff::C;
   mol.autodock_type(1) = mudock::autodock_ff::N;
   mol.autodock_type(2) = mudock::autodock_ff::OA;
   
   //Instantiate a vinardo layer
   mudock::vinardo_layer<mudock::static_containers> vinardo_layer(mol);


   
   return 0; 
   
}
