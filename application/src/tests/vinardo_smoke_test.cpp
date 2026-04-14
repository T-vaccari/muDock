#include <cassert>

#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/assign_autodock_types.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
// To compile&build
// cmake -S /workspaces/muDock -B /workspaces/muDock/build/linux-debug -DMUDOCK_ENABLE_TEST=ON
// cmake --build /workspaces/muDock/build/linux-debug --target test_vinardo_smoke

//To run 
// /workspaces/muDock/build/linux-debug/application/tests/test_vinardo_smoke



int main() {
  //Currently building on another branch a pipeline to test the atom typing 
  // and pre processing against the reference implementation(smina)
  

  return 0;
}
