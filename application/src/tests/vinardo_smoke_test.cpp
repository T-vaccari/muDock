#include <cassert>
#include <cmath>
#include <type_traits>

#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/assign_autodock_types.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/compute/vinardo_affinity.hpp>
#include <mudock/compute/vinardo_scoring_function.hpp>
#include <mudock/chem/vinardo_smina_helpers.hpp>
#include <mudock/compute/vinardo_score.hpp>
#include <mudock/compute/vinardo_score_kernel.hpp>
#include <mudock/compute/pipeline.hpp>



// To compile&build
// cmake -S /workspaces/muDock -B /workspaces/muDock/build/linux-debug -DMUDOCK_ENABLE_TEST=ON
// cmake --build /workspaces/muDock/build/linux-debug --target test_vinardo_smoke

//To run 
// /workspaces/muDock/build/linux-debug/application/tests/test_vinardo_smoke



int main() {
  static_assert(std::is_same_v<mudock::vinardo_score_pipeline,
                               mudock::scoring_pipeline<mudock::vinardo_score>>);
  static_assert(std::is_same_v<mudock::genetic_vinardo_pipeline,
                               mudock::genetic_scoring_pipeline<mudock::vinardo_score>>);

  return 0;
}
