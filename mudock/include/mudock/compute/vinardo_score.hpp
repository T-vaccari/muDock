#pragma once
#include <mudock/type_alias.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>




//Here we should expose the interface of the vinardo scoring function
//the scoring function should return a number the for now we assume being a fp

namespace mudock {

   fp_type vinardo_score(vinardo_layer<dynamic_containers>& protein,
                        vinardo_layer<static_containers>& ligand,
                        vinardo_preprocessed_pairs& preprocessed_pairs);

} // namespace mudock