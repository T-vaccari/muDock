#pragma once
#include <mudock/type_alias.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>




//Here we should expose the interface of the vinardo scoring function
//the scoring function should return a number the for now we assume being a fp
//We also maintains some API used during the test to receive more usefull data rather than the final score, such as the breakdown of the
// score into different terms, and the contribution of each pair to the score, etc.
namespace mudock {

   //Debug Structures
   struct vinardo_pair_terms {
      fp_type gauss;
      fp_type repulsion;
      fp_type hydrophobic;
      fp_type hbond;
      fp_type weighted_gauss;
      fp_type weighted_repulsion;
      fp_type weighted_hydrophobic;
      fp_type weighted_hbond;
      fp_type weighted_energy;
   };

   struct vinardo_score_breakdown {
      fp_type protein_ligand;
      fp_type ligand_ligand;
      fp_type total;
   };


   //Debug API's
   vinardo_pair_terms compute_vinardo_pair_terms(fp_type surface_distance,
                                                 bool hydrophobic_possible,
                                                 bool hbond_possible);

   vinardo_score_breakdown compute_vinardo_score_breakdown(vinardo_layer<dynamic_containers>& protein,
                                                           vinardo_layer<static_containers>& ligand,
                                                           vinardo_preprocessed_pairs& preprocessed_pairs);


   //Fast API's
   fp_type compute_vinardo_pair_energy(fp_type surface_distance,
                                       bool hydrophobic_possible,
                                       bool hbond_possible);

   fp_type vinardo_score(vinardo_layer<dynamic_containers>& protein,
                        vinardo_layer<static_containers>& ligand,
                        vinardo_preprocessed_pairs& preprocessed_pairs);

} // namespace mudock
