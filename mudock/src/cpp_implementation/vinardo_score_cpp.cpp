#include <mudock/type_alias.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <cmath>
#include <algorithm>



namespace mudock {


   inline fp_type gauss(const fp_type d) {
      //Const parameters of the function
      constexpr fp_type s1 = fp_type{0.8};
      
      const fp_type x = d / s1;
      return std::exp(-(x * x));
   }

   inline fp_type repulsion(const fp_type d) {
      //If it's greater than 0 there is no repulsion, otherwise we have a quadratic repulsion
      const fp_type clamped = std::min(d, fp_type{0});
      return clamped * clamped;
   }
   inline fp_type hydrophobic(const fp_type d) {
      constexpr fp_type p1 = fp_type{0};
      constexpr fp_type p2 = fp_type{2.5};

      if (d <= p1) {
         return fp_type{1};
      } else if (d < p2) {
         return fp_type{p2 - d};
      }else {
         return fp_type{0};
      }
  }

   inline fp_type hbond(const fp_type d) {
      //Const parameters of the function
      constexpr fp_type h1 = -0.6;

      if( d <= h1){
         return fp_type{1};
      } else if (d < fp_type{0}){
         return d / (h1);
      }else {
         return fp_type{0};
      }
   }


   fp_type vinardo_score(vinardo_layer<dynamic_containers>& protein,
                        vinardo_layer<static_containers>& ligand,
                        vinardo_preprocessed_pairs& preprocessed_pairs){

      fp_type score = 0.0;
      fp_type pl_score = fp_type{0};
      fp_type ll_score = fp_type{0};
      auto& pl_pairs = preprocessed_pairs.protein_ligand_pairs;
      auto& ll_pairs = preprocessed_pairs.ligand_ligand_pairs;


      //Now we have to process all the couples, the vinardo scoring function is characterized by being
      //an averaged sum of four terms, more precisely given a couple i,j we have:
      // f_{i,j} = w1 * gauss(d) + w2 * repulsion(d) + w3 * hydrophobic + w4 * hbond
      // where d is defined as surface distance and can be obtained by:
      // d = r_i_j - (R_i + R_j)
      //Where thr sum R_i and R_j is already precoumped for every couple and the distance between the two atoms has to be calculated

      //Here we define the costants weights of the scoring function
      constexpr fp_type w1 = -0.045;
      constexpr fp_type w2 = 0.800;
      constexpr fp_type w3 = -0.035;
      constexpr fp_type w4 = -0.600;

      //Now we can process the protein-ligand pairs

      for(const auto& pair: pl_pairs){
         //I need to retrieve the distance between the two atoms
         
         const fp_type dx = protein().x(pair.protein_atom_idx) - ligand().x(pair.ligand_atom_idx);
         const fp_type dy = protein().y(pair.protein_atom_idx) - ligand().y(pair.ligand_atom_idx);
         const fp_type dz = protein().z(pair.protein_atom_idx) - ligand().z(pair.ligand_atom_idx);
         const fp_type r_i_j = std::sqrt(dx * dx + dy * dy + dz * dz);
         
         const fp_type d = r_i_j - pair.radius_sum;
         
         //This skip is done is smina, but it's not grounded on the paper
         if (r_i_j >= fp_type{8}) {
            continue;
         }

         //We use this to avoid the if in the inner loop 
         const fp_type hydro_mask = static_cast<fp_type>(pair.hydrophobic_possible);
         const fp_type hbond_mask = static_cast<fp_type>(pair.hbond_possible);

         pl_score += w1 * gauss(d)+ w2 * repulsion(d) + hydro_mask * w3 * hydrophobic(d)+ hbond_mask * w4 * hbond(d);
      }


      for(const auto& pair: ll_pairs){
         //I need to retrieve the distance between the two atoms
         
         const fp_type dx = ligand().x(pair.ligand_atom_i_idx) - ligand().x(pair.ligand_atom_j_idx);
         const fp_type dy = ligand().y(pair.ligand_atom_i_idx) - ligand().y(pair.ligand_atom_j_idx);
         const fp_type dz = ligand().z(pair.ligand_atom_i_idx) - ligand().z(pair.ligand_atom_j_idx);
         const fp_type r_i_j = std::sqrt(dx * dx + dy * dy + dz * dz);
         
         const fp_type d = r_i_j - pair.radius_sum;
         
         //This skip is done is smina, but it's not grounded on the paper
         if (r_i_j >= fp_type{8}) {
            continue;
         }

         //We use this to avoid the if in the inner loop 
         const fp_type hydro_mask = static_cast<fp_type>(pair.hydrophobic_possible);
         const fp_type hbond_mask = static_cast<fp_type>(pair.hbond_possible);

         ll_score += w1 * gauss(d)+ w2 * repulsion(d) + hydro_mask * w3 * hydrophobic(d)+ hbond_mask * w4 * hbond(d);
      }

      score = pl_score + ll_score;

      return score;
   }

} // namespace mudock
