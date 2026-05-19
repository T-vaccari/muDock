#pragma once

#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/molecule.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/vinardo_smina_helpers.hpp>


namespace mudock {

struct vinardo_ligand : public vinardo_layer<static_containers> {
    vinardo_ligand(static_molecule& ligand): vinardo_layer<static_containers>(ligand) {
      prepare();
   }

   [[nodiscard]] const auto& get_ligand_ligand_pairs() const {
      return ll_pairs;
   }

   [[nodiscard]] unsigned get_num_tors() const {
      return num_tors;
   }

   private:
      //All the data are in the base molecule
      std::vector<vinardo_ligand_ligand_pair> ll_pairs;
      unsigned int num_tors = 0; //Number of torsional degrees of freedom, used for score normalization(e.g. affinity)
      

      void prepare() {
         auto& ligand = this->get_base_molecule();

         ll_pairs = preprocess_ligand_vinardo(*this, ligand.pdbqt_ligand_data.mobility_matrix);
         //Here I am getting a warning on this->get_vinardo_type(), maybe related to containers? Implicit conversion changes signedness: 'const int' to 'size_type' (aka 'unsigned long')
         num_tors = smina_num_tors(ligand, ligand.pdbqt_ligand_data.rotors,this->get_vinardo_type());
      }
      
};

}//namespace mudock
