#pragma once
#include <mudock/chem/molecule_layer.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/type_alias.hpp>
#include <cstdint>
#include <mudock/chem/assign_vinardo_type.hpp>

/*
This is a layer built on top of molecule_layer. It stores the atom-level static properties
needed to build a Vinardo view of a molecule and to support later preprocessing and scoring.
It does not store pairwise or pose-dependent data.
*/

namespace mudock{
   //As usual we have defined the data layout in the molecule, but the policy on how data are store in the backend it's policy depedent.
   //We have two main containers : 1)Static Container(Tipically used for ligand), in which we use array so we use fixed size containers,
   //the number of atoms for the ligand 2)Dinamic container(tipically used for protein)
   template<class container_aliases>
   //Concept used to check if the container alias is part of the one already defined
   requires is_container_specification<container_aliases>



   //Creating the new layer with inheritance from the base molecule layer
   struct vinardo_layer: public molecule_layer<container_aliases>{

      //Using the policy provided by the container to define an memory layout custom for the defined type(the memory shape with container aliases)
      //And what is inside with the T
      template<typename T>
      using atoms_array_type = container_aliases::template atoms_size<T>;
      using molecule_type = molecule<container_aliases>; //Simplified the firm of the method

      //Constructor (firstly call super() then initialize internal attributes)
      vinardo_layer(molecule_type& molecule): molecule_layer<container_aliases>(molecule) {
         const auto num_atoms = molecule.num_atoms();
         mudock::resize(atom_vinardo_type, num_atoms);
         mudock::resize(atom_radius, num_atoms);
         mudock::resize(atom_is_hydrophobic, num_atoms);
         mudock::resize(atom_is_hbond_donor, num_atoms);
         mudock::resize(atom_is_hbond_acceptor, num_atoms);
         prepare(molecule);
      }

      private: 
      void prepare(molecule_type& molecule){
         // Assign the Vinardo type and populate properties
         assign_vinardo_type(
            molecule,
            this->atom_vinardo_type,
            this->atom_radius,
            this->atom_is_hydrophobic,
            this->atom_is_hbond_donor,
            this->atom_is_hbond_acceptor
         );
      }
      
      public:
      //Minimal getters
      [[nodiscard]] inline auto get_vinardo_type() const {
         return make_span(atom_vinardo_type, this->get_base_molecule().num_atoms());
      }
      [[nodiscard]] inline auto get_is_hydrophobic() const {
         return make_span(atom_is_hydrophobic, this->get_base_molecule().num_atoms());
      }

      [[nodiscard]] inline auto get_is_hbond_donor() const {
         return make_span(atom_is_hbond_donor, this->get_base_molecule().num_atoms());
      }

      [[nodiscard]] inline auto get_is_hbond_acceptor() const {
         return make_span(atom_is_hbond_acceptor, this->get_base_molecule().num_atoms());
      }
      
      [[nodiscard]] inline auto get_radius() const {
         return make_span(atom_radius, this->get_base_molecule().num_atoms());
      }

      [[nodiscard]] inline auto num_atoms() const { 
         return this->get_base_molecule().num_atoms(); 
      }




      private:
      //Here i need to declare the attributes needed 
      //ft_type it's an alias for floating point type
      //TODO : Maybe in molecule?
      atoms_array_type<vinardo_atom_type> atom_vinardo_type;
      atoms_array_type<fp_type> atom_radius;
      
      atoms_array_type<std::uint8_t> atom_is_hydrophobic;
      atoms_array_type<std::uint8_t> atom_is_hbond_donor;
      atoms_array_type<std::uint8_t> atom_is_hbond_acceptor;

   };

}//mudock namespace
