#pragma once
#include <mudock/chem/molecule_layer.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/type_alias.hpp>
#include <cstdint>

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

      //TODO  :Implement constructor with a final call to prepare

      //TODO: Implement prepare that calls the helper to convert the type AD4->Vinardo

      //TODO : Implement getter

      private:
      //Here i need to declare the attributes needed 
      //ft_type it's an alias for floating point type
      atoms_array_type<vinardo_type> atom_vinardo_type;
      atoms_array_type<fp_type> atom_radius;

      //Given that it's only a flag i can use uint_8 
      atoms_array_type<std::uint8_t> atom_is_hydrophobic;
      atoms_array_type<std::uint8_t> atom_is_hbond_donor;
      atoms_array_type<std::uint8_t> atom_is_hbond_acceptor;


      
   };

}//mudock namespace
