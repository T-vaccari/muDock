#include <mudock/molecule/containers.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <cassert>

// Used to compile the new files: clang++ -std=c++20 -I mudock/include -c smoketest.cpp
// TODO : Add to the cmake this target
// Simple smoke test to compile the new files.
int main() { 
   //Instantiate a molecule
   using ligand = mudock::molecule<mudock::static_containers>;
   using layer_ligand = mudock::vinardo_layer<mudock::static_containers>;
   using protein = mudock::molecule<mudock::dynamic_containers>;
   using layer_protein = mudock::vinardo_layer<mudock::dynamic_containers>;

   
   //Test for the types with ambiguity that needs to be resolved using the two step approach

   // A carbon bonded to a heteroatom should not remain hydrophobic,
   // the second pass move the carbon from the hydrophobic aliphatic-carbon to the non-hydrophobic one.
   // Expected: AliphaticCarbonXSNonHydrophobe, hydrophobic == false.
   {
      ligand mol;
      mol.resize(2, 1);

      mol.autodock_type(0) = mudock::autodock_ff::C;
      mol.autodock_type(1) = mudock::autodock_ff::OA;

      mol.bonds(0).source = 0;
      mol.bonds(0).dest = 1;

      layer_ligand layer(mol);
      const auto vinardo_types = layer.get_vinardo_type();

      assert(vinardo_types[0] == mudock::vinardo_atom_type::AliphaticCarbonXSNonHydrophobe);
      assert(vinardo_types[1] == mudock::vinardo_atom_type::OxygenXSAcceptor);

      const auto hydrophobic = layer.get_is_hydrophobic();
      assert(hydrophobic[0] == false);

   }

   // Expected: NitrogenXSDonor, donor == true, acceptor == false.
   {
      ligand mol;
      mol.resize(2, 1);

      mol.autodock_type(0) = mudock::autodock_ff::N;
      mol.autodock_type(1) = mudock::autodock_ff::HD;

      mol.bonds(0).source = 0;
      mol.bonds(0).dest = 1;

      layer_ligand layer(mol);
      const auto vinardo_types = layer.get_vinardo_type();

      assert(vinardo_types[0] == mudock::vinardo_atom_type::NitrogenXSDonor);
      assert(vinardo_types[1] == mudock::vinardo_atom_type::PolarHydrogen);

      const auto donor = layer.get_is_hbond_donor();
      const auto acceptor = layer.get_is_hbond_acceptor();

      assert(donor[0] == true);
      assert(acceptor[0] == false);

   }

   // Expected: OxygenXSDonorAcceptor, donor == true, acceptor == true.
   {
      ligand mol;
      mol.resize(2, 1);

      mol.autodock_type(0) = mudock::autodock_ff::OA;
      mol.autodock_type(1) = mudock::autodock_ff::HD;

      mol.bonds(0).source = 0;
      mol.bonds(0).dest = 1;

      layer_ligand layer(mol);
      const auto vinardo_types = layer.get_vinardo_type();

      assert(vinardo_types[0] == mudock::vinardo_atom_type::OxygenXSDonorAcceptor);
      assert(vinardo_types[1] == mudock::vinardo_atom_type::PolarHydrogen);

      const auto donor = layer.get_is_hbond_donor();
      const auto acceptor = layer.get_is_hbond_acceptor();

      assert(donor[0] == true);
      assert(acceptor[0] == true);

   }


   return 0;
}
