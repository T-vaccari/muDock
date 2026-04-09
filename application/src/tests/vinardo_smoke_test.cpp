#include <cassert>
#include <filesystem>
#include <iostream>

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





void test_vinardo_type_resolution() {
  using ligand = mudock::molecule<mudock::static_containers>;
  using layer_ligand = mudock::vinardo_layer<mudock::static_containers>;

  // A carbon bonded to a heteroatom should not remain hydrophobic.
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
}

//Here we test the vinardo layer on actual parsed structures to ensure that the type resolution and property assignments work correctly 
// in a realistic scenario. 
void test_vinardo_layer_on_parsed_structures() {
  //Fetch input files
  const auto source_dir = std::filesystem::path{MUDOCK_SOURCE_DIR};
  const auto ligand_path = source_dir / "data/3udd/3udd_ligand.pdbqt";
  const auto protein_path = source_dir / "data/3udd/3udd_pocket.pdb";

  //Parse the structures
  auto parsed_ligand = mudock::parser<mudock::static_molecule>(ligand_path);
  auto parsed_protein = mudock::parser<mudock::dynamic_molecule>(protein_path);

  //Assign autodock types to the parsed structures, which is necessary for the vinardo layer to determine the correct vinardo types and properties.
  mudock::assign_autodock_types(parsed_ligand);
  mudock::assign_autodock_types(parsed_protein);

  //Create vinardo layers for both the ligand and the protein using the parsed structures.
  mudock::vinardo_layer<mudock::static_containers> ligand_layer(parsed_ligand);
  mudock::vinardo_layer<mudock::dynamic_containers> protein_layer(parsed_protein);

  //Check that the number of vinardo types assigned matches the number of atoms in both the ligand and the protein.
  const auto ligand_types = ligand_layer.get_vinardo_type();
  const auto protein_types = protein_layer.get_vinardo_type();

  assert(ligand_types.size() == static_cast<std::size_t>(parsed_ligand.num_atoms()));
  assert(protein_types.size() == static_cast<std::size_t>(parsed_protein.num_atoms()));

  //TODO : Add more specific assertion based on the specific expected types and properties

}

void test_vinardo_preprocessing() {
  const auto source_dir = std::filesystem::path{MUDOCK_SOURCE_DIR};
  const auto ligand_path = source_dir / "data/3udd/3udd_ligand.pdbqt";
  const auto protein_path = source_dir / "data/3udd/3udd_pocket.pdb";

  auto parsed_ligand = mudock::parser<mudock::static_molecule>(ligand_path);
  auto parsed_protein = mudock::parser<mudock::dynamic_molecule>(protein_path);

  mudock::assign_autodock_types(parsed_ligand);
  mudock::assign_autodock_types(parsed_protein);

  mudock::vinardo_layer<mudock::static_containers> ligand_layer(parsed_ligand);
  mudock::vinardo_layer<mudock::dynamic_containers> protein_layer(parsed_protein);

  const auto preprocessed = mudock::preprocess_for_vinardo(protein_layer, ligand_layer);

  const auto expected_pl_pairs =
      static_cast<std::size_t>(parsed_protein.num_atoms()) *
      static_cast<std::size_t>(parsed_ligand.num_atoms());

  const auto max_ll_pairs =
      static_cast<std::size_t>(parsed_ligand.num_atoms()) *
      static_cast<std::size_t>(parsed_ligand.num_atoms() - 1) / 2;

  assert(preprocessed.protein_ligand_pairs.size() == expected_pl_pairs);
  assert(preprocessed.ligand_ligand_pairs.size() < max_ll_pairs);

  for (const auto& pair : preprocessed.ligand_ligand_pairs) {
    assert(pair.ligand_atom_i_idx != pair.ligand_atom_j_idx);
  }
}

//Trivial test with a rotable bond
void test_vinardo_preprocessing_on_linear_ligand() {
  using ligand_type = mudock::molecule<mudock::static_containers>;
  using protein_type = mudock::molecule<mudock::dynamic_containers>;

  ligand_type ligand;
  ligand.resize(5, 4);

  for (int i = 0; i < 5; ++i) {
    ligand.autodock_type(i) = mudock::autodock_ff::C;
  }

  ligand.bonds(0).source = 0;
  ligand.bonds(0).dest = 1;
  ligand.bonds(0).can_rotate = false;

  ligand.bonds(1).source = 1;
  ligand.bonds(1).dest = 2;
  ligand.bonds(1).can_rotate = true;

  ligand.bonds(2).source = 2;
  ligand.bonds(2).dest = 3;
  ligand.bonds(2).can_rotate = false;

  ligand.bonds(3).source = 3;
  ligand.bonds(3).dest = 4;
  ligand.bonds(3).can_rotate = false;

  protein_type protein;
  protein.resize(0, 0);

  mudock::vinardo_layer<mudock::static_containers> ligand_layer(ligand);
  mudock::vinardo_layer<mudock::dynamic_containers> protein_layer(protein);

  const auto preprocessed = mudock::preprocess_for_vinardo(protein_layer, ligand_layer);
  const auto& ll_pairs = preprocessed.ligand_ligand_pairs;

  bool found_0_4 = false;
  bool found_0_3 = false;
  bool found_0_2 = false;

  assert(preprocessed.protein_ligand_pairs.empty());

  for (const auto& pair : ll_pairs) {
    const int i = pair.ligand_atom_i_idx;
    const int j = pair.ligand_atom_j_idx;

    if (i == 0 && j == 4) {
      found_0_4 = true;
    }
    if (i == 0 && j == 3) {
      found_0_3 = true;
    }
    if (i == 0 && j == 2) {
      found_0_2 = true;
    }

    assert(i != j);
  }

  assert(found_0_4);
  assert(!found_0_3);
  assert(!found_0_2);
}



int main() {
  test_vinardo_type_resolution();
  std::cout << "Vinardo type resolution passed\n";

  test_vinardo_layer_on_parsed_structures();
  std::cout << "Vinardo parsed-structure test passed\n";

  test_vinardo_preprocessing();
  std::cout << "Vinardo preprocessing test passed\n";

  test_vinardo_preprocessing_on_linear_ligand();
  std::cout << "Vinardo linear-ligand preprocessing test passed\n";

  return 0;
}
