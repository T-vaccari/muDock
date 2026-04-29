#include <boost/program_options.hpp>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <format>
#include <iostream>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/vinardo_smina_helpers.hpp>
#include <mudock/compute/vinardo_affinity.hpp>
#include <mudock/compute/vinardo_score.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/pdbqt_torsion_tree.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <stdexcept>

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;

  std::filesystem::path receptor_path = std::filesystem::path{MUDOCK_SOURCE_DIR} / "data/3udd/3udd_pocket.pdbqt";
  std::filesystem::path ligand_path = std::filesystem::path{MUDOCK_SOURCE_DIR} / "data/3udd/3udd_ligand.pdbqt";

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("receptor,r",
                                      po::value(&receptor_path)->default_value(receptor_path),
                                      "Path to the receptor PDBQT file");
  arguments_description.add_options()("ligand,l",
                                      po::value(&ligand_path)->default_value(ligand_path),
                                      "Path to the ligand PDBQT file");

  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  if (vm.contains("help")) {
    std::cout << all << '\n';
    return EXIT_SUCCESS;
  }

  mudock::info("Reading and parsing protein ", receptor_path, " ...");
  auto protein = mudock::parser<mudock::dynamic_molecule>(receptor_path);

  mudock::info("Reading and parsing ligand ", ligand_path, " ...");
  auto ligand = mudock::parser<mudock::static_molecule>(ligand_path, &mudock::pdbqt_rotate_check);

  [[maybe_unused]] mudock::autodock_dynamic_layer protein_autodock{
      protein,
      [receptor_path](mudock::autodock_dynamic_layer& layer) {
        mudock::apply_autodock_forcefield_pdbqt(layer, receptor_path);
      }};
  [[maybe_unused]] mudock::autodock_static_layer ligand_autodock{
      ligand,
      [ligand_path](mudock::autodock_static_layer& layer) {
        mudock::apply_autodock_forcefield_pdbqt(layer, ligand_path);
      }};

  mudock::vinardo_layer<mudock::dynamic_containers> protein_vinardo{protein};
  mudock::vinardo_layer<mudock::static_containers> ligand_vinardo{ligand};

  const auto torsion_tree = mudock::parse_pdbqt_torsion_tree(ligand_path);
  const auto smina_mobility = mudock::build_smina_mobility_matrix(ligand, torsion_tree);
  auto preprocessed = mudock::preprocess_for_vinardo(
      protein_vinardo,
      ligand_vinardo,
      mudock::vinardo_preprocess_options{std::span<const std::uint8_t>{smina_mobility}});

  const auto breakdown = mudock::compute_vinardo_score_breakdown(protein_vinardo, ligand_vinardo, preprocessed);
  const auto score = mudock::vinardo_score(protein_vinardo, ligand_vinardo, preprocessed);
  const auto num_tors = mudock::smina_num_tors(ligand, torsion_tree, ligand_vinardo.get_vinardo_type());
  const auto affinity = mudock::vinardo_affinity(breakdown.protein_ligand, num_tors);

  if (std::abs(score - breakdown.total) > mudock::fp_type{1e-5}) {
    throw std::runtime_error("Vinardo score and breakdown total diverge");
  }

  mudock::info(std::format("Protein-ligand score: {}", breakdown.protein_ligand));
  mudock::info(std::format("Ligand-ligand score: {}", breakdown.ligand_ligand));
  mudock::info(std::format("Vinardo raw score: {}", score));
  mudock::info(std::format("smina-like num tors: {}", num_tors));
  mudock::info(std::format("Vinardo affinity: {}", affinity));

  return EXIT_SUCCESS;
}
