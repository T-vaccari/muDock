#include <boost/program_options.hpp>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
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
#include <span>
#include <stdexcept>
#include <sstream>

namespace {

mudock::fp_type compute_vinardo_affinity(const std::filesystem::path& receptor_path,
                                         const std::filesystem::path& ligand_path) {
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
      std::span<const std::uint8_t>{smina_mobility});

  const auto breakdown = mudock::compute_vinardo_score_breakdown(protein_vinardo, ligand_vinardo, preprocessed);
  const auto num_tors = mudock::smina_num_tors(ligand, torsion_tree, ligand_vinardo.get_vinardo_type());
  return mudock::vinardo_affinity(breakdown.protein_ligand, num_tors);
}

} // namespace

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;

  std::filesystem::path receptor_path;
  std::filesystem::path ligand_path;
  std::filesystem::path csv_path = std::filesystem::path{MUDOCK_SOURCE_DIR} / "test/reference/vinardo_smina_affinity.csv";
  mudock::fp_type tolerance = mudock::fp_type{1e-3};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("receptor,r",
                                      po::value(&receptor_path),
                                      "Path to the receptor PDBQT file");
  arguments_description.add_options()("ligand,l",
                                      po::value(&ligand_path),
                                      "Path to the ligand PDBQT file");
  arguments_description.add_options()("csv", po::value(&csv_path)->default_value(csv_path), "Path to reference CSV");
  arguments_description.add_options()(
      "tolerance", po::value(&tolerance)->default_value(tolerance), "Absolute tolerance");

  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  if (vm.contains("help")) {
    std::cout << all << '\n';
    return EXIT_SUCCESS;
  }

  if (!receptor_path.empty() || !ligand_path.empty()) {
    if (receptor_path.empty() || ligand_path.empty()) {
      throw std::runtime_error("Both receptor and ligand paths are required");
    }

    mudock::info(std::format("Vinardo affinity: {}", compute_vinardo_affinity(receptor_path, ligand_path)));
    return EXIT_SUCCESS;
  }

  std::ifstream input{csv_path};
  if (!input) {
    throw std::runtime_error(std::format("Cannot open reference CSV {}", csv_path.string()));
  }

  bool failed = false;
  std::string line;
  std::getline(input, line);
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }

    std::stringstream ss{line};
    std::string name;
    std::string receptor;
    std::string ligand;
    std::string smina_affinity_field;
    std::getline(ss, name, ',');
    std::getline(ss, receptor, ',');
    std::getline(ss, ligand, ',');
    std::getline(ss, smina_affinity_field, ',');

    const auto receptor_csv_path = std::filesystem::path{MUDOCK_SOURCE_DIR} / receptor;
    const auto ligand_csv_path = std::filesystem::path{MUDOCK_SOURCE_DIR} / ligand;
    const auto smina_affinity = static_cast<mudock::fp_type>(std::stod(smina_affinity_field));
    const auto mudock_affinity = compute_vinardo_affinity(receptor_csv_path, ligand_csv_path);
    const auto diff = mudock_affinity - smina_affinity;
    mudock::info(std::format("{} muDock={} smina={} diff={}",
                             name,
                             mudock_affinity,
                             smina_affinity,
                             diff));

    if (std::abs(diff) > tolerance) {
      failed = true;
      mudock::error(std::format("{} exceeds tolerance {} with abs diff {}",
                                name,
                                tolerance,
                                std::abs(diff)));
    }
  }

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
