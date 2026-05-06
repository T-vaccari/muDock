#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/vinardo_smina_helpers.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/compute/vinardo_affinity.hpp>
#include <mudock/compute/vinardo_score.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/pdbqt_torsion_tree.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/type_alias.hpp>

namespace {

std::string clean(std::string_view text) {
  const auto begin = text.find_first_not_of(' ');
  if (begin == std::string_view::npos) {
    return {};
  }
  const auto end = text.find_last_not_of(' ');
  return std::string{text.substr(begin, end - begin + 1)};
}

std::string csv_escape(const std::string& value) {
  if (value.find_first_of(",\"") == std::string::npos) {
    return value;
  }

  std::string escaped = "\"";
  for (const char c : value) {
    if (c == '"') {
      escaped += "\"\"";
    } else {
      escaped += c;
    }
  }
  escaped += "\"";
  return escaped;
}

struct pdbqt_atom_metadata {
  int id;
  std::string name;
};

std::vector<pdbqt_atom_metadata> read_atom_metadata(const std::filesystem::path& pdbqt_path) {
  std::ifstream input{pdbqt_path};
  if (!input.good()) {
    throw std::runtime_error("Cannot open PDBQT file");
  }

  std::vector<pdbqt_atom_metadata> atoms;
  std::string line;
  while (std::getline(input, line)) {
    if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0) {
      continue;
    }
    atoms.push_back(pdbqt_atom_metadata{std::stoi(line.substr(6, 5)), clean(std::string_view{line}.substr(12, 4))});
  }
  return atoms;
}

struct scoring_terms {
  mudock::fp_type distance;
  mudock::fp_type r2;
  mudock::vinardo_pair_terms pair;
};

scoring_terms compute_terms(const mudock::fp_type dx,
                            const mudock::fp_type dy,
                            const mudock::fp_type dz,
                            const mudock::fp_type radius_sum,
                            const bool hydrophobic_possible,
                            const bool hbond_possible) {
  scoring_terms terms{};
  terms.r2 = dx * dx + dy * dy + dz * dz;
  terms.distance = std::sqrt(terms.r2);
  const auto surface_distance = terms.distance - radius_sum;
  terms.pair = mudock::compute_vinardo_pair_terms(surface_distance, hydrophobic_possible, hbond_possible);
  return terms;
}

void write_header(std::ofstream& output) {
  output << "kind,ligand_atom_id,other_atom_id,ligand_atom_index,other_atom_index,other_role,"
            "ligand_atom_name,other_atom_name,ligand_atom_type,other_atom_type,distance,r2,"
            "hydrophobic_possible,hbond_possible,radius_sum,gauss,repulsion,hydrophobic,"
            "non_dir_h_bond,weighted_gauss,weighted_repulsion,weighted_hydrophobic,"
            "weighted_non_dir_h_bond,weighted_energy\n";
}

void write_row(std::ofstream& output,
               const char* kind,
               const int ligand_atom_id,
               const int other_atom_id,
               const int ligand_atom_index,
               const int other_atom_index,
               const char* other_role,
               const std::string& ligand_atom_name,
               const std::string& other_atom_name,
               const mudock::vinardo_atom_type ligand_atom_type,
               const mudock::vinardo_atom_type other_atom_type,
               const bool hydrophobic_possible,
               const bool hbond_possible,
               const mudock::fp_type radius_sum,
               const scoring_terms& terms) {
  output << kind << ','
         << ligand_atom_id << ','
         << other_atom_id << ','
         << ligand_atom_index << ','
         << other_atom_index << ','
         << other_role << ','
         << csv_escape(ligand_atom_name) << ','
         << csv_escape(other_atom_name) << ','
         << mudock::to_string(ligand_atom_type) << ','
         << mudock::to_string(other_atom_type) << ','
         << terms.distance << ','
         << terms.r2 << ','
         << static_cast<int>(hydrophobic_possible) << ','
         << static_cast<int>(hbond_possible) << ','
         << radius_sum << ','
         << terms.pair.gauss << ','
         << terms.pair.repulsion << ','
         << terms.pair.hydrophobic << ','
         << terms.pair.hbond << ','
         << terms.pair.weighted_gauss << ','
         << terms.pair.weighted_repulsion << ','
         << terms.pair.weighted_hydrophobic << ','
         << terms.pair.weighted_hbond << ','
         << terms.pair.weighted_energy << '\n';
}

void write_score_summary(const std::filesystem::path& pairs_output_path,
                         const unsigned smina_like_tors,
                         const mudock::fp_type raw_c_inter,
                         const mudock::fp_type raw_c_intra,
                         const mudock::vinardo_score_breakdown& breakdown,
                         const mudock::fp_type affinity) {
  const auto summary_path = pairs_output_path.parent_path() / "score_summary.csv";
  std::ofstream output{summary_path};
  output << std::setprecision(10);
  output << "smina_like_num_tors,raw_C_inter,raw_C_intra,raw_C_inter_plus_intra,"
            "production_C_inter,production_C_intra,production_C_inter_plus_intra,affinity\n";
  output << smina_like_tors << ','
         << raw_c_inter << ','
         << raw_c_intra << ','
         << raw_c_inter + raw_c_intra << ','
         << breakdown.protein_ligand << ','
         << breakdown.ligand_ligand << ','
         << breakdown.total << ','
         << affinity << '\n';
}

void dump_scoring_pairs(const std::filesystem::path& protein_path,
                        const std::filesystem::path& ligand_path,
                        const std::filesystem::path& output_path) {
  auto protein = mudock::parser<mudock::dynamic_molecule>(protein_path);
  auto ligand = mudock::parser<mudock::static_molecule>(ligand_path, &mudock::pdbqt_rotate_check);
  const auto protein_atoms = read_atom_metadata(protein_path);
  const auto ligand_atoms = read_atom_metadata(ligand_path);

  [[maybe_unused]] mudock::autodock_dynamic_layer protein_autodock{
      protein,
      [protein_path](mudock::autodock_dynamic_layer& layer) {
        mudock::apply_autodock_forcefield_pdbqt(layer, protein_path);
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
  const auto protein_types = protein_vinardo.get_vinardo_type();
  const auto ligand_types = ligand_vinardo.get_vinardo_type();
  const auto breakdown = mudock::compute_vinardo_score_breakdown(protein_vinardo, ligand_vinardo, preprocessed);
  const auto smina_like_tors = mudock::smina_num_tors(ligand, torsion_tree, ligand_types);
  const auto affinity = mudock::vinardo_affinity(breakdown.protein_ligand, smina_like_tors);
  mudock::fp_type raw_c_inter = mudock::fp_type{0};
  mudock::fp_type raw_c_intra = mudock::fp_type{0};

  std::ofstream output{output_path};
  output << std::setprecision(10);
  write_header(output);

  for (const auto& pair : preprocessed.protein_ligand_pairs) {
    const auto protein_idx = static_cast<std::size_t>(pair.protein_atom_idx);
    const auto ligand_idx = static_cast<std::size_t>(pair.ligand_atom_idx);
    const auto terms = compute_terms(protein.x(pair.protein_atom_idx) - ligand.x(pair.ligand_atom_idx),
                                     protein.y(pair.protein_atom_idx) - ligand.y(pair.ligand_atom_idx),
                                     protein.z(pair.protein_atom_idx) - ligand.z(pair.ligand_atom_idx),
                                     pair.radius_sum,
                                     pair.hydrophobic_possible,
                                     pair.hbond_possible);
    raw_c_inter += terms.pair.weighted_energy;
    if (terms.distance >= mudock::fp_type{8}) {
      continue;
    }
    write_row(output,
              "protein_ligand",
              ligand_atoms.at(ligand_idx).id,
              protein_atoms.at(protein_idx).id,
              pair.ligand_atom_idx,
              pair.protein_atom_idx,
              "protein",
              ligand_atoms.at(ligand_idx).name,
              protein_atoms.at(protein_idx).name,
              ligand_types[ligand_idx],
              protein_types[protein_idx],
              pair.hydrophobic_possible,
              pair.hbond_possible,
              pair.radius_sum,
              terms);
  }

  for (const auto& pair : preprocessed.ligand_ligand_pairs) {
    const auto i = static_cast<std::size_t>(pair.ligand_atom_i_idx);
    const auto j = static_cast<std::size_t>(pair.ligand_atom_j_idx);
    const auto terms = compute_terms(ligand.x(pair.ligand_atom_i_idx) - ligand.x(pair.ligand_atom_j_idx),
                                     ligand.y(pair.ligand_atom_i_idx) - ligand.y(pair.ligand_atom_j_idx),
                                     ligand.z(pair.ligand_atom_i_idx) - ligand.z(pair.ligand_atom_j_idx),
                                     pair.radius_sum,
                                     pair.hydrophobic_possible,
                                     pair.hbond_possible);
    raw_c_intra += terms.pair.weighted_energy;
    if (terms.distance >= mudock::fp_type{8}) {
      continue;
    }
    write_row(output,
              "ligand_ligand",
              ligand_atoms.at(i).id,
              ligand_atoms.at(j).id,
              pair.ligand_atom_i_idx,
              pair.ligand_atom_j_idx,
              "ligand",
              ligand_atoms.at(i).name,
              ligand_atoms.at(j).name,
              ligand_types[i],
              ligand_types[j],
              pair.hydrophobic_possible,
              pair.hbond_possible,
              pair.radius_sum,
              terms);
  }

  write_score_summary(output_path, smina_like_tors, raw_c_inter, raw_c_intra, breakdown, affinity);
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "usage: dump_scoring_pairs protein.pdbqt ligand.pdbqt out.csv\n";
    return 1;
  }

  const auto protein_path = std::filesystem::path{argv[1]};
  const auto ligand_path = std::filesystem::path{argv[2]};
  const auto output_path = std::filesystem::path{argv[3]};
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
  }

  dump_scoring_pairs(protein_path, ligand_path, output_path);
  return 0;
}
