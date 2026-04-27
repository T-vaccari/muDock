#include <filesystem>
#include <fstream>
#include <iostream>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/reader.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

//remove leading and trailing spaces from the type in the pdbqt
std::string clean(std::string_view text) {
  const auto begin = text.find_first_not_of(' ');
  if (begin == std::string_view::npos) {
    return {};
  }
  const auto end = text.find_last_not_of(' ');
  return std::string{text.substr(begin, end - begin + 1)};
}


//This should not be necessary given that are simply type names, but It's just to not break the format 
// in case of erros during the parsing
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

//Needed to retreive the original atom names given that we do not store it 
// in the molecule
std::vector<std::string> read_atom_names(const std::filesystem::path& pdbqt_path) {
  std::ifstream input{pdbqt_path};
  if (!input.good()) {
    throw std::runtime_error("Cannot open PDBQT file");
  }

  std::vector<std::string> atom_names;
  std::string line;
  while (std::getline(input, line)) {
    if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0) {
      continue;
    }
    //A bit fragile
    if (line.size() < 16) {
      throw std::runtime_error("Malformed PDBQT");
    }
    atom_names.push_back(clean(std::string_view{line}.substr(12, 4)));
  }
  return atom_names;
}

//Helper to populate the csv
template<class container_aliases>
void dump_atoms_csv(std::ofstream& output,
                    const std::string_view source_name,
                    mudock::molecule<container_aliases>& molecule,
                    const std::vector<std::string>& atom_names) {
  if (atom_names.size() != static_cast<std::size_t>(molecule.num_atoms())) {
    throw std::runtime_error("The number of atom names does not match the number of atoms in the molecule");
  }

  mudock::vinardo_layer<container_aliases> vinardo{molecule};
  const auto types = vinardo.get_vinardo_type();
  const auto radii = vinardo.get_radius();
  const auto hydrophobic = vinardo.get_is_hydrophobic();
  const auto donor = vinardo.get_is_hbond_donor();
  const auto acceptor = vinardo.get_is_hbond_acceptor();

  output << "source,id,atom_name,smina_type,xs_radius,xs_hydrophobe,xs_donor,xs_acceptor\n";
  for (int i = 0; i < molecule.num_atoms(); ++i) {
    const auto index = static_cast<std::size_t>(i);
    output << source_name << ','
           << molecule.atom_id_at(i) << ','
           << csv_escape(atom_names[index]) << ','
           << mudock::to_string(types[index]) << ','
           << radii[index] << ','
           << static_cast<int>(hydrophobic[index]) << ','
           << static_cast<int>(donor[index]) << ','
           << static_cast<int>(acceptor[index]) << '\n';
  }
}

void dump_protein(const std::filesystem::path& protein_path, const std::filesystem::path& output_path) {
  //Pars & Build the protein using the usual pipeline
  auto protein = mudock::parser<mudock::dynamic_molecule>(protein_path);
  const auto atom_names = read_atom_names(protein_path);
  //I need the apply forcefield to populate the molecule with the autodock types
  [[maybe_unused]] mudock::autodock_dynamic_layer autodock{
      protein,
      [protein_path](mudock::autodock_dynamic_layer& layer) {
        mudock::apply_autodock_forcefield_pdbqt(layer, protein_path);
      }};
  std::ofstream output{output_path};
  dump_atoms_csv(output, "protein", protein, atom_names);
}

//same here for the ligand
void dump_ligand(const std::filesystem::path& ligand_path, const std::filesystem::path& output_path) {
  auto ligand = mudock::parser<mudock::static_molecule>(ligand_path, &mudock::pdbqt_rotate_check);
  const auto atom_names = read_atom_names(ligand_path);
  [[maybe_unused]] mudock::autodock_static_layer autodock{
      ligand,
      [ligand_path](mudock::autodock_static_layer& layer) {
        mudock::apply_autodock_forcefield_pdbqt(layer, ligand_path);
      }};
  std::ofstream output{output_path};
  dump_atoms_csv(output, "ligand", ligand, atom_names);
}

} // namespace


// cmake --build --preset linux-debug --target dump_vinardo_atom_types
//To use the dumper pass the protein and the ligand and then the outdir,
//Note : Under the benchmark folder we have a folder for each test in which
// we have the *.pdbqt files for the protein and the ligand, and then the folder in which we have
//the dump. 
//Then to see the difference use the notebook
int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "How to use it: dump_vinardo_atom_types protein.pdbqt ligand.pdbqt outdir";
    return 1;
  }

  const auto protein_path = std::filesystem::path{argv[1]};
  const auto ligand_path = std::filesystem::path{argv[2]};
  const auto outdir = std::filesystem::path{argv[3]};
  std::filesystem::create_directories(outdir);

  dump_protein(protein_path, outdir / "protein_atoms.csv");
  dump_ligand(ligand_path, outdir / "ligand_atoms.csv");
  return 0;
}
