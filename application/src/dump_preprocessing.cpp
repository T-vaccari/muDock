#include <filesystem>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/chem/vinardo_preprocessing.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/reader.hpp>
#include <stdexcept>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

//Struct to represent the tree of branches in the pdbqt, I need it to apply the smina mobility semantics to the ligand-ligand pairs
struct pdbqt_branch_node;

struct pdbqt_branch {
  std::uint32_t from;
  std::uint32_t to;
  pdbqt_branch_node* child;
};

struct pdbqt_branch_node {
  std::vector<std::uint32_t> atom_ids;
  std::vector<pdbqt_branch> children;
};


//as usual remove leading and trailing spaces from the type in the pdbqt
std::string clean(std::string_view text) {
  const auto begin = text.find_first_not_of(' ');
  if (begin == std::string_view::npos) {
    return {};
  }
  const auto end = text.find_last_not_of(' ');
  return std::string{text.substr(begin, end - begin + 1)};
}

// just to not break the format in case of erros during the parsing
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
    if (line.size() < 16) {
      throw std::runtime_error("Malformed PDBQT");
    }
    atom_names.push_back(clean(std::string_view{line}.substr(12, 4)));
  }
  return atom_names;
}

std::uint32_t atom_id_from_pdbqt_line(const std::string& line) {
  return static_cast<std::uint32_t>(std::stoul(line.substr(6, 5)));
}

// Recursively parse branches and build the representation of the torsion tree.
pdbqt_branch parse_branch_tree(const std::vector<std::string>& lines, std::size_t& pos, const std::uint32_t from, const std::uint32_t to,
                              std::vector<pdbqt_branch_node>& storage) {

  storage.emplace_back(); //Allocate space for the new node
  auto& node = storage.back();

  while (pos < lines.size()) {
    const auto& line = lines[pos];
    pos++;
    if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) { //(line.starts_with("ATOM") || line.starts_with("HETATM"))
      node.atom_ids.push_back(atom_id_from_pdbqt_line(line));
    } else if (line.rfind("BRANCH", 0) == 0) {
      //We have found a new recursive branch , we need to parse:
      // Branch from(parent) to(child)
      std::istringstream stream{line};
      std::string token;
      std::uint32_t child_from;
      std::uint32_t child_to;
      stream >> token >> child_from >> child_to;
      auto child = parse_branch_tree(lines, pos, child_from, child_to, storage);
      node.atom_ids.push_back(child.to);
      node.children.push_back(child);
    } else if (line.rfind("ENDBRANCH", 0) == 0) {
      //We have found the closing keyword of the current branch we can return to the parent
      break;
    }
  }
  // The child-side axis atom is kept in the branch metadata, not as a normal child-node atom.
  node.atom_ids.erase(std::remove(node.atom_ids.begin(), node.atom_ids.end(), to), node.atom_ids.end());
  return pdbqt_branch{from, to, &node}; //return the branch with the pointer to the node that contains the children and the atoms in the branch
}


//Same parsing as in smina
pdbqt_branch_node parse_pdbqt_tree(const std::filesystem::path& pdbqt_path, std::vector<pdbqt_branch_node>& storage) {
  //This is the top level function in which we build the three
  std::ifstream input{pdbqt_path};
  if (!input.good()) {
    throw std::runtime_error("Cannot open PDBQT file");
  }

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  storage.reserve(lines.size());

  pdbqt_branch_node root;
  std::size_t pos = 0;
  while (pos < lines.size() && lines[pos].rfind("ROOT", 0) != 0) {
    ++pos;
  }
  if (pos == lines.size()) {
    throw std::runtime_error("Missing ROOT in PDBQT file");
  }
  ++pos;
  //Here we parse the root block
  while (pos < lines.size()) {
    const auto& current = lines[pos];
    pos++;
    if (current.rfind("ATOM", 0) == 0 || current.rfind("HETATM", 0) == 0) {
      root.atom_ids.push_back(atom_id_from_pdbqt_line(current));
    } else if (current.rfind("ENDROOT", 0) == 0) {
      break;
    }
  }

  //Now we need to parse the branch section
  while (pos < lines.size()) {
    const auto& current = lines[pos];
    pos++;
    if (current.rfind("BRANCH", 0) != 0) {
      continue;
    }
    std::istringstream stream{current};
    std::string token;
    std::uint32_t from;
    std::uint32_t to;
    stream >> token >> from >> to;
    auto child = parse_branch_tree(lines, pos, from, to, storage);
    //The invariant here is that we maintain the "to" atom in the parent,
    //from and to are maintained also in the struct as they are needed to apply the mobility semantics
    root.atom_ids.push_back(child.to);
    root.children.push_back(child);
  }

  return root;
}

void mark_not_variable(std::vector<std::uint8_t>& variable, const std::size_t num_atoms, const int i, const int j) {
  // Variable means the distance can vary so the two atoms are relatovely movable
  variable[i * num_atoms + j] = 0;
  variable[j * num_atoms + i] = 0;
}

void apply_smina_fixed_marks(const pdbqt_branch_node& node, const std::unordered_map<std::uint32_t, int>& atom_id_to_index, std::vector<std::uint8_t>& variable,
                            const std::size_t num_atoms) {
  //mark as non movable all the atoms inside the same branch
  //but i need to convert the atom id into the relative position in the matrix
  for (std::size_t i = 0; i < node.atom_ids.size(); ++i) {
    for (std::size_t j = i + 1; j < node.atom_ids.size(); ++j) {
      mark_not_variable(variable, num_atoms, atom_id_to_index.at(node.atom_ids[i]),atom_id_to_index.at(node.atom_ids[j]));
    }
  }
  //Iterate over the child branches
  for (const auto& child : node.children) {
    const auto from_idx = atom_id_to_index.at(child.from);
    const auto to_idx = atom_id_to_index.at(child.to);
    mark_not_variable(variable, num_atoms, from_idx, to_idx);
    for (const auto atom_id : child.child->atom_ids) {
      //!!! Note here that with the smina semantics we are excluding as relativley movable
      //couples of atoms that relate to from&to, instead with the fragments
      //this were included
      const auto atom_idx = atom_id_to_index.at(atom_id);
      mark_not_variable(variable, num_atoms, from_idx, atom_idx);
      mark_not_variable(variable, num_atoms, to_idx, atom_idx);
    }
    apply_smina_fixed_marks(*child.child, atom_id_to_index, variable, num_atoms);
  }
}

std::vector<std::uint8_t> pdbqt_smina_variable_matrix(const mudock::static_molecule& ligand, const std::filesystem::path& ligand_path) {
  const auto num_atoms = static_cast<std::size_t>(ligand.num_atoms());
  std::vector<std::uint8_t> variable(num_atoms * num_atoms, 1);
  //Adjust the diagonal
  for (std::size_t i = 0; i < num_atoms; ++i) {
    variable[i * num_atoms + i] = 0;
  }

  std::unordered_map<std::uint32_t, int> atom_id_to_index;
  for (int i = 0; i < ligand.num_atoms(); ++i) {
    atom_id_to_index.emplace(ligand.atom_id_at(i), i);
  }

  std::vector<pdbqt_branch_node> storage;
  auto root = parse_pdbqt_tree(ligand_path, storage);
  apply_smina_fixed_marks(root, atom_id_to_index, variable, num_atoms);
  return variable;
}

void dump_preprocessing(const std::filesystem::path& protein_path, const std::filesystem::path& ligand_path, const std::filesystem::path& output_path, const bool use_smina_mobility){
  //Parse & Build the protein and the ligand using the usual pipeline
  auto protein = mudock::parser<mudock::dynamic_molecule>(protein_path);
  auto ligand = mudock::parser<mudock::static_molecule>(ligand_path, &mudock::pdbqt_rotate_check);
  const auto protein_names = read_atom_names(protein_path);
  const auto ligand_names = read_atom_names(ligand_path);

  //Apply the forcefield so I have the types in molecule
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

  //Build the vinardo protein and ligand
  mudock::vinardo_layer<mudock::dynamic_containers> protein_vinardo{protein};
  mudock::vinardo_layer<mudock::static_containers> ligand_vinardo{ligand};
  const auto preprocessed = mudock::preprocess_for_vinardo(protein_vinardo, ligand_vinardo);
  const auto protein_types = protein_vinardo.get_vinardo_type();
  const auto ligand_types = ligand_vinardo.get_vinardo_type();
  const auto smina_variable = use_smina_mobility
                                  ? pdbqt_smina_variable_matrix(ligand, ligand_path)
                                  : std::vector<std::uint8_t>{};
  std::ofstream output{output_path};
  output << "kind,ligand_atom_id,other_atom_id,ligand_atom_name,other_atom_name,"
            "ligand_atom_type,other_atom_type,other_role,hydrophobic_possible,hbond_possible,radius_sum\n";

  for (const auto& pair : preprocessed.protein_ligand_pairs) {
    const auto protein_idx = static_cast<std::size_t>(pair.protein_atom_idx);
    const auto ligand_idx = static_cast<std::size_t>(pair.ligand_atom_idx);
    output << "protein_ligand,"
           << ligand.atom_id_at(pair.ligand_atom_idx) << ','
           << protein.atom_id_at(pair.protein_atom_idx) << ','
           << csv_escape(ligand_names[ligand_idx]) << ','
           << csv_escape(protein_names[protein_idx]) << ','
           << mudock::to_string(ligand_types[ligand_idx]) << ','
           << mudock::to_string(protein_types[protein_idx]) << ','
           << "protein,"
           << static_cast<int>(pair.hydrophobic_possible) << ','
           << static_cast<int>(pair.hbond_possible) << ','
           << pair.radius_sum << '\n';
  }

  for (const auto& pair : preprocessed.ligand_ligand_pairs) {
    const auto i = static_cast<std::size_t>(pair.ligand_atom_i_idx);
    const auto j = static_cast<std::size_t>(pair.ligand_atom_j_idx);
    //Here we skip the couples accordingly to the smina semantic, if requested
    if (use_smina_mobility && !smina_variable[i * static_cast<std::size_t>(ligand.num_atoms()) + j]) {
      continue;
    }
    output << "ligand_ligand,"
           << ligand.atom_id_at(pair.ligand_atom_i_idx) << ','
           << ligand.atom_id_at(pair.ligand_atom_j_idx) << ','
           << csv_escape(ligand_names[i]) << ','
           << csv_escape(ligand_names[j]) << ','
           << mudock::to_string(ligand_types[i]) << ','
           << mudock::to_string(ligand_types[j]) << ','
           << "ligand,"
           << static_cast<int>(pair.hydrophobic_possible) << ','
           << static_cast<int>(pair.hbond_possible) << ','
           << pair.radius_sum << '\n';
  }
}

}  // namespace


//To build
// cmake --preset linux-debug
//cmake --build --preset linux-debug --target dump_preprocessing

// to run
// ./build/linux-debug/application/dump_preprocessing protein.pdbqt ligand.pdbqt out.csv


int main(int argc, char* argv[]) {
  if (argc != 4 && argc != 5) {
    std::cerr << "usage: dump_preprocessing protein.pdbqt ligand.pdbqt out.csv "
                 " to use the smina semantic of rotable bonds:[--pdbqt-smina-mobility]";
    return 1;
  }
  const auto protein_path = std::filesystem::path{argv[1]};
  const auto ligand_path = std::filesystem::path{argv[2]};
  const auto output_path = std::filesystem::path{argv[3]};
  const bool use_smina_mobility = argc == 5 && std::string_view{argv[4]} == "--pdbqt-smina-mobility";
  if (argc == 5 && !use_smina_mobility) {
    std::cerr << "unknown option: " << argv[4];
    return 1;
  }
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
  }

  dump_preprocessing(protein_path, ligand_path, output_path, use_smina_mobility);
  return 0;
}
