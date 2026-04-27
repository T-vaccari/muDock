#pragma once  

#include <mudock/molecule/fragments.hpp>
#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/type_alias.hpp>
#include <queue>
#include <vector>
#include <cstdint>
#include <utility>



namespace mudock{
//TODO : Insert utility 
struct vinardo_protein_ligand_pair {
   int protein_atom_idx;
   int ligand_atom_idx;
   std::uint8_t hydrophobic_possible;
   std::uint8_t hbond_possible;
   mudock::fp_type radius_sum;
};

struct vinardo_ligand_ligand_pair {
   int ligand_atom_i_idx;
   int ligand_atom_j_idx;
   std::uint8_t hydrophobic_possible;
   std::uint8_t hbond_possible;
   mudock::fp_type radius_sum;
};

struct vinardo_preprocessed_pairs {
   std::vector<vinardo_protein_ligand_pair> protein_ligand_pairs;
   std::vector<vinardo_ligand_ligand_pair> ligand_ligand_pairs;
};

struct vinardo_preprocess_options {
   bool include_hydrogens = false;
};

//Template functions to append pairs to the output vector, with the option to include or exclude hydrogens.
// The template approach is used because it allows to avoid
// the if in the inner loop and can be directly resolved at compile time

template<bool include_hydrogens>
inline void append_protein_ligand_pairs(std::vector<vinardo_protein_ligand_pair>& pl_pairs,
                                        const vinardo_layer<mudock::dynamic_containers>& protein_layer,
                                        const vinardo_layer<mudock::static_containers>& ligand_layer,
                                        const auto& protein_radii,
                                        const auto& ligand_radii,
                                        const auto& protein_donor,
                                        const auto& protein_acceptor,
                                        const auto& protein_hydro,
                                        const auto& ligand_donor,
                                        const auto& ligand_acceptor,
                                        const auto& ligand_hydro) {
   const auto protein_atoms = static_cast<std::size_t>(protein_layer.num_atoms());
   const auto ligand_atoms = static_cast<std::size_t>(ligand_layer.num_atoms());

   pl_pairs.reserve(protein_atoms * ligand_atoms);
   for (std::size_t i = 0; i < protein_atoms; ++i) {
      for (std::size_t j = 0; j < ligand_atoms; ++j) {
         if constexpr (!include_hydrogens) {
            //Only hydrogen and polar hydrogen have radius 0 so I can check this property
            //Maybe comparing zeros in floating point it's a bit fragile, so I can insert a tolerance if needed.
            if (protein_radii[i] == mudock::fp_type{0} || ligand_radii[j] == mudock::fp_type{0}) {
               continue;
            }
         }
         vinardo_protein_ligand_pair pair;
         pair.protein_atom_idx      = static_cast<int>(i);
         pair.ligand_atom_idx       = static_cast<int>(j);
         pair.hydrophobic_possible  = protein_hydro[i] && ligand_hydro[j];
         pair.hbond_possible        = (protein_donor[i] && ligand_acceptor[j]) ||
                               (protein_acceptor[i] && ligand_donor[j]);
         pair.radius_sum            = protein_radii[i] + ligand_radii[j];
         pl_pairs.push_back(pair);
      }
   }
}

template<bool include_hydrogens>
inline void append_ligand_ligand_pairs(std::vector<vinardo_ligand_ligand_pair>& ll_pairs,
                                       const vinardo_layer<mudock::static_containers>& ligand_layer,
                                       const auto& ligand_radii,
                                       const auto& ligand_donor,
                                       const auto& ligand_acceptor,
                                       const auto& ligand_hydro,
                                       const std::vector<std::uint8_t>& relatively_movable,
                                       const std::vector<std::uint8_t>& within_three_bonds) {
   const auto num_atoms = static_cast<std::size_t>(ligand_layer.num_atoms());

   ll_pairs.reserve((num_atoms * (num_atoms - 1)) / 2);
   for (std::size_t i = 0; i < num_atoms; ++i) {
      for (std::size_t j = i + 1; j < num_atoms; ++j) {
         if constexpr (!include_hydrogens) {
            if (ligand_radii[i] == mudock::fp_type{0} || ligand_radii[j] == mudock::fp_type{0}) {
               continue;
            }
         }
         if (relatively_movable[i * num_atoms + j] && !(within_three_bonds[i * num_atoms + j])) {
            vinardo_ligand_ligand_pair pair;
            pair.ligand_atom_i_idx    = static_cast<int>(i);
            pair.ligand_atom_j_idx    = static_cast<int>(j);
            pair.hydrophobic_possible = ligand_hydro[i] && ligand_hydro[j];
            pair.hbond_possible       = (ligand_donor[i] && ligand_acceptor[j]) ||
                                  (ligand_acceptor[i] && ligand_donor[j]);
            pair.radius_sum           = ligand_radii[i] + ligand_radii[j];
            ll_pairs.push_back(pair);
         }
      }
   }
}

// To compute the relatively movable pairs, for each atom pair (i, j)
// check whether there exists at least one rotatable bond whose fragment mask
// places i and j on opposite sides.
// If such a bond exists, the pair is relatively movable.
// Otherwise, it is not relatively movable.
// Note: ligand_fragments contains all the ligand cuts(partitions) induced by removing,
// one at a time, each rotatable bond from the molecular graph.
inline std::vector<std::uint8_t> precompute_rotatable_pairs(const mudock::fragments<mudock::static_containers>& ligand_fragments, const std::size_t num_atoms){
   std::vector<std::uint8_t> relatively_movable(num_atoms * num_atoms, 0);
   for (std::size_t i = 0; i < num_atoms; ++i) {
      for (std::size_t j = i + 1; j < num_atoms; ++j) {
         bool movable = false;
         for (int rot_idx = 0; rot_idx < ligand_fragments.get_num_rotatable_bonds(); ++rot_idx) {
            auto mask = ligand_fragments.get_mask(rot_idx);

            if (!mudock::same_fragment(mask, i, j)) {
               movable = true;
               break;
            }
         }

         if (movable) {
            relatively_movable[i * num_atoms + j] = 1;
            relatively_movable[j * num_atoms + i] = 1;
         }
      }
   }
   return relatively_movable;
}

inline std::vector<std::uint8_t> precompute_within_n_bonds(const auto& g, const std::size_t num_atoms, const std::size_t n){
   std::vector<std::uint8_t> within_n_bonds(num_atoms * num_atoms, 0);

   for(size_t i =0 ; i < num_atoms; ++i){
      //For each atom i, we perform a breadth first search up to depth n to find all atoms that are within n bonds from i.
      std::queue<mudock::vertex_type> current_queue;
      std::queue<mudock::vertex_type> next_queue;

      current_queue.push(static_cast<mudock::vertex_type>(i));
      std::vector<std::uint8_t> visited(num_atoms,0);
      visited[i] = 1;
      // Note: this BFS assumes atom indices match graph vertex descriptors.
      // This is currently true because make_graph adds one vertex per atom in index order, but imo it's fragile.
      for(size_t depth = 0; depth < n; ++depth){
         while(!current_queue.empty()){
            auto visiting_node = current_queue.front();  //This return the value
            current_queue.pop();                                       //This remove the value
            const auto [begin, end] = boost::adjacent_vertices(visiting_node, g);

            for(auto iterator = begin; iterator != end; ++iterator){
               //The iterator points towards a vertex descriptr
               auto neighbor = *iterator;
               auto neighbor_idx = g[neighbor].atom_index;
               if(visited[neighbor_idx] != 1){
                  next_queue.push(neighbor);
                  visited[neighbor_idx] = 1;
               }

            }

         }
         //Swap the queue 
         std::swap(current_queue, next_queue);
      }
      
      // visited[j] == 1 means that atom j is reachable from source atom i
      // within at most n BFS levels, i.e. within n covalent bonds.

      for(size_t j = i+1 ; j< num_atoms; ++j){
         if(visited[j]==1){
            within_n_bonds[i*num_atoms + j] = 1;
            within_n_bonds[j*num_atoms + i] = 1;
         }
      }

   }
   return within_n_bonds;
}

inline vinardo_preprocessed_pairs preprocess_for_vinardo(
   vinardo_layer<mudock::dynamic_containers>& protein_layer,
   vinardo_layer<mudock::static_containers>& ligand_layer,
   const vinardo_preprocess_options options = {}){
   

   //Here I need to build the list of pairs of atoms that will be used for the vinardo scoring. 
   // I need to consider both protein-ligand pairs and ligand-ligand pairs.

   const auto protein_radii = protein_layer.get_radius();
   const auto ligand_radii  = ligand_layer.get_radius();
   const auto protein_donor    = protein_layer.get_is_hbond_donor();
   const auto protein_acceptor = protein_layer.get_is_hbond_acceptor();
   const auto protein_hydro    = protein_layer.get_is_hydrophobic();
   const auto ligand_donor    = ligand_layer.get_is_hbond_donor();
   const auto ligand_acceptor = ligand_layer.get_is_hbond_acceptor();
   const auto ligand_hydro    = ligand_layer.get_is_hydrophobic();

   vinardo_preprocessed_pairs result;
   //Protein-Ligand Pairs: Cartesian product
   std::vector<vinardo_protein_ligand_pair> pl_pairs;
   if (options.include_hydrogens) {
      append_protein_ligand_pairs<true>(pl_pairs,
                                        protein_layer,
                                        ligand_layer,
                                        protein_radii,
                                        ligand_radii,
                                        protein_donor,
                                        protein_acceptor,
                                        protein_hydro,
                                        ligand_donor,
                                        ligand_acceptor,
                                        ligand_hydro);
   } else {
      append_protein_ligand_pairs<false>(pl_pairs,
                                         protein_layer,
                                         ligand_layer,
                                         protein_radii,
                                         ligand_radii,
                                         protein_donor,
                                         protein_acceptor,
                                         protein_hydro,
                                         ligand_donor,
                                         ligand_acceptor,
                                         ligand_hydro);
   }

   //Assign the result to the output struct with move so i don't have to copy the vector.
   result.protein_ligand_pairs = std::move(pl_pairs);

   //Ligand-Ligand Pairs: included only if topological distance >= 4 bonds (excluding 1-2, 1-3, 1-4 interactions, inherited from Vina in Vinardo)
   std::vector<vinardo_ligand_ligand_pair> ll_pairs; 
   //I need to compute the graph so I can asses the distance and the relative movement
   auto& ligand = ligand_layer.get_base_molecule();
   auto graph = mudock::make_graph(ligand.get_bonds(), ligand.num_atoms());
   const auto num_atoms = static_cast<std::size_t>(ligand.num_atoms());

   //Utils for the rotable edges
   mudock::fragments<mudock::static_containers> ligand_fragments(graph, ligand.get_bonds(),ligand.num_atoms());

   auto relatively_movable = precompute_rotatable_pairs(ligand_fragments, num_atoms);
   auto within_three_bonds = precompute_within_n_bonds(graph, num_atoms, 3);

   if (options.include_hydrogens) {
      append_ligand_ligand_pairs<true>(ll_pairs,
                                       ligand_layer,
                                       ligand_radii,
                                       ligand_donor,
                                       ligand_acceptor,
                                       ligand_hydro,
                                       relatively_movable,
                                       within_three_bonds);
   } else {
      append_ligand_ligand_pairs<false>(ll_pairs,
                                        ligand_layer,
                                        ligand_radii,
                                        ligand_donor,
                                        ligand_acceptor,
                                        ligand_hydro,
                                        relatively_movable,
                                        within_three_bonds);
   }
   result.ligand_ligand_pairs = std::move(ll_pairs);



   return result;
}

   
}//mudock
