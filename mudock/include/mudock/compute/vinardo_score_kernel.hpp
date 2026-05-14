#pragma once

#include <concepts>
#include <cstdint>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct vinardo_score_kernel {
    static constexpr char vinardo_region_name[] = "vinardo_score_kernel";
    vinardo_score_kernel(
                        const int scores_per_ligand_, 
                        const int batch_ligands_,
                        const int batch_atoms_,
                        //Ligands positions
                        const fp_type* x_scratch_b_,
                        const fp_type* y_scratch_b_,
                        const fp_type* z_scratch_b_,
                        //Protein position
                        const fp_type* prot_x_b_,
                        const fp_type* prot_y_b_,
                        const fp_type* prot_z_b_,

                        const int* vinardo_num_tors_b_,
                        //PL pairs
                        const int* pl_offsets_b_,
                        const int* pl_counts_b_,
                        const int* pl_protein_atom_idx_b_,
                        const int* pl_ligand_atom_idx_b_,
                        const fp_type* pl_radius_sum_b_,
                        const std::uint8_t* pl_hydrophobic_possible_b_,
                        const std::uint8_t* pl_hbond_possible_b_,
                        //LL pairs
                        const int* ll_offsets_b_,
                        const int* ll_counts_b_,
                        const int* ll_atom_i_idx_b_,
                        const int* ll_atom_j_idx_b_,
                        const fp_type* ll_radius_sum_b_,
                        const std::uint8_t* ll_hydrophobic_possible_b_,
                        const std::uint8_t* ll_hbond_possible_b_,

                        fp_type* scores_b_,
                        std::shared_ptr<queue_type> q_):
                        
                        scores_per_ligand(scores_per_ligand_),
                        batch_ligands(batch_ligands_),
                        batch_atoms(batch_atoms_),
                        x_scratch_b(x_scratch_b_),
                        y_scratch_b(y_scratch_b_),
                        z_scratch_b(z_scratch_b_),
                        prot_x_b(prot_x_b_),
                        prot_y_b(prot_y_b_),
                        prot_z_b(prot_z_b_),
                        vinardo_num_tors_b(vinardo_num_tors_b_),
                        pl_offsets_b(pl_offsets_b_),
                        pl_counts_b(pl_counts_b_),
                        pl_protein_atom_idx_b(pl_protein_atom_idx_b_),
                        pl_ligand_atom_idx_b(pl_ligand_atom_idx_b_),
                        pl_radius_sum_b(pl_radius_sum_b_),
                        pl_hydrophobic_possible_b(pl_hydrophobic_possible_b_),
                        pl_hbond_possible_b(pl_hbond_possible_b_),
                        ll_offsets_b(ll_offsets_b_),
                        ll_counts_b(ll_counts_b_),
                        ll_atom_i_idx_b(ll_atom_i_idx_b_),
                        ll_atom_j_idx_b(ll_atom_j_idx_b_),
                        ll_radius_sum_b(ll_radius_sum_b_),
                        ll_hydrophobic_possible_b(ll_hydrophobic_possible_b_),
                        ll_hbond_possible_b(ll_hbond_possible_b_),
                        scores_b(scores_b_),
                        q(q_) {}

    void operator()();

    vinardo_score_kernel(const vinardo_score_kernel &)            = default;
    vinardo_score_kernel(vinardo_score_kernel &&)                 = default;
    vinardo_score_kernel &operator=(const vinardo_score_kernel &) = delete;
    vinardo_score_kernel &operator=(vinardo_score_kernel &&)      = delete;

    ~vinardo_score_kernel() = default;

  private:
    const int scores_per_ligand;
    const int batch_ligands;
    const int batch_atoms;

    const fp_type* x_scratch_b;
    const fp_type* y_scratch_b;
    const fp_type* z_scratch_b;

    const fp_type* prot_x_b;
    const fp_type* prot_y_b;
    const fp_type* prot_z_b;

    const int* vinardo_num_tors_b;

    const int* pl_offsets_b;
    const int* pl_counts_b;
    const int* pl_protein_atom_idx_b;
    const int* pl_ligand_atom_idx_b;
    const fp_type* pl_radius_sum_b;
    const std::uint8_t* pl_hydrophobic_possible_b;
    const std::uint8_t* pl_hbond_possible_b;

    const int* ll_offsets_b;
    const int* ll_counts_b;
    const int* ll_atom_i_idx_b;
    const int* ll_atom_j_idx_b;
    const fp_type* ll_radius_sum_b;
    const std::uint8_t* ll_hydrophobic_possible_b;
    const std::uint8_t* ll_hbond_possible_b;

    fp_type* scores_b;

    std::shared_ptr<queue_type> q;
  };

} // namespace mudock
