#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>






namespace mudock{
   // Vinardo atom types and associated lookup data derived from Smina
   enum class vinardo_atom_type : std::uint8_t {
      Hydrogen,
      PolarHydrogen,
      AliphaticCarbonXSHydrophobe,
      AliphaticCarbonXSNonHydrophobe,
      AromaticCarbonXSHydrophobe,
      AromaticCarbonXSNonHydrophobe,
      Nitrogen,
      NitrogenXSDonor,
      NitrogenXSDonorAcceptor,
      NitrogenXSAcceptor,
      Oxygen,
      OxygenXSDonor,
      OxygenXSDonorAcceptor,
      OxygenXSAcceptor,
      Sulfur,
      SulfurAcceptor,
      Phosphorus,
      Fluorine,
      Chlorine,
      Bromine,
      Iodine,
      Magnesium,
      Manganese,
      Zinc,
      Calcium,
      Iron,
      GenericMetal,
      NumTypes
   };

   //This is a subset of autodock_ff that is used as input for the lookup table
   enum class ad_lookup_table_input_type : std::uint8_t {
      H,
      HD,
      C,
      A,
      N,
      O,
      NA,
      OA,
      S,
      SA,
      P,
      F,
      Cl,
      Br,
      I,
      Mg,
      Mn,
      Zn,
      Ca,
      Fe,
      GenericMetal,
      Unsupported,
      NumTypes
   };


   struct VinardoTypeInfo {
      vinardo_atom_type vd_type;
      ad_lookup_table_input_type ad_type; //Needed for the first conversion
      float xs_radius;
      bool xs_hydrophobe;
      bool xs_donor;
      bool xs_acceptor;
   };
   
   //This tabel is based on the one in smina
   const VinardoTypeInfo vinardo_data[static_cast<std::size_t>(vinardo_atom_type::NumTypes)] = {
      {vinardo_atom_type::Hydrogen,                       ad_lookup_table_input_type::H,  0.0f, false, false, false},
      {vinardo_atom_type::PolarHydrogen,                  ad_lookup_table_input_type::HD, 0.0f, false, false, false},
      {vinardo_atom_type::AliphaticCarbonXSHydrophobe,    ad_lookup_table_input_type::C,  1.9f, true,  false, false},
      {vinardo_atom_type::AliphaticCarbonXSNonHydrophobe, ad_lookup_table_input_type::C,  1.9f, false, false, false},
      {vinardo_atom_type::AromaticCarbonXSHydrophobe,     ad_lookup_table_input_type::A,  1.9f, true,  false, false},
      {vinardo_atom_type::AromaticCarbonXSNonHydrophobe,  ad_lookup_table_input_type::A,  1.9f, false, false, false},
      {vinardo_atom_type::Nitrogen,                       ad_lookup_table_input_type::N,  1.8f, false, false, false},
      {vinardo_atom_type::NitrogenXSDonor,                ad_lookup_table_input_type::N,  1.8f, false, true,  false},
      {vinardo_atom_type::NitrogenXSDonorAcceptor,        ad_lookup_table_input_type::NA, 1.8f, false, true,  true },
      {vinardo_atom_type::NitrogenXSAcceptor,             ad_lookup_table_input_type::NA, 1.8f, false, false, true },
      {vinardo_atom_type::Oxygen,                         ad_lookup_table_input_type::O,  1.7f, false, false, false},
      {vinardo_atom_type::OxygenXSDonor,                  ad_lookup_table_input_type::O,  1.7f, false, true,  false},
      {vinardo_atom_type::OxygenXSDonorAcceptor,          ad_lookup_table_input_type::OA, 1.7f, false, true,  true },
      {vinardo_atom_type::OxygenXSAcceptor,               ad_lookup_table_input_type::OA, 1.7f, false, false, true },
      {vinardo_atom_type::Sulfur,                         ad_lookup_table_input_type::S,  2.0f, false, false, false},
      {vinardo_atom_type::SulfurAcceptor,                 ad_lookup_table_input_type::SA, 2.0f, false, false, false},
      {vinardo_atom_type::Phosphorus,                     ad_lookup_table_input_type::P,  2.1f, false, false, false},
      {vinardo_atom_type::Fluorine,                       ad_lookup_table_input_type::F,  1.5f, true,  false, false},
      {vinardo_atom_type::Chlorine,                       ad_lookup_table_input_type::Cl, 1.8f, true,  false, false},
      {vinardo_atom_type::Bromine,                        ad_lookup_table_input_type::Br, 2.0f, true,  false, false},
      {vinardo_atom_type::Iodine,                         ad_lookup_table_input_type::I,  2.2f, true,  false, false},
      {vinardo_atom_type::Magnesium,                      ad_lookup_table_input_type::Mg, 1.2f, false, true,  false},
      {vinardo_atom_type::Manganese,                      ad_lookup_table_input_type::Mn, 1.2f, false, true,  false},
      {vinardo_atom_type::Zinc,                           ad_lookup_table_input_type::Zn, 1.2f, false, true,  false},
      {vinardo_atom_type::Calcium,                        ad_lookup_table_input_type::Ca, 1.2f, false, true,  false},
      {vinardo_atom_type::Iron,                           ad_lookup_table_input_type::Fe, 1.2f, false, true,  false},
      {vinardo_atom_type::GenericMetal,                   ad_lookup_table_input_type::GenericMetal,  1.2f, false, true,  false}
   };

   //Needed to dump the type as a string
   [[nodiscard]] inline std::string_view to_string(const vinardo_atom_type type) {
      switch (type) {
         case vinardo_atom_type::Hydrogen: return "Hydrogen";
         case vinardo_atom_type::PolarHydrogen: return "PolarHydrogen";
         case vinardo_atom_type::AliphaticCarbonXSHydrophobe: return "AliphaticCarbonXSHydrophobe";
         case vinardo_atom_type::AliphaticCarbonXSNonHydrophobe: return "AliphaticCarbonXSNonHydrophobe";
         case vinardo_atom_type::AromaticCarbonXSHydrophobe: return "AromaticCarbonXSHydrophobe";
         case vinardo_atom_type::AromaticCarbonXSNonHydrophobe: return "AromaticCarbonXSNonHydrophobe";
         case vinardo_atom_type::Nitrogen: return "Nitrogen";
         case vinardo_atom_type::NitrogenXSDonor: return "NitrogenXSDonor";
         case vinardo_atom_type::NitrogenXSDonorAcceptor: return "NitrogenXSDonorAcceptor";
         case vinardo_atom_type::NitrogenXSAcceptor: return "NitrogenXSAcceptor";
         case vinardo_atom_type::Oxygen: return "Oxygen";
         case vinardo_atom_type::OxygenXSDonor: return "OxygenXSDonor";
         case vinardo_atom_type::OxygenXSDonorAcceptor: return "OxygenXSDonorAcceptor";
         case vinardo_atom_type::OxygenXSAcceptor: return "OxygenXSAcceptor";
         case vinardo_atom_type::Sulfur: return "Sulfur";
         case vinardo_atom_type::SulfurAcceptor: return "SulfurAcceptor";
         case vinardo_atom_type::Phosphorus: return "Phosphorus";
         case vinardo_atom_type::Fluorine: return "Fluorine";
         case vinardo_atom_type::Chlorine: return "Chlorine";
         case vinardo_atom_type::Bromine: return "Bromine";
         case vinardo_atom_type::Iodine: return "Iodine";
         case vinardo_atom_type::Magnesium: return "Magnesium";
         case vinardo_atom_type::Manganese: return "Manganese";
         case vinardo_atom_type::Zinc: return "Zinc";
         case vinardo_atom_type::Calcium: return "Calcium";
         case vinardo_atom_type::Iron: return "Iron";
         case vinardo_atom_type::GenericMetal: return "GenericMetal";
         case vinardo_atom_type::NumTypes: return "NumTypes";
      }
      return "Unknown";
   }
};
