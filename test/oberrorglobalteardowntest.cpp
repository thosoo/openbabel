/**********************************************************************
oberrorglobalteardowntest.cpp - Process-level global obErrorLog teardown test.
***********************************************************************/

#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/mol.h>

#include <iostream>
#include <string>

using namespace OpenBabel;

int main()
{
  const std::string cif =
    "data_oberror_teardown\n"
    "_cell_length_a 6.9827\n"
    "_cell_length_b 11.8748\n"
    "_cell_length_c 15.1296\n"
    "_cell_angle_alpha 90\n"
    "_cell_angle_beta 98.397\n"
    "_cell_angle_gamma 90\n"
    "_space_group_name_H-M_alt 'P 1'\n"
    "\n"
    "loop_\n"
    "_atom_site_label\n"
    "_atom_site_type_symbol\n"
    "_atom_site_fract_x\n"
    "_atom_site_fract_y\n"
    "_atom_site_fract_z\n"
    "_atom_site_U_iso_or_equiv\n"
    "_atom_site_adp_type\n"
    "_atom_site_occupancy\n"
    "C1 C 0.8169 0.4433 0.71499 0.0194 Uani 1\n"
    "\n"
    "loop_\n"
    "_atom_site_aniso_label\n"
    "_atom_site_aniso_U_11\n"
    "_atom_site_aniso_U_22\n"
    "_atom_site_aniso_U_33\n"
    "_atom_site_aniso_U_23\n"
    "_atom_site_aniso_U_13\n"
    "_atom_site_aniso_U_12\n"
    "C1 0.0231 0.0235 0.0118 0.0011 0.0036 0.0017\n";

  if (!obErrorLog.StartErrorWrap())
    return 1;

  OBConversion conversion;
  if (!conversion.SetInFormat("cif"))
    return 2;

  OBMol molecule;
  if (!conversion.ReadString(&molecule, cif))
    return 3;

  if (molecule.NumAtoms() != 1)
    return 4;

  // Keep text buffered in the wrapping streambuf so the global teardown path
  // exercises destruction while cerr is still wrapped.
  std::cerr << std::nounitbuf;
  std::cerr << "OpenBabel global teardown regression probe";

  // Deliberately do not call obErrorLog.StopErrorWrap(). Returning normally
  // lets C++ global teardown invoke obErrorLog's destructor.
  return 0;
}
