/**********************************************************************
cifadptest.cpp - Tests for CIF/mmCIF anisotropic displacement parameters.
***********************************************************************/

#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/generic.h>
#include <openbabel/elements.h>

#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>

using namespace OpenBabel;
using namespace std;

static double pairAsDouble(OBAtom* atom, const string& attr)
{
  OBPairData* pd = dynamic_cast<OBPairData*>(atom->GetData(attr));
  OB_REQUIRE(pd != nullptr);
  return atof(pd->GetValue().c_str());
}

static string pairAsString(OBAtom* atom, const string& attr)
{
  OBPairData* pd = dynamic_cast<OBPairData*>(atom->GetData(attr));
  OB_REQUIRE(pd != nullptr);
  return pd->GetValue();
}

static void setPair(OBAtom* atom, const string& attr, const string& value)
{
  OBPairData* pd = new OBPairData;
  pd->SetAttribute(attr);
  pd->SetValue(value);
  atom->SetData(pd);
}

static bool near(double a, double b, double tol = 1.0e-8)
{
  return fabs(a - b) < tol;
}

static void testMMCIFAnisotropIdMatching()
{
  const string mmcif =
    "data_mmcif_anisotrop_test\n"
    "loop_\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n"
    "_atom_site.occupancy\n"
    "_atom_site.adp_type\n"
    "A001 C C1 1.0 2.0 3.0 1.0 Uani\n"
    "A002 O O1 4.0 5.0 6.0 1.0 Uani\n"
    "\n"
    "loop_\n"
    "_atom_site_anisotrop.id\n"
    "_atom_site_anisotrop.type_symbol\n"
    "_atom_site_anisotrop.U[1][1]\n"
    "_atom_site_anisotrop.U[1][2]\n"
    "_atom_site_anisotrop.U[1][3]\n"
    "_atom_site_anisotrop.U[2][2]\n"
    "_atom_site_anisotrop.U[2][3]\n"
    "_atom_site_anisotrop.U[3][3]\n"
    "A002 O 0.040 0.001 0.002 0.050 0.003 0.060\n"
    "A001 C 0.010 0.004 0.005 0.020 0.006 0.030\n";

  OBConversion conv;
  OB_REQUIRE(conv.SetInFormat("mmcif"));
  OBMol mol;
  OB_REQUIRE(conv.ReadString(&mol, mmcif));
  OB_COMPARE(mol.NumAtoms(), 2u);

  OBAtom* c = mol.GetAtom(1);
  OBAtom* o = mol.GetAtom(2);
  OB_ASSERT(near(pairAsDouble(c, "adp_U_11"), 0.010));
  OB_ASSERT(near(pairAsDouble(c, "adp_U_22"), 0.020));
  OB_ASSERT(near(pairAsDouble(c, "adp_U_33"), 0.030));
  OB_ASSERT(near(pairAsDouble(c, "adp_U_12"), 0.004));
  OB_ASSERT(near(pairAsDouble(c, "adp_U_13"), 0.005));
  OB_ASSERT(near(pairAsDouble(c, "adp_U_23"), 0.006));
  OB_ASSERT(near(pairAsDouble(o, "adp_U_11"), 0.040));
  OB_ASSERT(near(pairAsDouble(o, "adp_U_22"), 0.050));
  OB_ASSERT(near(pairAsDouble(o, "adp_U_33"), 0.060));
  OB_ASSERT(near(pairAsDouble(o, "adp_U_12"), 0.001));
  OB_ASSERT(near(pairAsDouble(o, "adp_U_13"), 0.002));
  OB_ASSERT(near(pairAsDouble(o, "adp_U_23"), 0.003));
  OB_COMPARE(pairAsString(c, "adp_source"), string("mmcif_atom_site_anisotrop"));
  OB_COMPARE(pairAsString(c, "adp_basis"), string("cif"));
  OB_COMPARE(pairAsString(c, "adp_type"), string("Uani"));
}


static void testMMCIFAnisotropBTensorConversion()
{
  const double pi = 3.141592653589793238462643383279502884;
  const double bToU = 1.0 / (8.0 * pi * pi);
  const string mmcif =
    "data_mmcif_anisotrop_b_test\n"
    "loop_\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n"
    "1 C C1 1.0 2.0 3.0\n"
    "\n"
    "loop_\n"
    "_atom_site_anisotrop.id\n"
    "_atom_site_anisotrop.B[1][1]\n"
    "_atom_site_anisotrop.B[1][2]\n"
    "_atom_site_anisotrop.B[1][3]\n"
    "_atom_site_anisotrop.B[2][2]\n"
    "_atom_site_anisotrop.B[2][3]\n"
    "_atom_site_anisotrop.B[3][3]\n"
    "1 0.80 0.08 0.16 1.60 0.24 2.40\n";

  OBConversion conv;
  OB_REQUIRE(conv.SetInFormat("mmcif"));
  OBMol mol;
  OB_REQUIRE(conv.ReadString(&mol, mmcif));
  OB_COMPARE(mol.NumAtoms(), 1u);

  OBAtom* atom = mol.GetAtom(1);
  OB_ASSERT(near(pairAsDouble(atom, "adp_U_11"), 0.80 * bToU));
  OB_ASSERT(near(pairAsDouble(atom, "adp_U_22"), 1.60 * bToU));
  OB_ASSERT(near(pairAsDouble(atom, "adp_U_33"), 2.40 * bToU));
  OB_ASSERT(near(pairAsDouble(atom, "adp_U_12"), 0.08 * bToU));
  OB_ASSERT(near(pairAsDouble(atom, "adp_U_13"), 0.16 * bToU));
  OB_ASSERT(near(pairAsDouble(atom, "adp_U_23"), 0.24 * bToU));
  OB_COMPARE(pairAsString(atom, "adp_input_type"), string("B"));
}

static void testCIFADPWriterUsesAtomSiteLabels()
{
  OBMol mol;
  mol.BeginModify();
  OBAtom* atom = mol.NewAtom();
  atom->SetAtomicNum(6);
  atom->SetVector(0.1, 0.2, 0.3);
  mol.EndModify();
  setPair(atom, "adp_U_11", "0.010");
  setPair(atom, "adp_U_22", "0.020");
  setPair(atom, "adp_U_33", "0.030");
  setPair(atom, "adp_U_12", "0.004");
  setPair(atom, "adp_U_13", "0.005");
  setPair(atom, "adp_U_23", "0.006");

  OBConversion conv;
  OB_REQUIRE(conv.SetOutFormat("cif"));
  string cif = conv.WriteString(&mol);
  size_t atomLabel = cif.find("C0");
  size_t adpLabel = cif.find("C0", cif.find("_atom_site_aniso_U_23"));
  OB_ASSERT(atomLabel != string::npos);
  OB_ASSERT(adpLabel != string::npos);

  OBConversion readConv;
  OB_REQUIRE(readConv.SetInFormat("cif"));
  OBMol reread;
  OB_REQUIRE(readConv.ReadString(&reread, cif));
  OB_COMPARE(reread.NumAtoms(), 1u);
  OB_ASSERT(near(pairAsDouble(reread.GetAtom(1), "adp_U_11"), 0.010));
}

int cifadptest(int argc, char* argv[])
{
#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif
  int choice = 1;
  if (argc > 1 && sscanf(argv[1], "%d", &choice) != 1)
    return -1;
  switch (choice) {
  case 1: testMMCIFAnisotropIdMatching(); break;
  case 2: testCIFADPWriterUsesAtomSiteLabels(); break;
  case 3: testMMCIFAnisotropBTensorConversion(); break;
  default: return -1;
  }
  return 0;
}
