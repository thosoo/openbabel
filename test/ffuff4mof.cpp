/**********************************************************************
ffuff4mof.cpp - Regression tests for the UFF4MOF force field
***********************************************************************/

#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
#include <dlfcn.h>
#endif

#include <openbabel/atom.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/plugin.h>

using namespace std;
using namespace OpenBabel;

namespace {

OBMol MakeZnAcetate()
{
  OBMol mol;

  OBAtom* zn = mol.NewAtom();
  zn->SetAtomicNum(30);
  zn->SetVector(0.0, 0.0, 0.0);

  OBAtom* o1 = mol.NewAtom();
  o1->SetAtomicNum(8);
  o1->SetVector(1.95, 0.0, 0.0);
  OBAtom* c1 = mol.NewAtom();
  c1->SetAtomicNum(6);
  c1->SetVector(3.15, 0.0, 0.0);
  OBAtom* o2 = mol.NewAtom();
  o2->SetAtomicNum(8);
  o2->SetVector(4.10, 0.90, 0.0);
  OBAtom* c2 = mol.NewAtom();
  c2->SetAtomicNum(6);
  c2->SetVector(3.45, -1.45, 0.0);

  OBAtom* o3 = mol.NewAtom();
  o3->SetAtomicNum(8);
  o3->SetVector(-1.95, 0.0, 0.0);
  OBAtom* c3 = mol.NewAtom();
  c3->SetAtomicNum(6);
  c3->SetVector(-3.15, 0.0, 0.0);
  OBAtom* o4 = mol.NewAtom();
  o4->SetAtomicNum(8);
  o4->SetVector(-4.10, -0.90, 0.0);
  OBAtom* c4 = mol.NewAtom();
  c4->SetAtomicNum(6);
  c4->SetVector(-3.45, 1.45, 0.0);

  OBAtom* o5 = mol.NewAtom();
  o5->SetAtomicNum(8);
  o5->SetVector(0.0, 2.20, 0.0);

  mol.AddBond(1, 2, 1);
  mol.AddBond(2, 3, 1);
  mol.AddBond(3, 4, 2);
  mol.AddBond(3, 5, 1);
  mol.AddBond(1, 6, 1);
  mol.AddBond(6, 7, 1);
  mol.AddBond(7, 8, 2);
  mol.AddBond(7, 9, 1);
  mol.AddBond(1, 10, 1);

  mol.SetTitle("Zn acetate fragment");
  return mol;
}

OBMol MakeZnThioether()
{
  OBMol mol;

  OBAtom* zn = mol.NewAtom();
  zn->SetAtomicNum(30);
  zn->SetVector(0.0, 0.0, 0.0);

  OBAtom* s1 = mol.NewAtom();
  s1->SetAtomicNum(16);
  s1->SetVector(2.25, 0.0, 0.0);
  OBAtom* c1 = mol.NewAtom();
  c1->SetAtomicNum(6);
  c1->SetVector(3.75, 0.90, 0.0);
  OBAtom* c2 = mol.NewAtom();
  c2->SetAtomicNum(6);
  c2->SetVector(3.75, -0.90, 0.0);

  OBAtom* s2 = mol.NewAtom();
  s2->SetAtomicNum(16);
  s2->SetVector(-2.25, 0.0, 0.0);
  OBAtom* c3 = mol.NewAtom();
  c3->SetAtomicNum(6);
  c3->SetVector(-3.75, 0.90, 0.0);
  OBAtom* c4 = mol.NewAtom();
  c4->SetAtomicNum(6);
  c4->SetVector(-3.75, -0.90, 0.0);

  mol.AddBond(1, 2, 1);
  mol.AddBond(2, 3, 1);
  mol.AddBond(2, 4, 1);
  mol.AddBond(1, 5, 1);
  mol.AddBond(5, 6, 1);
  mol.AddBond(5, 7, 1);

  mol.SetTitle("Zn thioether fragment");
  return mol;
}

void* LoadForceFieldPluginHandle()
{
#if defined(__unix__) || defined(__APPLE__)
  static void* handle = nullptr;
  if (handle != nullptr)
    return handle;

  const char* libDir = getenv("BABEL_LIBDIR");
  if (libDir == nullptr || *libDir == '\0')
    return nullptr;

  stringstream ss(libDir);
  string dir;
  while (getline(ss, dir, ':')) {
    if (dir.empty())
      continue;
    const string pluginPath = dir + "/plugin_forcefields.so";
    handle = dlopen(pluginPath.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if (handle != nullptr)
      return handle;
  }
  return nullptr;
#else
  return nullptr;
#endif
}

OBForceField* GetForceFieldBySymbol(const char* id)
{
  OBForceField* ff = OBForceField::FindForceField(id);
  if (ff != nullptr)
    return ff;

#if defined(__unix__) || defined(__APPLE__)
  void* handle = LoadForceFieldPluginHandle();
  if (handle == nullptr)
    return nullptr;

  const char* symbol = strcmp(id, "UFF4MOF") == 0
    ? "_ZN9OpenBabel20theForceFieldUFF4MOFE"
    : "_ZN9OpenBabel16theForceFieldUFFE";
  void* sym = dlsym(handle, symbol);
  return reinterpret_cast<OBForceField*>(sym);
#else
  return nullptr;
#endif
}

bool CheckType(OBMol& mol, int idx, const char* expected)
{
  OBAtom* atom = mol.GetAtom(idx);
  return atom != nullptr && strcmp(atom->GetType(), expected) == 0;
}

bool IsFiniteEnergy(OBForceField* ff)
{
  const double energy = ff->Energy(false);
  return std::isfinite(energy);
}

bool TestZnCarboxylateSetup()
{
  OBMol mol = MakeZnAcetate();

  OBForceField* ff = GetForceFieldBySymbol("UFF4MOF");
  if (ff == nullptr)
    return false;
  ff->SetLogLevel(OBFF_LOGLVL_NONE);

  if (!ff->Setup(mol) || !ff->GetAtomTypes(mol))
    return false;

  return IsFiniteEnergy(ff);
}

bool TestAcetateEnergyDiffersFromUFF()
{
  OBMol molUFF = MakeZnAcetate();
  OBMol molUFF4MOF = MakeZnAcetate();

  OBForceField* uff = GetForceFieldBySymbol("UFF");
  OBForceField* uff4mof = GetForceFieldBySymbol("UFF4MOF");
  if (uff == nullptr || uff4mof == nullptr)
    return false;
  uff->SetLogLevel(OBFF_LOGLVL_NONE);
  uff4mof->SetLogLevel(OBFF_LOGLVL_NONE);

  if (!uff->Setup(molUFF) || !uff4mof->Setup(molUFF4MOF) || !uff->GetAtomTypes(molUFF) || !uff4mof->GetAtomTypes(molUFF4MOF))
    return false;

  const double uffEnergy = uff->Energy(false);
  const double uff4mofEnergy = uff4mof->Energy(false);
  return std::isfinite(uffEnergy) && std::isfinite(uff4mofEnergy) && fabs(uffEnergy - uff4mofEnergy) > 1.0e-6;
}

bool TestSulfurEnergyDiffersFromUFF()
{
  OBMol mol = MakeZnThioether();

  OBMol molUFF = mol;
  OBForceField* uff = GetForceFieldBySymbol("UFF");
  OBForceField* uff4mof = GetForceFieldBySymbol("UFF4MOF");
  if (uff == nullptr || uff4mof == nullptr)
    return false;
  uff->SetLogLevel(OBFF_LOGLVL_NONE);
  uff4mof->SetLogLevel(OBFF_LOGLVL_NONE);

  if (!uff->Setup(molUFF) || !uff4mof->Setup(mol))
    return false;

  const double uffEnergy = uff->Energy(false);
  const double uff4mofEnergy = uff4mof->Energy(false);
  return std::isfinite(uffEnergy) && std::isfinite(uff4mofEnergy) && fabs(uffEnergy - uff4mofEnergy) > 1.0e-6;
}

int RunPart(int part)
{
  int currentTest = 0;
  bool ok = false;

  LoadForceFieldPluginHandle();
  OBPlugin::LoadAllPlugins();
  cout << "# Testing UFF4MOF force field..." << endl;

  switch (part) {
  case 1:
    ok = GetForceFieldBySymbol("UFF4MOF") != nullptr;
    cout << (ok ? "ok " : "not ok ") << ++currentTest << "# plugin registration" << endl;
    ok = TestZnCarboxylateSetup();
    cout << (ok ? "ok " : "not ok ") << ++currentTest << "# zinc carboxylate setup yields finite UFF4MOF energy" << endl;
    break;
  case 2:
    ok = TestAcetateEnergyDiffersFromUFF();
    cout << (ok ? "ok " : "not ok ") << ++currentTest << "# UFF4MOF energy differs from UFF on a zinc carboxylate fragment" << endl;
    ok = TestSulfurEnergyDiffersFromUFF();
    cout << (ok ? "ok " : "not ok ") << ++currentTest << "# UFF4MOF energy differs from UFF on a zinc thioether fragment" << endl;
    break;
  default:
    cout << "Test number " << part << " does not exist!" << endl;
    return -1;
  }

  cout << "1.." << currentTest << endl;
  return 0;
}

}

int ffuff4mof(int argc, char* argv[])
{
#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif

  int choice = 1;
  if (argc > 1) {
    if (sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  return RunPart(choice);
}
