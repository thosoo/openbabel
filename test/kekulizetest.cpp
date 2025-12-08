/**********************************************************************
kekulizetest.cpp - Regression tests for aromatic kekulization

Copyright (C) 2025 Open Babel contributors

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include "obtest.h"

#include <openbabel/kekulize.h>
#include <openbabel/obconversion.h>

#include <string>
#include <vector>

using namespace std;
using namespace OpenBabel;

static bool KekulizeFromSmiles(const std::string &smiles)
{
  OBMol mol;
  OBConversion conv;

  OB_REQUIRE(conv.SetInFormat("smi"));
  OB_REQUIRE(conv.ReadString(&mol, smiles));

  return OBKekulize(&mol);
}

int kekulizetest(int argc, char *argv[])
{
  cout << endl << "# Testing aromatic kekulization stability...  " << endl;

  vector<string> aromaticSmiles = {
      "c1ccccc1",
      "c1ccccc1N",
      "Cc1ccccc1",
  };

  for (const string &smiles : aromaticSmiles) {
    OB_ASSERT(KekulizeFromSmiles(smiles));
  }

  return 0;
}

