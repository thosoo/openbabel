#include "obtest.h"
#include <openbabel/obconversion.h>
#include <openbabel/kekulize.h>

using namespace OpenBabel;
using namespace std;

int avogadrotest(int argc, char* argv[])
{
#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif

  cout << endl << "# Kekulizing Avogadro examples..." << endl;

  const char* lines[] = {
    "c1ccccc1  phenyl",
    "c1ccccc1N  aniline",
    "Cc1ccccc1  toluene",
    nullptr
  };

  OBConversion conv;
  OB_REQUIRE(conv.SetInFormat("smi"));

  OBMol mol;
  unsigned int testCount = 0;

  for(const char** line = lines; *line; ++line) {
    mol.Clear();
    OB_REQUIRE(conv.ReadString(&mol, *line));
    if (OBKekulize(&mol))
      cout << "ok " << ++testCount << "\n";
    else
      cout << "not ok " << ++testCount << " # failed to kekulize\n";
  }

  return 0;
}
