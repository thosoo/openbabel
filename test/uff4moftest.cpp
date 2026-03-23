#include "obtest.h"

#include <openbabel/forcefield.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace OpenBabel;

namespace {

struct ParamRow
{
  vector<double> dpars;
  int coord;
};

int CoordFromType(const string &type)
{
  char coord = type.size() > 2 ? type[2] : '\0';
  switch (coord) {
  case '1': return 1;
  case '2':
  case 'R': return 2;
  case '3': return 3;
  case '4': return 4;
  case '5': return 5;
  case '6': return 6;
  case '7': return 7;
  case '8': return 8;
  default: return 1;
  }
}

string GetDataPath(const string &filename)
{
  const char *dataDir = getenv("BABEL_DATADIR");
  OB_REQUIRE(dataDir != nullptr);
  return string(dataDir) + "/" + filename;
}

map<string, ParamRow> ReadParamFile(const string &filename)
{
  ifstream ifs(GetDataPath(filename).c_str());
  OB_REQUIRE(ifs.good());

  map<string, ParamRow> params;
  string line;
  while (getline(ifs, line)) {
    istringstream iss(line);
    string tag;
    iss >> tag;
    if (tag != "param")
      continue;

    string type;
    ParamRow row;
    iss >> type;
    row.coord = CoordFromType(type);
    for (int i = 0; i < 11; ++i) {
      double value;
      iss >> value;
      row.dpars.push_back(value);
    }
    params[type] = row;
  }
  return params;
}

void CheckParam(const map<string, ParamRow> &params, const string &type,
                double r1, double theta0, int coord)
{
  map<string, ParamRow>::const_iterator it = params.find(type);
  OB_REQUIRE(it != params.end());
  OB_ASSERT(fabs(it->second.dpars[0] - r1) < 1.0e-6);
  OB_ASSERT(fabs(it->second.dpars[1] - theta0) < 1.0e-6);
  OB_COMPARE(it->second.coord, coord);
}

void CheckInheritedTerms(const map<string, ParamRow> &base,
                         const map<string, ParamRow> &derived,
                         const string &baseType,
                         const string &derivedType)
{
  map<string, ParamRow>::const_iterator baseIt = base.find(baseType);
  map<string, ParamRow>::const_iterator derivedIt = derived.find(derivedType);
  OB_REQUIRE(baseIt != base.end());
  OB_REQUIRE(derivedIt != derived.end());
  for (std::size_t i = 2; i < derivedIt->second.dpars.size(); ++i)
    OB_ASSERT(fabs(baseIt->second.dpars[i] - derivedIt->second.dpars[i]) < 1.0e-12);
}

void testSeparateParameterFiles()
{
  const map<string, ParamRow> uff = ReadParamFile("UFF.prm");
  const map<string, ParamRow> uff4mof = ReadParamFile("UFF4MOF.prm");

  OB_ASSERT(uff.find("O_3_f") == uff.end());
  OB_ASSERT(uff.find("S_3_f") == uff.end());
  OB_ASSERT(uff.find("U_8f4") == uff.end());
  OB_ASSERT(uff4mof.find("C_3") == uff4mof.end());
  OB_ASSERT(uff4mof.find("O_3") == uff4mof.end());
  OB_ASSERT(uff4mof.find("Zn3+2") == uff4mof.end());

  CheckParam(uff4mof, "O_3_f", 0.634, 109.47, 3);
  CheckParam(uff4mof, "O_2_z", 0.528, 120.0, 2);
  CheckParam(uff4mof, "S_3_f", 0.854, 109.47, 3);
  CheckParam(uff4mof, "Mn8f4", 1.52, 109.47, 8);
  CheckParam(uff4mof, "U_8f4", 1.73, 109.47, 8);
  CheckParam(uff4mof, "Cu1f1", 1.24, 180.0, 1);
  CheckParam(uff4mof, "Cu2f2", 1.11, 120.0, 2);
  CheckParam(uff4mof, "Cu3f2", 1.19, 109.47, 3);
  CheckParam(uff4mof, "Pt4f2", 1.125, 90.0, 4);
  CheckParam(uff4mof, "Pd6f3", 1.19, 90.0, 6);
}

void testInheritanceAndRegistration()
{
  const map<string, ParamRow> uff = ReadParamFile("UFF.prm");
  const map<string, ParamRow> uff4mof = ReadParamFile("UFF4MOF.prm");

  CheckInheritedTerms(uff, uff4mof, "O_3", "O_3_f");
  CheckInheritedTerms(uff, uff4mof, "O_2", "O_2_z");
  CheckInheritedTerms(uff, uff4mof, "Zn3+2", "Zn3f2");
  CheckInheritedTerms(uff, uff4mof, "Pt4+2", "Pt4f2");
  CheckInheritedTerms(uff, uff4mof, "U_6+4", "U_8f4");
  CheckInheritedTerms(uff, uff4mof, "S_2", "S_3_f");

}

void testUFFRegressionView()
{
  const map<string, ParamRow> uff = ReadParamFile("UFF.prm");
  CheckParam(uff, "C_3", 0.757, 109.47, 3);
  CheckParam(uff, "C_R", 0.729, 120.0, 2);
  CheckParam(uff, "O_3", 0.658, 104.51, 3);
}

}

int uff4moftest(int argc, char* argv[])
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

  switch (choice) {
  case 1:
    testSeparateParameterFiles();
    testInheritanceAndRegistration();
    testUFFRegressionView();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
