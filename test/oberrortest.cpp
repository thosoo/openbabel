/**********************************************************************
oberrortest.cpp - Tests for Open Babel error handling.
***********************************************************************/

#include "obtest.h"
#include <openbabel/oberror.h>

#include <iostream>
#include <sstream>

using namespace OpenBabel;
using namespace std;

namespace
{
  class InspectableMessageHandler : public OBMessageHandler
  {
  public:
    std::streambuf* InWrapStreamBuf() const { return _inWrapStreamBuf; }
    std::streambuf* FilterStreamBuf() const { return _filterStreamBuf; }
  };
}

static void testConstructDestroyWithoutWrapping()
{
  OBMessageHandler handler;
  OB_COMPARE(handler.IsDestructing(), false);
}

static void testStartStopErrorWrapIsIdempotent()
{
  InspectableMessageHandler handler;
  std::streambuf* original = cerr.rdbuf();

  OB_REQUIRE(handler.StartErrorWrap());
  OB_ASSERT(handler.InWrapStreamBuf() == original);
  OB_ASSERT(handler.FilterStreamBuf() != nullptr);
  OB_ASSERT(cerr.rdbuf() == handler.FilterStreamBuf());

  OB_REQUIRE(handler.StopErrorWrap());
  OB_ASSERT(handler.InWrapStreamBuf() == nullptr);
  OB_ASSERT(cerr.rdbuf() == original);

  OB_REQUIRE(handler.StopErrorWrap());
  OB_ASSERT(cerr.rdbuf() == original);
}

static void testStartStopStartErrorWrap()
{
  InspectableMessageHandler handler;
  std::streambuf* original = cerr.rdbuf();

  OB_REQUIRE(handler.StartErrorWrap());
  std::streambuf* filter = handler.FilterStreamBuf();
  OB_REQUIRE(handler.StopErrorWrap());
  OB_ASSERT(cerr.rdbuf() == original);

  OB_REQUIRE(handler.StartErrorWrap());
  OB_ASSERT(handler.InWrapStreamBuf() == original);
  OB_ASSERT(handler.FilterStreamBuf() == filter);
  OB_ASSERT(cerr.rdbuf() == filter);

  OB_REQUIRE(handler.StopErrorWrap());
  OB_ASSERT(cerr.rdbuf() == original);
}

static void testActiveWrapDestructorRestoresCerr()
{
  std::streambuf* original = cerr.rdbuf();
  std::streambuf* filter = nullptr;

  {
    InspectableMessageHandler handler;
    OB_REQUIRE(handler.StartErrorWrap());
    filter = handler.FilterStreamBuf();
    OB_ASSERT(filter != nullptr);
    OB_ASSERT(cerr.rdbuf() == filter);
  }

  OB_ASSERT(cerr.rdbuf() == original);
  OB_ASSERT(cerr.rdbuf() != filter);
}

static void testOrdinaryErrorReportingStillWorks()
{
  OBMessageHandler handler;
  std::ostringstream output;

  handler.SetOutputStream(&output);
  handler.SetOutputLevel(obInfo);
  handler.ThrowError("oberrortest", "ordinary reporting", obInfo);

  OB_COMPARE(handler.GetInfoMessageCount(), 1u);
  OB_ASSERT(output.str().find("ordinary reporting") != std::string::npos);
}

int oberrortest(int argc, char* argv[])
{
  int choice = 1;
  if (argc > 1 && sscanf(argv[1], "%d", &choice) != 1)
    return -1;
  switch (choice) {
  case 1: testConstructDestroyWithoutWrapping(); break;
  case 2: testStartStopErrorWrapIsIdempotent(); break;
  case 3: testStartStopStartErrorWrap(); break;
  case 4: testActiveWrapDestructorRestoresCerr(); break;
  case 5: testOrdinaryErrorReportingStillWorks(); break;
  default: return -1;
  }
  return 0;
}
