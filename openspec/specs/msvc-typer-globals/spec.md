# msvc-typer-globals Specification

## Purpose

Defines how the global typer singletons (`atomtyper`, `aromtyper`) must be declared and defined so that MSVC shared-library (DLL) builds compile without C2492 while every other platform keeps thread-local storage, and so that building the library never requires mutating checked-out source files.

## Requirements

### Requirement: Platform-conditional typer global declarations
The `atomtyper` and `aromtyper` declarations in `include/openbabel/typer.h` SHALL be ordinary DLL-interface declarations when compiled by MSVC, and SHALL retain the `THREAD_LOCAL OB_EXTERN` form on all other compilers. The `aromtyper` and `atomtyper` definitions in `src/typer.cpp` SHALL follow the same platform split.

#### Scenario: MSVC DLL build compiles
- **WHEN** the library is built as a shared library with MSVC
- **THEN** compilation succeeds with no C2492 diagnostics in `typer.h` or `typer.cpp`

#### Scenario: Non-MSVC builds keep thread-local storage
- **WHEN** the library is built with GCC, Clang, or another non-MSVC compiler
- **THEN** `atomtyper` and `aromtyper` retain thread-local storage semantics identical to the pre-change declarations and definitions

### Requirement: Pristine-source MSVC build
Building the library from a pristine checkout SHALL NOT require any post-checkout modification of source files on any platform, including MSVC.

#### Scenario: Fresh clone builds without patching
- **WHEN** a fresh clone of this repository is configured and built with MSVC as a shared library
- **THEN** the build succeeds with no source file modified by the build system

### Requirement: Globals remain DLL-visible single instances on MSVC
On MSVC builds the two typer globals SHALL remain DLL-interface exported single instances (one per process), matching what the Avogadro build-time patch produced.

#### Scenario: Symbols exported from the DLL
- **WHEN** the MSVC shared build is inspected (e.g. with `dumpbin /exports`)
- **THEN** the `atomtyper` and `aromtyper` symbols are exported from the OpenBabel DLL
