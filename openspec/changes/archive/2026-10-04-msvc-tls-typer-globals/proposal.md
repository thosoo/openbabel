# Proposal: MSVC-compatible typer globals (TLS fix in-source)

## Why

The global typer singletons (`atomtyper`, `aromtyper` in `include/openbabel/typer.h` and `src/typer.cpp`) are declared `THREAD_LOCAL OB_EXTERN ...`. On MSVC, `OB_EXTERN` expands to a DLL interface declaration and MSVC rejects DLL interface decorations on thread-local globals (C2492), so shared (DLL) Windows builds cannot compile this source as-is. The downstream Avogadro fork currently works around this by regex-patching these two files at build time (`PatchOpenBabelMsvcTls.cmake`), which mutates the working tree after checkout. Landing the fix natively in this fork removes that build-time mutation and is the prerequisite for Avogadro pinning this fork via a git submodule.

## What Changes

- In `include/openbabel/typer.h`: wrap the `atomtyper` and `aromtyper` declarations in `#if defined(_MSC_VER)` blocks — ordinary `OB_EXTERN` declarations on MSVC, unchanged `THREAD_LOCAL OB_EXTERN` declarations elsewhere.
- In `src/typer.cpp`: wrap the `aromtyper`/`atomtyper` definitions in an `#if defined(_MSC_VER)` block — ordinary definitions on MSVC, unchanged `THREAD_LOCAL` definitions elsewhere.
- No behavioral change on non-MSVC platforms; the globals remain single DLL-interface instances on MSVC (exported, not thread-local).
- After verification, tag the commit (e.g. annotated tag `avogadro-2026.10`) so Avogadro can pin it durably.

## Capabilities

### New Capabilities
- `msvc-typer-globals`: how the global typer singletons are declared and defined so MSVC DLL builds compile (avoiding C2492) while thread-local semantics are preserved on all other platforms, with no post-checkout source mutation required.

### Modified Capabilities
- (none — first capability in this repository)

## Impact

- **Code**: `include/openbabel/typer.h` (two declaration sites), `src/typer.cpp` (one definition block).
- **Consumers**: Avogadro's `openbabel-submodule` change pins this fork; once this change is committed and tagged, that change can proceed and `PatchOpenBabelMsvcTls.cmake` is deleted on the Avogadro side.
- **ABI note**: on MSVC the two globals become non-thread-local ordinary globals; Windows consumers that relied on per-thread typer state (none known — these typers are initialized once per process) see process-global state, matching what the Avogadro patch already produces today.
