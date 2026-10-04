# Design: MSVC-compatible typer globals

## Context

`include/openbabel/typer.h:70` and `:84` declare `THREAD_LOCAL OB_EXTERN OBAtomTyper atomtyper;` / `THREAD_LOCAL OB_EXTERN OBAromaticTyper aromtyper;`; `src/typer.cpp:43-44` defines both as `THREAD_LOCAL ...`. `OB_EXTERN` expands to a DLL interface attribute when consuming/building the OpenBabel DLL; MSVC rejects DLL interface decorations on thread-local globals (C2492). `THREAD_LOCAL` is defined in `include/openbabel/mol.h:32-37` (empty unless a C++11 `thread_local` path applies). The exact transformation is already specified by Avogadro's `PatchOpenBabelMsvcTls.cmake` (in the avogadro repo), which has been verifying itself with extensive content assertions on every Avogadro MSVC build. This change ports that proven transformation into the source, so the Avogadro-side script can be deleted and the fork pinned via submodule.

## Goals / Non-Goals

**Goals:**
- MSVC DLL builds compile from pristine source (no C2492)
- Byte-for-byte equivalence with what `PatchOpenBabelMsvcTls.cmake` produces, so the Avogadro-side swap is provably a no-op for the patched content
- Preserve thread-local behavior on every non-MSVC platform

**Non-Goals:**
- Rethinking whether these globals should be thread-local at all (upstream design question, out of scope)
- Fixing other MSVC issues, if any
- Changing visibility/export macros (`OB_EXTERN`, `OBAPI`) beyond the conditional forms
- Release/versioning activities beyond the single pin tag

## Decisions

1. **`#if defined(_MSC_VER)` conditionals, not macro removal.** Keep the existing `THREAD_LOCAL OB_EXTERN` form verbatim in the `#else` branch so non-MSVC translation units are bit-identical to today. Alternative (dropping TLS everywhere) would change semantics for non-Windows consumers; rejected.

2. **Mirror the patch script's exact text.** The header declarations become:
   ```c
   #if defined(_MSC_VER)
   OB_EXTERN OBAtomTyper      atomtyper;
   #else
   THREAD_LOCAL OB_EXTERN OBAtomTyper      atomtyper;
   #endif
   ```
   (same for `aromtyper`), and `typer.cpp` wraps both definitions in one `#if defined(_MSC_VER)` / `#else` / `#endif` block preserving the two-space indentation. This makes `diff` between a scratch-patched tree and this commit empty — the verification used by task 1.1 of the Avogadro change.

3. **Tag after verification, not before.** The annotated tag (e.g. `avogadro-2026.10`) is cut at the verified commit so Avogadro's submodule pin (and any future archive-URL escape hatch) references a known-good state.

## Risks / Trade-offs

- [Windows consumers lose per-thread typer state] → Accepted: the Avogadro patch (in production use) already produces exactly this; typers are init-once process-global resources.
- [Upstream drift: future upstream changes to these lines will conflict] → Accepted for a companion fork; the conditional blocks are small and localized.
- [A compiler other than MSVC defines `_MSC_VER` (e.g. clang-cl)] → clang-cl intentionally masquerades as MSVC and accepts non-TLS DLL-interface globals; matches the patch script's behavior, unchanged.

## Migration Plan

1. Apply the two-file edit (spec + design above).
2. Verify content equivalence against a scratch run of Avogadro's `PatchOpenBabelMsvcTls.cmake`.
3. Build-verify: at minimum a non-MSVC build (this machine, GCC) to prove the `#else` path is untouched; MSVC verification happens downstream in Avogadro CI (the submodule pin will exercise it).
4. Commit, push, cut annotated tag `avogadro-2026.10`, push the tag.
5. Downstream: Avogadro `openbabel-submodule` change pins the tag.

Rollback: revert the commit; nothing else depends on it until the Avogadro submodule lands.

## Open Questions

- None.
