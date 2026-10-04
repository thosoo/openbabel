# Tasks: MSVC-compatible typer globals

## 1. Source change

- [x] 1.1 Edit `include/openbabel/typer.h`: wrap the `atomtyper` declaration (line ~70) and the `aromtyper` declaration (line ~84) each in `#if defined(_MSC_VER)` (plain `OB_EXTERN`) / `#else` (`THREAD_LOCAL OB_EXTERN`, verbatim) / `#endif`; verify with `git diff` that the `#else` branches preserve the original lines byte-for-byte — DONE. `git diff` shows both `#else` lines unchanged (context); MSVC branches add plain `OB_EXTERN` declarations.
- [x] 1.2 Edit `src/typer.cpp`: wrap the `aromtyper`/`atomtyper` definitions (lines ~43-44) in a single `#if defined(_MSC_VER)` (plain definitions) / `#else` (`THREAD_LOCAL` definitions, verbatim) / `#endif` block; verify with `git diff` — DONE. Single block around both definitions; `#else` branch byte-identical to original.

## 2. Equivalence verification

- [x] 2.1 Run Avogadro's patch script against a pristine scratch copy of this tree (`cmake -DOPENBABEL_SOURCE_DIR=<scratch> -P PatchOpenBabelMsvcTls.cmake` from the avogadro repo) and confirm `diff -r` between the scratch-patched files and this working tree shows no differences in `typer.h`/`typer.cpp`; verify the script reports "already present" against the new tree — DONE. Ran against a `git archive` pristine copy: script patched both files; `diff` vs working tree is IDENTICAL for both. Re-run against the working tree reports "OpenBabel MSVC TLS patch already present" with no mutation.
- [x] 2.2 Configure + build this repo with GCC (non-MSVC path) from a clean build dir and confirm it compiles and `git status` is clean (no source mutation); verify the `THREAD_LOCAL` `#else` branch compiled — DONE. Clean `/tmp/ob-build` configure + `cmake --build --target openbabel` succeeded (linked libopenbabel.so). `typer.cpp.o` built. `git status` shows only the two intended edits (no build mutation).

## 3. Commit and tag

- [x] 3.1 Commit the two-file change with a message describing the C2492 rationale; push to `master` — DONE. Commit `f30fefc71a914f98ba4bcb2e2f283d1f9c68ea47` on `master`; pushed to `origin master` (HEAD).
- [x] 3.2 Cut annotated tag `avogadro-2026.10` at the commit, push the tag; verify with `git rev-parse avogadro-2026.10^{commit}` — DONE. Annotated tag cut at `f30fefc71`, pushed to `origin`; `git rev-parse avogadro-2026.10^{commit}` = `f30fefc71a914f98ba4bcb2e2f283d1f9c68ea47`.
- [x] 3.3 Downstream handoff: report the tag/commit to the avogadro repo so the `openbabel-submodule` change can pin it (its tasks 2.1+) — DONE. Recorded the available pin (commit `f30fefc71`, tag `avogadro-2026.10`) plus fork-side verification evidence in the avogadro `openbabel-submodule` change design (`openspec/changes/openbabel-submodule/design.md`) so tasks 2.1+ can pin directly.
