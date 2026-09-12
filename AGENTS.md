# Project instructions

These instructions apply throughout the primecount repository.

## Project overview

primecount is a command-line program and C/C++ library for computing π(x), the number of primes less than or equal to x. It uses highly optimized combinatorial prime-counting algorithms with OpenMP parallelization. Gourdon's algorithm is the default.

Read `README.md` once at the beginning of each new session before making code changes to understand the project's capabilities, usage, and algorithm overview. Consult relevant sections again as needed. Use the build and verification instructions in this `AGENTS.md` when developing.

## Core algorithms and required reading

primecount implements a family of combinatorial prime-counting algorithms. The Gourdon algorithm in `src/gourdon` improves on the Deléglise-Rivat algorithm in `src/deleglise-rivat`, which improves on the Lagarias-Miller-Odlyzko (LMO) algorithm in `src/lmo`.

These algorithms combine many formulas. Two of the most important contributions are the easy special leaves and the hard special leaves. Read the related mathematical papers before editing core algorithm code. Use the following mappings:

| Component | Implementation files | Required papers |
| --- | --- | --- |
| Easy special leaves | `src/deleglise-rivat/S2_easy.cpp`, `src/gourdon/AC.cpp`, `src/gourdon/AC_*.hpp` | `doc/Easy-Special-Leaves.pdf` |
| Hard special leaves | `src/deleglise-rivat/S2_hard.cpp`, `src/gourdon/D.cpp`, `src/gourdon/D_*.hpp` | `doc/Hard-Special-Leaves.pdf`, `doc/Hard-Special-Leaves-SIMD-Filtering.pdf` |
| Partial sieve function | `src/phi.cpp`, `src/phi_vector.hpp`, `src/phi_vector.cpp` | `doc/Partial-Sieve-Function.pdf` |

For other core algorithm code, identify and read the related mathematical paper before editing, using references in the source files and project documentation.

Read each relevant paper once per Codex session/chat, when work on the corresponding algorithm first begins and before the first code edit. For example, work on `src/gourdon/D.cpp` requires reading both `doc/Hard-Special-Leaves.pdf` and `doc/Hard-Special-Leaves-SIMD-Filtering.pdf`. Reuse that reading for subsequent edits in the same session; do not reread a paper for every edit or message. When work moves to another component, read any required papers that have not yet been read in the session. In a new session/chat, read the relevant papers again before editing.

Use the papers to understand the mathematics, but do not blindly copy their notation into variable names. Choose names that make the source code readable and fit existing naming conventions. Prefer shorter names when they remain clear.

## Vendored primesieve library

primecount uses the primesieve C++ library, vendored in `lib/primesieve`. Exclude this directory from project-wide changes, including requests to change all occurrences of a pattern. Do not make code changes there as part of such work.

Fixes to primesieve belong in the upstream primesieve repository. Afterwards, update the vendored copy in primecount using `scripts/update_libprimesieve.sh`. The user usually performs this update manually; leave it to them unless they ask you to perform it.

## Preserve user changes

Preserve existing user changes, including staged changes. Do not revert, overwrite, or recreate code the user intentionally changed or removed unless the requested task requires it.

## Compatibility

Keep primecount's implementation compatible with C++14 and preserve existing platform support. C++14 is only a requirement for building primecount itself: applications using its public C++ headers and linking against the library must continue to work with C++11 or later. Keep public headers compatible with C++11 and avoid propagating a C++14 requirement to library consumers through build or package metadata. Do not introduce newer language requirements, dependencies, or public API changes unless the task calls for them.

## Coding conventions

primecount has no official coding style guide that can be enforced by a tool. Infer the coding style from the file being edited and follow its existing formatting. If the file is too small or lacks examples of the code construct being written, inspect one or two other source files to determine how to format it.

Keep changes limited to the requested task. Avoid unrelated refactoring, renaming, formatting, or whitespace changes. Preserve the existing file encoding and line endings.

- Do not break a variable initialization immediately after `=` except in rare cases, such as a complex boolean initializer with many conditions.
- Prefer a multi-line `if` condition over introducing a boolean variable used only for that condition.
- Split compound preprocessor conditions in `#if` and `#elif` directives across multiple lines, with one condition per line. Use `\` line continuations and align the continued conditions.
- For SIMD instruction sets that have both `ENABLE_<ISA>` and `ENABLE_MULTIARCH_<ISA>` macros, the native `ENABLE_<ISA>` case takes precedence when both are defined. Native builds such as `-march=native` should use the SIMD implementation directly without runtime CPU-feature checks; the `ENABLE_MULTIARCH_<ISA>` case is for portable builds that require runtime dispatch. Structure the preprocessor logic with the native case first and the multiarch case in `#elif`.
- Source files that use `ENABLE_<ISA>` macros must include `<cpu_arch_macros.hpp>`. This is not required for files that use only `ENABLE_MULTIARCH_<ISA>`, since those macros are defined by the build system.
- Split overly complicated expressions, especially nested `min()`/`max()` calls combined with table lookups, into simpler intermediate calculations.
- Do not use `UINT64_C(...)` or similar integer-constant macros.

## Integer types and casts

When choosing between 64-bit and 128-bit integer arithmetic, establish the mathematical bounds using primecount's supported input range of `x <= 10^31`. Prefer 64-bit arithmetic when the result and all intermediate operations are proven to fit the chosen signed or unsigned type. Use 128-bit arithmetic where needed.

Generally keep signedness consistent within a function: use signed or unsigned integers as appropriate and avoid unnecessary mixing of the two. Mixing 64-bit and 128-bit widths is common and does not require mixing signedness.

Avoid unnecessary explicit casts and integer literal suffixes. Rely on implicit conversions and the usual arithmetic conversions when the surrounding expression already establishes the desired integer type and the conversion is safe and unambiguous. For example, prefer `65537 + j * 2` when `j` is a `uint64_t`, and prefer `svwhilelt_b64(0, active)` when the intended overload is unambiguous.

Use an explicit cast or integer literal suffix only when it affects the semantics, prevents an unsafe conversion, or is needed to select the intended overload. For example, use `1ull << 63` when the literal itself must be 64-bit before the shift, and cast before an operation when widening must occur before that operation, as in `(uint128_t(12345) << 64) | 987654321`.

When a cast is necessary, follow the convention used in the surrounding function or file. Do not introduce named C++ casts such as `static_cast<uint64_t>(x)`, `reinterpret_cast`, or similar forms. primecount uses both function-style casts such as `uint64_t(x)` and C-style casts such as `(uint64_t) x`; in C++ code, prefer `uint64_t(x)` when there is no nearby precedent. Function-style casts such as `uint64_t(x)` are C++-only, so in C files use C-style casts such as `(uint64_t) x`. Preserve nearby pointer-cast style as well, for example `(uint8_t*) sieve_.data()` when that matches the surrounding code.

## Internal utilities and the C++ standard library

Prefer primecount's internal utilities over C++ standard library equivalents. Check the project's internal headers for existing functionality before using the standard library. Generally limit standard library use to low-level operations for which primecount has no equivalent; C++ atomics and printing with `std::cout` are examples of allowed uses. Use primecount's `Vector` and `Array` instead of `std::vector` and `std::array`.

- The code frequently mixes `uint128_t`, `int128_t`, `uint64_t`, and `int64_t`. Use the internal `min()` and `max()` functions in `include/min.hpp` to avoid unnecessary explicit casts. For supported mixed integral types, put the wider type first: `min()` returns the second argument's type and `max()` returns the first argument's type. For example, `int64_t res64 = min(var_i128, var_i64);` needs no explicit cast. Respect the helpers' signedness and value-range requirements.
- Preserve existing `std::min()` and `std::max()` calls when they remain suitable for the modified code. Do not replace them merely because `min.hpp` is included or internal `min()`/`max()` is used elsewhere in the file.
- For new min/max calls, follow the convention of the surrounding snippet or function when the arguments have the same type and require no casts. Use the internal helpers when mixed types would otherwise require explicit casts. Do not add `min.hpp` solely to replace suitable standard-library calls.
- `include/min.hpp` also provides `min3(a, b, c)` and `max3(a, b, c)` for three arguments with similar mixed-type support; `min3()` returns the last argument's type and `max3()` returns the first argument's type.
- Use `in_between(lower, value, upper)` from `include/imath.hpp` instead of `std::clamp()`. Note that the lower bound is the first argument and the value being clamped is the second.

## Code comments

- Comment sparsely. Add a comment without being asked only when it helps a human reader understand an important, non-obvious aspect of the code, such as an algorithmic choice or a correctness constraint. Avoid comments that merely restate the code.
- When modifying code with an existing comment, update the comment if the code changes make it incorrect.
- When the user asks to copy a code section, copy its comments along with the code.
- When writing a multi-line comment, check the width of nearby comments and use a similar wrapping width.
- Keep comments short and compact. There is no fixed line count; use nearby comments as a guide to the expected length, especially when there are many examples nearby.
- All `*.cpp`, `*.c`, `*.hpp`, and `*.h` files contain a top-level comment describing the file and providing license and copyright information. Whenever updating one of these files, update the copyright year in that comment to the current year. If the copyright uses a year range, preserve the starting year and update the ending year.
- When creating a new `*.cpp`, `*.c`, `*.hpp`, or `*.h` file, add the standard top-level file description, license, and copyright comment using existing project files as a template and the current copyright year. This standard header is required regardless of the preference for sparse comments.

## Refactoring

After a significant code change, perform a dedicated refactoring and cleanup pass on the newly added or modified code before considering the task complete. A change is considered significant if the total number of newly added and modified existing lines is at least 30.

During this pass, inspect a few closely related source files, up to a maximum of 5, and analyze their coding conventions. Prefer the most relevant nearby or analogous implementations. Use those files to match the project's existing coding style, formatting, naming, comments, and code structure as closely as possible.

Do not refactor code outside the newly added or modified code unless it is necessary for the requested change. Do not refactor merely to satisfy this requirement. If the implementation is already simple and consistent with the surrounding code, leave it unchanged. Avoid unrelated cleanup outside the scope of the task.

When working on performance-critical core algorithms, avoid refactoring that could deteriorate performance merely to reduce code duplication or the number of lines of code. In hot inner loops, prefer keeping performance-critical code inline rather than extracting it into functions or abstractions that the compiler might fail to inline. Performance takes precedence over reducing code size or duplication in such cases.

## ChangeLog

Whenever making a notable change in a file, add a short description, ideally one line, to the root `ChangeLog` under the next primecount version. Follow the existing entry format and update the date in that version's title to the current date using `YYYY-MM-DD` format.

The user may revise or remove ChangeLog entries, including entries they consider insufficiently notable. Respect these editorial decisions: preserve their revisions and do not restore removed entries unless explicitly asked.

## Building on Windows

Use MSYS2 MinGW-w64 as the primary Windows compiler. On this PC the tools are installed in `C:\msys64\mingw64\bin`.

Prepend that directory to `PATH` in every new shell used to configure, build, or run tests. This makes both the build tools and their runtime DLLs available. Set the compiler and make program explicitly as shown below to avoid tool discovery delays.

The following PowerShell commands were verified on this PC using the `Unix Makefiles` generator from `doc/Build.md`, CMake 4.4.0, and GCC 16.1.0. Run them from the repository root. `build-agents` is an ignored build directory; reuse it when its compiler and generator match, or choose a fresh build directory when they differ.

```powershell
$env:Path = 'C:\msys64\mingw64\bin;' + $env:Path
New-Item -ItemType Directory -Force -Path build-agents | Out-Null
Set-Location build-agents
cmake .. -G 'Unix Makefiles' -DCMAKE_MAKE_PROGRAM=C:/msys64/mingw64/bin/mingw32-make.exe -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe -DCMAKE_C_COMPILER=C:/msys64/mingw64/bin/gcc.exe '-DCMAKE_CXX_FLAGS=-Wall -Wextra -pedantic -Werror' '-DCMAKE_C_FLAGS=-Wall -Wextra -pedantic -Werror' -DBUILD_TESTS=ON
if ($LASTEXITCODE -ne 0) { throw 'CMake configuration failed' }
cmake --build . --parallel 14
if ($LASTEXITCODE -ne 0) { throw 'Build failed' }
ctest
if ($LASTEXITCODE -ne 0) { throw 'Unit tests failed' }
.\primecount.exe --test
if ($LASTEXITCODE -ne 0) { throw 'Prime-counting implementation tests failed' }
.\primecount.exe 1e18 -s
if ($LASTEXITCODE -ne 0) { throw 'Large computation failed' }
```

The build defaults to `Release`. The parallel build above uses this PC's 14 logical processors; adjust the job count on other machines.

This configuration uses GNU's linker, with `CMAKE_LINKER` resolving to `C:/msys64/mingw64/bin/ld.exe`. Microsoft Visual Studio's `link.exe` is not required. Do not add a Visual Studio linker path or force `link.exe` for this MinGW build.

## Required verification

Always build and run all of the following checks after code changes, using a separate build directory:

1. Configure with tests enabled and warnings treated as errors for both C++ and C:

   ```sh
   cmake .. -DCMAKE_CXX_FLAGS="-Wall -Wextra -pedantic -Werror" -DCMAKE_C_FLAGS="-Wall -Wextra -pedantic -Werror" -DBUILD_TESTS=ON
   ```

   On Windows, use the explicit tool paths and generator in the PowerShell recipe above.

2. Build using `cmake --build . --parallel`.
3. Run `ctest` from the build directory. If a test fails, inspect its output using `ctest --output-on-failure`.
4. Run `./primecount --test` to test all prime-counting implementations.
5. Run `./primecount 1e18 -s` and verify that the result is exactly `24739954287740860`.

In PowerShell, use `.\primecount.exe` as shown above. Check every command's exit status and report any failures or checks that could not be completed.

See `doc/Build.md` for general build instructions and `test/README.md` for additional testing options.

## Windows Python process cleanup

- Do not use the Microsoft Store `python.exe` launcher for long-running or potentially cancelled Python commands. Resolve and invoke the real Python interpreter executable directly.
- Record the process IDs created by every Python command. If the command is cancelled, times out, or is otherwise terminated early, stop its entire process tree immediately and verify that it is gone.
- Never finish a turn after invoking Python without checking for Python processes created during that turn and stopping any that are still running.
- Do not terminate a running Python tool cell without performing the process-tree cleanup and verification in the same turn.
