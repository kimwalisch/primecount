# Retuning the StatusS2 progress model (WorkModelS2)

This guide is for an AI agent on the dual-socket AMD EPYC 7642 server
(Zen 2, 2 × 48 cores). It explains how to measure the work distribution
of the hard special leaves algorithms, fit the constants of primecount's
progress model and transfer them into the source code. The algorithms
are D (Gourdon, `src/gourdon/D.cpp`), S2_hard (Deléglise-Rivat,
`src/deleglise-rivat/S2_hard.cpp`) and the S2 formula of LMO
(Lagarias-Miller-Odlyzko, `src/lmo/pi_lmo_parallel.cpp`).

The guide is self-contained: the scripts and the logging patch are in
the appendices. The method and all settings are those of the tuning run
of 2026-11-04, whose results are in [section 13](#13-reference-results-2026-11-04).

Workflow:

1. Check the server: SMT off, idle, tools installed ([section 3](#3-server-setup)).
2. Build an instrumented binary in a separate checkout ([section 4](#4-build-the-instrumented-binary)).
3. For each algorithm, in the order D, S2_hard, LMO: collect data
   ([section 7](#7-collect-data)), check it ([section 8](#8-check-the-data)),
   run two fits ([section 9](#9-fit)), report and **wait for the user's
   confirmation** ([section 10](#10-checkpoint-report)).
4. After the user has confirmed all three: transfer the constants
   ([section 11](#11-transfer-the-constants)) and validate them
   ([section 12](#12-validate)).

## Contents

1. [Overview](#1-overview)
2. [Rules](#2-rules)
3. [Server setup](#3-server-setup)
4. [Build the instrumented binary](#4-build-the-instrumented-binary)
5. [NUMA: CPU pinning](#5-numa-cpu-pinning)
6. [Long jobs in screen](#6-long-jobs-in-screen)
7. [Collect data](#7-collect-data)
8. [Check the data](#8-check-the-data)
9. [Fit](#9-fit)
10. [Checkpoint report](#10-checkpoint-report)
11. [Transfer the constants](#11-transfer-the-constants)
12. [Validate](#12-validate)
13. [Reference results (2026-11-04)](#13-reference-results-2026-11-04)
14. [Appendix A: collect.py](#appendix-a-collectpy)
15. [Appendix B: fit.py](#appendix-b-fitpy)
16. [Appendix C: logging patch](#appendix-c-logging-patch)

## 1. Overview

**Purpose of the model.** The three algorithms sieve the interval
`[0, sieve_limit[` in chunks assigned by `LoadBalancerS2`.
`StatusS2::getPercent(low)` estimates the share of the total work below
position `low`. `LoadBalancerS2::remaining_secs()` uses it to size the
chunks near the end of the computation, clamping the percentage to
>= 10%. The `--status` option prints it.

**Model** (`src/StatusS2.cpp`), with `L = sieve_limit`:

```
s          = in_between(1, leaves_start, L)
sieve      = sieve_work * (L / s)^sieve_exponent
work(low)  = sieve * low / L + (low > s ? 1 - (1 + (low - s) / (leaves_scale * s))^-leaves_shape : 0)
percent    = 100 * work(low) / work(L), clamped to [0, 100]
```

The first term is evenly distributed sieving work. The second is the
special leaves: their work density decays like a Lomax distribution,
O(n^-(1 + leaves_shape)), with a total of 1. `leaves_start` is not the
position of the first leaf but where the leaves' work becomes
significant; about 0.1% to a few percent of the work is below it.

**Ground truth.** The logging patch records every chunk `[low, high[`
with its CPU time. For a finished run, the true progress at position
`pos` is the CPU time of all work below `pos` divided by the total CPU
time, interpolating linearly within chunks. It does not depend on
scheduling, so multi-threaded (MT) runs measure it. Repetitions are
combined by the per-position median.

**Parameters** (names as in fit.py):

| Parameter | Meaning |
|---|---|
| `start_*` | Multipliers of the `leaves_start` formula. D: `max(min(a * x^(2/5), b * x / y^(4/3)), c * x / y^(3/2))`, S2_hard: `min(a * x / y^(3/2), b * y)`, LMO: `a * x / y^2` |
| `leaves_scale` | Lomax scale, relative to `leaves_start` |
| `leaves_shape` | Lomax shape |
| `sieve_work` | Sieving work at the default alpha, relative to the leaves' work of 1 |
| `sieve_exponent` | Dependence of the sieving share on `sieve_limit / leaves_start` |
| `alpha_exponent` | `sieve_work` is multiplied by `(default_alpha / alpha)^alpha_exponent` |

**Where the constants live.** D, as of 2026-11-04:

```cpp
// src/gourdon/D.cpp
double leaves_start = 1.6473 * std::pow(xd, 2 / 5.0);
leaves_start = std::min(leaves_start, 1.1504 * xd / std::pow(yd, 4 / 3.0));
leaves_start = std::max(leaves_start, 2.0473 * xd / (yd * std::sqrt(yd)));
// ...
double sieve_work = 0.0149 * std::pow(default_alpha / alpha, 1.082);
WorkModelS2 model(leaves_start, sieve_work, 1.9307, 0.41271, 0.34966);
```

S2_hard and LMO have the same structure with their own formula. The
constructor order is `WorkModelS2(leaves_start, sieve_work,
leaves_scale, leaves_shape, sieve_exponent)`, which differs from
fit.py's parameter order. fit.py prints the C++ form of its result.

The same constants are in two more places, keep all three in sync:

- `test/StatusS2.cpp`: the `models` array, one row per algorithm in
  constructor order with `leaves_start = 0`.
- `CURRENT` in fit.py ([Appendix B](#appendix-b-fitpy)). The replica
  checks compare it with the values the instrumented binary logs.

## 2. Rules

1. **Checkpoints.** Tune one algorithm at a time, in the order D,
   S2_hard, LMO. After the checks and fits of an algorithm, report
   ([section 10](#10-checkpoint-report)) and wait for the user's explicit
   confirmation before collecting data for the next algorithm. Don't
   change primecount's source code before the user has confirmed the
   constants of all three algorithms.
2. **Also stop and report, without waiting for the algorithm to
   finish,** if SMT is active (unless the user approved the `--no-smt`
   fallback, [section 3.1](#31-smt-off)), a data check fails
   ([section 8](#8-check-the-data)), the collection runs out of budget
   before the largest default-alpha size has its repetitions, runs fail
   for any reason other than the time limit, or a build or patch fails.
3. **Long jobs in screen.** Run every job that may take more than a few
   minutes in a detached screen session ([section 6](#6-long-jobs-in-screen)).
4. **Idle server.** Nothing else may run on the server during collection
   and validation.
5. **No installs.** Don't install or update packages, report missing
   tools to the user.
6. **No commits.** For code changes follow `AGENTS.md`.
7. **Don't mix data** from different setups in one fit: SMT on vs off,
   different builds (e.g. with and without `-mpopcnt`), or different
   versions of the logging patch.
8. **Short reports.** The user prefers a short summary with a clear
   recommendation over long tables.

## 3. Server setup

### 3.1 SMT off

```bash
cat /sys/devices/system/cpu/smt/active   # must print 0
nproc                                    # must print 96
```

If SMT is active, ask the user to turn it off (needs root):

```bash
echo off | sudo tee /sys/devices/system/cpu/smt/control
```

SMT must stay off until the validation ([section 12](#12-validate)) is
done. Then the user can turn it back on:
`echo on | sudo tee /sys/devices/system/cpu/smt/control`.

Why: with SMT, two threads share each core. D's 192-thread runs used
1.84× the CPU time per unit of work, and their work distribution differed
from single-threaded runs by 2.44 pp. The constants must suit CPUs
without SMT (Apple M-series, Intel Arrow Lake, most AWS instances).
Data collected with SMT must not be used.

Fallback without root: `collect.py --no-smt` pins MT runs to one hardware
thread per physical core. This mode was not used in the 2026 run, so
prefer turning SMT off. Use it only if the user can't turn SMT off and
approves the fallback; then active SMT is not a reason to stop, but
every `collect.py` run must print `MT runs: pinned, 96 threads` in its
first line ([section 7](#7-collect-data)), and the reports must say
that `--no-smt` was used.

### 3.2 Idle server

```bash
uptime                                         # 1-minute load average < 1
screen -ls                                     # no running tuning jobs
ps -eo pcpu,pid,user,etime,comm --sort=-pcpu | head
```

### 3.3 Tools

Python 3 with NumPy and SciPy (used: Python 3.14.4, NumPy 2.3.5,
SciPy 1.16.3), git, CMake, GCC with OpenMP, screen, taskset and lscpu.

```bash
python3 -c "import numpy, scipy; print(numpy.__version__, scipy.__version__)"
which git cmake g++ screen taskset lscpu
```

### 3.4 Directories

```bash
export PC=~/Desktop/primecount           # primecount repository (the code to tune)
export TUNE=~/Desktop/primecount-tuning  # separate checkout: binary, scripts, data
```

If `$TUNE` exists, ask the user whether to delete it or use another
name. Shell variables may not persist between your commands: set both
in every new shell. The screen jobs below read `$TUNE` from the
environment, hence the `export`.

## 4. Build the instrumented binary

Use a separate checkout so the logging patch can never end up in a
commit. The clone contains the committed state of `$PC`.

```bash
git clone -q "$PC" "$TUNE" && cd "$TUNE"
git log --oneline -1    # the commit being tuned, note it in the reports

# Extract collect.py, fit.py and the logging patch from this guide
python3 - "$PC/doc/StatusS2-tuning.md" "$TUNE" <<'EOF'
import os, re, sys
fence = "`" * 3
text = open(sys.argv[1]).read()
pattern = r"<!-- file: (\S+) -->\n%s[a-z]*\n(.*?)\n%s\n" % (fence, fence)
for name, body in re.findall(pattern, text, re.S):
    if name.endswith(".patch"):
        # restore blank context lines if an editor removed their space
        body = re.sub(r"(?m)^$", " ", body)
    with open(os.path.join(sys.argv[2], name), "w") as f:
        f.write(body + "\n")
    print("wrote", name)
EOF

git apply primecount-s2-tuning-log.patch
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-mpopcnt"
cmake --build build --parallel
```

The extraction prints three `wrote` lines. `-mpopcnt` slightly improves
performance on x64; use the same build for all three algorithms.

Quick check: the following must print D(x, y), and the log must start
with a header line like the one below. At the default alpha,
`sieve_work` equals the constant in D.cpp because the alpha factor is 1.

```bash
PRIMECOUNT_S2_LOG=$TUNE/quick.log ./build/primecount 1e18 -D
head -c 230 $TUNE/quick.log
# x 1e+18 y 67144000 sieve_limit 14893363517 leaves_start 26107945.559411947 sieve_work 0.0149 leaves_scale 1.9307000000000001 ...
```

**What the patch does**, in case it must be re-created for newer code:

- `LoadBalancerS2.hpp`: adds a destructor and the members `log_chunks_`,
  `log_path_`, `log_x_`, `log_model_` (a copy of the `WorkModelS2`),
  `chunk_frac_` and `x15_`.
- Constructor: reads `PRIMECOUNT_S2_LOG=<file>`; logging is off unless
  it is set. `PRIMECOUNT_S2_CHUNK_FRAC=<f>` is honoured only with 1
  thread and switches to fixed chunks of `f * low`.
- `get_work()`: records each finished chunk (`low`, `high`, init
  seconds, seconds, start and stop time, `getPercent(low)`, thread
  number). In fixed-chunk mode it also assigns the chunks.
- Destructor: appends one block per `LoadBalancerS2` to the log file. A
  header line `# x .. y .. sieve_limit .. leaves_start .. sieve_work ..
  leaves_scale .. leaves_shape .. sieve_exponent .. threads .. wall ..
  chunks ..` is followed by one line per chunk:
  `low high init_secs secs start stop percent thread`. fit.py reads the
  header by key and the chunk lines by column.

## 5. NUMA: CPU pinning

**The issue.** collect.py runs sizes whose all-thread run takes less
than 5 s with fewer threads ([section 7](#7-collect-data)). In the 2026
run these runs were pinned with taskset to physical cores interleaved
over both sockets, e.g. 18 threads on 9 + 9 cores. LMO's 1e17
default-alpha runs then became bimodal, reproducibly:

- About half of the runs used 8–27% more CPU time, all of it in the
  first ~35% of the work.
- The extra time was spread over most threads (1.05–1.8× each), as
  expected from remote-memory (NUMA) penalties.
- Noise rose to 2.6–2.8 pp, and the median curve shifted by 2.5 pp.

The same runs unpinned were all clean: noise 0.01 pp, CPU time
58.2–58.3 s.

**The fix** (in Appendix A):

- Reduced-thread MT runs are not pinned. With `--no-smt` they are pinned
  to the physical cores of one socket first.
- Single-threaded runs stay pinned, one process per physical core.
- All-thread runs are not pinned; with `--no-smt` they are pinned to all
  physical cores.

**Detection.** fit.py's bimodality check ([section 8](#8-check-the-data))
flags runs that use more than 5% more CPU time than the fastest run of
their configuration, and prints how much the median curve moves without
them. A single slow run among 6 or more is harmless: it moved the median
by at most 0.02 pp in 2026. Replace a configuration if its median moves
by 0.2 pp or more.

**Replacing a bad configuration** without touching the original logs.
The example is LMO's `lmo_1e17_def_mt`; run it from `$TUNE`, in screen
if it takes more than a few minutes.

**Step 1.** Replay the configuration's recorded commands into a new
directory, without CPU pinning for MT runs. Every rerun must print
status 0 and the original result.

```bash
python3 - data_lmo lmo_1e17_def_mt data_lmo_rerun <<'EOF'
# Replay the recorded runs of one configuration into a new directory,
# without CPU pinning for MT runs (ST runs keep their single CPU)
import json, os, subprocess, sys
src, config, out = sys.argv[1:]
os.makedirs(out, exist_ok=True)
for line in open(os.path.join(src, "runs.jsonl")):
    r = json.loads(line)
    if r["tag"].rsplit("_r", 1)[0] != config or r["status"] != 0: continue
    cmd = r["cmd"]
    if cmd[0] == "taskset" and r["mode"] != "st": cmd = cmd[3:]
    env = dict(os.environ, PRIMECOUNT_S2_LOG=os.path.join(out, r["tag"] + ".log"))
    if r["mode"] == "st": env["PRIMECOUNT_S2_CHUNK_FRAC"] = "0.02"
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    rec = dict(r, cmd=cmd, status=p.returncode, result=p.stdout.strip())
    with open(os.path.join(out, "runs.jsonl"), "a") as f: f.write(json.dumps(rec) + "\n")
    print(r["tag"], "status", p.returncode, "result", rec["result"], "(original %s)" % r["result"])
EOF
```

**Step 2.** Create a merged directory of symlinks: all original logs
except the replaced configuration, plus the reruns.

```bash
mkdir data_lmo_fixed
for f in data_lmo/lmo_*.log; do
  case "$f" in */lmo_1e17_def_mt_r*) ;; *) ln -s ../$f data_lmo_fixed/;; esac
done
for f in data_lmo_rerun/*.log; do ln -s ../$f data_lmo_fixed/; done
```

**Step 3.** Use the merged directory for the checks and both fits,
including the MT-only directory in [section 9](#9-fit).

## 6. Long jobs in screen

Run every job that may take more than a few minutes in a detached screen
session: collection, cross-validated fits, builds, validation runs. The
jobs then survive if the user's SSH or VS Code connection drops.

```bash
screen -ls    # don't start a second copy of a running job
screen -dmS tune_D_collect bash -c 'cd "$TUNE" && python3 collect.py --algo D --binary build/primecount --out data_D --budget-min 110 --rep-scale 1.5 > collect_D.out 2>&1; echo "exit=$?" >> collect_D.out'
```

- **Names:** `tune_<algo>_<job>`, e.g. `tune_D_collect`, `tune_D_fit_all`,
  `tune_D_fit_mt`, `tune_build`, `tune_validate_D`.
- **Logs:** each job writes all output to a log file ending with
  `exit=<status>`.
- **Never attach** (`screen -r`); read the log files instead.
- **Reports:** every report lists the running sessions and their log
  files, so the user or a new session can find them.
- **Interrupted collections:** don't rerun collect.py into the same
  `--out` directory, because it deletes logs with the same names.
  Collect into a new directory (e.g. `data_D_2`) and merge with
  `fit.py --data data_D data_D_2`.

**Wait by polling.** Run this as a background command and check
occasionally, don't busy-wait:

```bash
while screen -ls | grep -q tune_D_collect; do sleep 60; done; tail -3 collect_D.out
```

**Watch for failures** while a collection runs. The first command must
print nothing, the second 0:

```bash
grep -E "TIMEOUT|FAILED|Traceback|aborted" collect_D.out
grep -vc '"status": 0' data_D/runs.jsonl
```

## 7. Collect data

One algorithm at a time, each in screen ([section 6](#6-long-jobs-in-screen)),
run from `$TUNE`:

```bash
python3 collect.py --algo D       --binary build/primecount --out data_D   --budget-min 110 --rep-scale 1.5 > collect_D.out 2>&1
python3 collect.py --algo S2-hard --binary build/primecount --out data_S2  --budget-min 110 --rep-scale 1.5 > collect_S2.out 2>&1
python3 collect.py --algo lmo     --binary build/primecount --out data_lmo --budget-min 60  --rep-scale 3   > collect_lmo.out 2>&1
```

The first output line must be:

```
SMT active: 0, physical cores: 96 (48/48 per socket), MT runs: primecount default threads, reduced-thread runs: not pinned
```

Otherwise stop the job and report. With the approved `--no-smt`
fallback ([section 3.1](#31-smt-off)) the line shows `SMT active: 1` and
`MT runs: pinned, 96 threads, reduced-thread runs: pinned to one socket first`
instead.

**Phases**, in priority order. The budget (`--budget-min`) is checked
before each repetition, each ST run and each non-default-alpha run. The
calibration runs and the first run of each reduced-thread and ST size
always run, bounded only by the per-run time limits, because the later
phases depend on their runtimes; the budget can therefore be exceeded
by up to a few calibration runs.

| Phase | What | Settings |
|---|---|---|
| 1. Calibration | Default alpha, all threads, x = 10^12, 10^13, ... | Stops when the next size is predicted to exceed 0.9 × 600 s (growth >= 3× per decade). Runs < 5 s are logged as mode `cal` and are not fitted |
| 2. Threads per size | All-thread run >= 5 s: all 96 threads. Shorter: k = 96 × runtime / 5 s threads, not pinned, corrected once if the first run takes < 3 s. k < 2: single-threaded | `--target-secs 5` |
| 3. Default alpha, MT | All-thread sizes largest first, then the reduced-thread sizes | Repetitions by runtime: < 2 s: 10, < 10 s: 8, < 30 s: 6, < 120 s: 4, < 300 s: 3, else 2, times `--rep-scale` |
| 4. Single-threaded (ST) | Sizes too small for 2 threads, default and all alpha factors, plus the smallest MT size for the ST/MT check | Fixed chunks of 2% of `low`, 48 concurrent processes each pinned to its own physical core. Repetitions 20/12/8/5 for ST runtimes < 0.1 s / 1 s / 10 s / longer, times `--rep-scale`; a third of that for non-default alphas |
| 5. Non-default alphas, MT | Same thread counts as the default alpha. Runtime ratio per factor measured at the smallest MT size, then the larger sizes largest first | Skips runs predicted (× 1.2) to exceed 0.9 × 300 s. Half the default repetitions |

**Time limits:** 600 s for default-alpha runs (`--limit-default`), 300 s
for non-default alphas (`--limit-alpha`), 120 s for the first ST run of
each size (`--st-max-secs`, which decides whether the size is run ST at
all). The parallel ST runs get 5 × 120 s = 600 s, because they include
the non-default alphas, which can be ~3.4× slower than the default
alpha. Runs over the limit are killed and their logs deleted, since an
unfinished run has no total.

**Failures:** if `primecount` exits with a nonzero status, collect.py
deletes the run's log, records the status and the end of stderr in
`runs.jsonl`, prints `FAILED`, lets running ST jobs finish without
starting new ones and exits with status 1. Then stop and report
([section 2](#2-rules)).

`runs.jsonl` records every run: tag, mode, threads, command, wall time,
status, result, stderr. A calibration run long enough to be used as MT
data is renamed from `_cal_r1.log` to `_mt_r1.log` and gets a second
record with the new tag and mode (field `promoted_from`).

**Alpha factors** (multiples of the default alpha, clamped to
[1, x^(1/6)], skipped if within 5% of the default):

- D: 0.15, 0.3, 0.6, 1.6, 3, 6. The `x / y^(3/2)` term applies for
  factors < ~0.21 at 1e20 (< ~0.15 at 1e23), the `x / y^(4/3)` term for
  factors > ~1.4.
- S2_hard: 0.25, 0.5, 0.75, 1.5, 2.5, 4. The `y` term applies for
  factors < ~0.93 at 1e17, ~0.75 at 1e20, ~0.68 at 1e22.
- LMO: 0.45, 0.7, 1.4, 2, 3.

**Expected results** (2026-11-04, SMT off, `-mpopcnt`):

| | D | S2_hard | LMO |
|---|---|---|---|
| Calibration, all threads | 1e19 1.41 s, 1e20 5.16 s, 1e21 19.4 s, 1e22 74.3 s, 1e23 283 s | 1e19 2.86 s, 1e20 11.0 s, 1e21 42.8 s, 1e22 168 s; 1e23 (~650 s) skipped | 1e16 0.24 s, 1e17 0.97 s, 1e18 4.39 s |
| Threads per size | 96: 1e20–1e23, 27: 1e19, 7: 1e18, 2: 1e17, ST: <= 1e16 | 96: 1e20–1e22, 54: 1e19, 14: 1e18, 4: 1e17, ST: <= 1e16 | 84: 1e18, 18: 1e17, 4: 1e16, ST: <= 1e15 |
| Default-alpha MT repetitions | 12 (1e17–1e20), 9, 6, 4 (1e21–1e23) | 12 (1e17–1e19), 9, 6, 4 (1e20–1e22) | 24 (1e16–1e18) |
| Non-default alphas, MT | 6 factors at 1e17–1e22 | 6 factors at 1e17–1e21, 3 at 1e22 | 5 factors at 1e16–1e18 |
| Runs, time | 583 runs, 84 min | 543 runs, 76 min | 675 runs, 30 min |

Progress lines in the output: the calibration line, one line per
reduced-thread size (`1e19: 27 threads, 4.41s per run`), at the end
`threads per size: ...` and `done, elapsed ... min`, then `exit=0`.

## 8. Check the data

```bash
python3 fit.py --algo D --data data_D --check-only
```

It prints, in order:

0. **Incomplete logs.** A log whose chunks don't cover
   `[0, sieve_limit[` without gaps, whose chunk count differs from the
   header's `chunks`, or with invalid or zero CPU times prints
   `warning: incomplete log <file>: <reason>` and is not used. Any such
   line is a failed data check: stop and report.
1. **Default-alpha summary** per size: repetitions, noise and wall time.
   The noise (median absolute deviation of single runs, max over the
   curve) must be <= ~1 pp at the largest sizes; it was <= 0.06 pp in
   2026. Otherwise add repetitions in a new directory and merge.
   MT configurations with a median wall time < 2 s (`--mt-min-wall`)
   are marked "too short, not fitted": thread startup dominates such
   runs.
2. **Replica checks.** These must print
   `replica checks OK (N configs), getPercent() max diff ...` with a
   difference of about 5e-11 pp (the log's precision). They
   compare the logged `WorkModelS2` values with fit.py's replica
   (`CURRENT`, `leaves_start` formula, default alpha) and the logged
   `getPercent()` of every chunk with StatusS2's formula. On a mismatch,
   fit.py's model doesn't match the binary: stop and report.
3. **Bimodality check.** Must print `bimodality check OK` or only
   `note:` lines. On a `warning:` (median moves >= 0.2 pp), replace the
   configuration ([section 5](#5-numa-cpu-pinning)).
4. **ST vs MT** median curves at the smallest MT size (runs of ~5 s):
   max diff <= ~2 pp. It was 0.21–0.29 pp in 2026. If larger, stop and
   report: SMT or NUMA contention changes the work distribution.

## 9. Fit

Two fits per algorithm, in parallel, each in screen. The (b′) fit uses
a directory of symlinks to all non-ST logs.

```bash
cd "$TUNE"
mkdir data_D_mt
for f in data_D/D_*.log; do case "$f" in *_st_*) ;; *) ln -s ../$f data_D_mt/;; esac; done
screen -dmS tune_D_fit_all bash -c 'cd "$TUNE" && python3 fit.py --algo D --data data_D --out fit_D.json > fit_D.txt 2>&1; echo "exit=$?" >> fit_D.txt'
screen -dmS tune_D_fit_mt bash -c 'cd "$TUNE" && python3 fit.py --algo D --data data_D_mt --mt-min-wall 2 --out fit_D_mt.json > fit_D_mt.txt 2>&1; echo "exit=$?" >> fit_D_mt.txt'
```

- **(a) all data:** ST runs at small x plus all MT runs, 7–12
  decades. This pins down the x-dependence and usually extrapolates
  better.
- **(b′) MT only:** the conditions of real use, but only 3–7 decades.
  Constants that control the x-scaling or rarely used branches may come
  out unstable.

For S2_hard use `--algo S2-hard`, `data_S2` and `S2-hard_*.log`; for
LMO `--algo lmo`, `data_lmo` and `lmo_*.log`. If configurations were
replaced ([section 5](#5-numa-cpu-pinning)), use the merged directory.
The fits took 0.5–6 minutes in 2026, with D's (a) fit the longest.
Since then each leave-decade-out fold is fitted from scratch (see
*Cross-validation* below), which makes the fits roughly 2–4× slower:
expect up to about 20 minutes.

**Method** (fit.py, NumPy and SciPy):

- **Objective** per configuration: work-weighted mean |model − truth|
  in pp, plus 10 × the work-weighted mean |log2(remaining-time ratio)|
  for 10–99% true progress (`--rem-weight 10`). The remaining time is
  what `remaining_secs()` uses.
- **Weighting:** (n − 10)^1 per decade (`--x-weight-power 1`); the
  default alpha counts 3× (`--default-weight 3`); each
  (n, default/non-default) group shares its weight.
- **Optimizer:** `scipy.optimize.minimize` with Nelder–Mead in
  log-parameter space, then Powell. It starts from `CURRENT` plus 6
  random restarts (`--starts 6`).
- **Cross-validation:** *leave-decade-out* predicts each decade from a
  fit without it. *Extrapolated* predicts the two largest decades from a
  fit on the smaller ones. Both start from `CURRENT` with the same
  restarts as the full fit, never from the full fit's result, so the
  held-out data can't influence them through the start point. (In the
  2026 run the leave-decade-out folds started from the full fit; the
  extrapolated fits, which led to the same decisions, did not.)

**Output:**

- The objective and the constants, current vs fitted.
- A table per configuration: max / mean pp / remaining-time log2 error
  for current, fitted, leave-decade-out and extrapolated; the
  wall-clock check `wall cur/fit`; the `lag`.
- The max error at the default alpha for the two largest sizes, overall
  and for true progress >= 10%.
- The C++ form of the fitted constants (5 significant digits).
- The JSON (`--out`): `params`, `current`, `extrapolation_params`, the
  objectives, `default_alpha_errors_max_mean_rem` and `cpp`.

**Decision, per algorithm:**

1. Compare (a) and (b′) at the default alpha for the two largest sizes,
   in both the leave-decade-out and the extrapolated columns. Choose
   (b′) only if its errors are at least as good as (a)'s in both;
   otherwise choose (a).
2. Accept the chosen fit only if:
   - its leave-decade-out errors at the two largest sizes are lower
     than the current constants' (the in-sample column is not enough);
   - its extrapolated errors are not much worse than in-sample, and
     `extrapolation_params` are close to `params`, especially
     `sieve_exponent`. A large difference means the x-dependence isn't
     determined, so accuracy beyond the largest size is uncertain;
   - the wall-clock check is only a weak warning: report it, but don't
     reject a fit because of it.

   If no fit passes, recommend keeping the current constants.
3. `remaining_secs()` clamps the *estimated* percentage to >= 10%, so
   errors don't affect the load balancer while the estimate stays
   <= 10%. Underestimates early in the computation (like D's in 2026)
   therefore only affect the status display; don't change the model to
   fix them. An overestimate that pushes the estimate above 10% while
   the true progress is still below 10% does change the remaining-time
   estimate: report it if it is larger than ~2 pp.
4. The `lag` (median and range over repetitions) is the load balancer's
   in-flight lead and doesn't depend on the model; it was +2.7 to +5.7 pp
   at the largest sizes in 2026. Report it, but don't tune the model
   around it.
5. Non-default alphas matter only when the user sets `--alpha`. Report
   the worst configurations and the largest regressions against the
   current constants.

## 10. Checkpoint report

Keep it short:

```markdown
## <Algorithm>: I recommend (a) | (b′) | keeping the current constants

**Data:** N runs, 0 failures, M min; sizes 1e12–1e23; threads per size {...}.
Checks: logs complete, replica OK, bimodality OK, ST vs MT 0.27 pp at 1e17, noise <= 0.03 pp.

**Mean error at the default alpha (pp):**

| | Current | (a) leave-decade-out | (a) extrapolated | (b′) leave-decade-out | (b′) extrapolated |
|---|---|---|---|---|---|
| 1e22 | | | | | |
| 1e23 | | | | | |

**Why:** one to three bullets (decision rule, stability of the constants).
**Cost:** the worst non-default alphas, the max error and where it occurs.
**Constants:** the C++ line printed by fit.py.
**Running jobs:** none, or screen session and log file.
```

Then stop and wait for the user's decision.

## 11. Transfer the constants

Only after the user has confirmed the constants of all three algorithms.

1. Follow `AGENTS.md`. Before editing D.cpp and S2_hard.cpp, read
   `doc/Hard-Special-Leaves.pdf` and
   `doc/Hard-Special-Leaves-SIMD-Filtering.pdf`.
2. In `$PC`, for each algorithm, replace the multipliers in the
   `leaves_start` formula, the constants in the `sieve_work` line and
   the last three `WorkModelS2(...)` arguments with the values fit.py
   printed. Update the multipliers in the existing `leaves_start`
   comment above the formula. Don't add comments, and don't name tuning
   files or scripts in comments.
3. `test/StatusS2.cpp`: update the `models` rows,
   `{ 0, sieve_work, leaves_scale, leaves_shape, sieve_exponent }`.
4. Update `CURRENT` in `$TUNE/fit.py` and in this guide
   ([Appendix B](#appendix-b-fitpy)), and the constants in
   [section 1](#1-overview). Don't change
   [section 13](#13-reference-results-2026-11-04): it documents the
   2026-11-04 run. Add a new dated section after it with the new run's
   chosen fits, constants and errors in the same format.
5. `ChangeLog`: a one-line entry under the next version, and set that
   version's date to today.
6. Build an instrumented binary of the new code and run the replica
   check (below).
7. Run the verification in `AGENTS.md` (separate build directory,
   `-Werror`, ctest, `./primecount --test`, `./primecount 1e18 -s`).

**Instrumented binary of the new code**, in a second checkout (in
screen):

```bash
git -C "$PC" diff HEAD > "$TUNE/new-constants.diff"   # staged and unstaged changes
git clone -q "$PC" "$TUNE-new" && cd "$TUNE-new"
git apply "$TUNE/new-constants.diff" && git apply "$TUNE/primecount-s2-tuning-log.patch"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-mpopcnt"
cmake --build build --parallel
```

**Replica check** of the new binary against the updated `CURRENT`
(takes about a minute). All three algorithms must print
`replica checks OK (...)` with a `getPercent()` difference of about
5e-11 pp.

```bash
cd "$TUNE" && rm -rf replica_new && mkdir replica_new
python3 - "$TUNE-new/build/primecount" replica_new <<'EOF'
import os, subprocess, sys
from collect import default_alpha, trunc3, iroot, ALGO_OPT
binary, out = sys.argv[1:]
for algo, ns, factors in (("D", (15, 17, 18), (0.15, 3)), ("S2-hard", (15, 17), (0.25, 2.5)), ("lmo", (14, 16), (0.45, 2))):
    for n in ns:
        x = 10 ** n
        for f in (None,) + factors:
            alpha = "def" if f is None else trunc3(min(max(default_alpha(algo, x) * f, 1), iroot(x, 6)))
            log = os.path.join(out, "%s_1e%d_%s_mt_r1.log" % (algo, n, alpha))
            cmd = [binary, "1e%d" % n, ALGO_OPT[algo]] + ([] if f is None else ["-a", str(alpha)])
            subprocess.run(cmd, env=dict(os.environ, PRIMECOUNT_S2_LOG=log), check=True, capture_output=True)
EOF
for a in D S2-hard lmo; do python3 fit.py --algo $a --data replica_new --check-only | grep -E "replica|warning"; done
```

## 12. Validate

Compare the old binary (`$TUNE/build/primecount`, old constants) with
the new one (`$TUNE-new/build/primecount`), running each in screen.

**Step 1.** Create the run script in `$TUNE`:

```bash
cat > validate.sh <<'EOF'
#!/bin/bash
# Usage: bash validate.sh <D|S2-hard|lmo> <old binary> <new binary> <reps> <n> [<n> ...]
# Default-alpha runs of the old and the new instrumented binary,
# alternating their order to avoid a run-order bias
algo="$1"; old="$2"; new="$3"; reps="$4"; shift 4
mkdir -p "val_${algo}_old" "val_${algo}_new"
for n in "$@"; do
  for r in $(seq 1 "$reps"); do
    order="old new"; [ $((r % 2)) -eq 0 ] && order="new old"
    for b in $order; do
      bin="$old"; [ "$b" = new ] && bin="$new"
      PRIMECOUNT_S2_LOG="val_${algo}_$b/${algo}_1e${n}_def_mt_r$r.log" "$bin" "1e$n" "--$algo" > /dev/null || exit 1
      echo "$algo 1e$n r$r $b done"
    done
  done
done
EOF
```

**Step 2.** Run it, e.g. for D. Sizes: D 1e20–1e23, S2_hard
1e21–1e22, LMO 1e18 with 5 repetitions. D takes 1.2 minutes per run at
1e22 and 5 minutes at 1e23, so ask the user before including them.

```bash
screen -dmS tune_validate_D bash -c 'cd "$TUNE" && bash validate.sh D build/primecount "$TUNE-new/build/primecount" 3 20 21 22 23 > validate_D.out 2>&1; echo "exit=$?" >> validate_D.out'
```

**Step 3.** Compare. The old logs need the old constants, which are the
`current` field of the fit's JSON:

```bash
python3 fit.py --algo D --data val_D_new --no-cv | tee validate_D_new.txt
python3 fit.py --algo D --data val_D_old --no-cv --current fit_D.json:current | tee validate_D_old.txt
```

Report:

- the position errors (`current` column), which should be close to the
  fit's in-sample errors;
- the wall-clock check and the lag;
- the runtime: the median `wall` per size from the summary lines. It
  must not regress beyond the noise (~0.5%), since the model affects the
  chunk sizes through `remaining_secs()`.

**Step 4.** Ask the user to turn SMT back on ([section 3.1](#31-smt-off)).

## 13. Reference results (2026-11-04)

Code tuned: commit `6ab61d7`, built with `-mpopcnt`, SMT off. The
constants were committed in "Tune StatusS2 work model per algorithm".

**Chosen fits and constants** (old → new):

| | D: (a) all data | S2_hard: (b′) MT only | LMO: (a) all data |
|---|---|---|---|
| `leaves_start` multipliers | x^(2/5): 1.5 → 1.6473; x/y^(4/3): 1.125 → 1.1504; x/y^(3/2): 1.5 → 2.0473 | x/y^(3/2): 1.5 → 1.4064; y: 0.6 → 0.45799 | x/y^2: 1.5 → 1.3352 |
| `leaves_scale` | 2.4667 → 1.9307 | 2.4667 → 4.3158 | 2.4667 → 1.1764 |
| `leaves_shape` | 0.5 → 0.41271 | 0.5 → 0.61257 | 0.5 → 0.36355 |
| `sieve_work` | 0.070424 → 0.0149 | 0.070424 → 0.076243 | 0.070424 → 0.47309 |
| `sieve_exponent` | 0.16 → 0.34966 | 0.16 → 0.1403 | 0.16 → 0.030408 |
| `alpha_exponent` | 1.05 → 1.082 | 0.8 → 0.9068 | 2.7 → 2.5849 |
| Objective | 3.12 → 2.46 | 4.77 → 3.45 | 7.70 → 3.36 |

**Mean error at the default alpha** (pp; LDO = leave-decade-out, ext = extrapolated):

| Size | Old constants | (a) LDO / ext | (b′) LDO / ext |
|---|---|---|---|
| D 1e22 | 1.48 | **0.78 / 0.95** | 0.45 / 0.94 |
| D 1e23 | 2.05 | **1.07 / 1.16** | 0.90 / 1.33 |
| S2_hard 1e21 | 2.26 | 1.59 / 2.31 | **1.16 / 1.58** |
| S2_hard 1e22 | 2.40 | 1.98 / 2.63 | **1.56 / 1.91** |
| LMO 1e17 | 2.59 | **0.60 / 1.11** | 0.52 / 4.93 |
| LMO 1e18 | 2.26 | **1.07 / 1.54** | 0.86 / 9.10 |

**Notes and lessons:**

- D: (a) was chosen because (b′) extrapolated worse at 1e23 and its
  `x / y^(3/2)` multiplier was unstable (0.13 fitted, 1.70
  extrapolated). Cost: alpha factors 0.15 and 0.3 got 2–4 pp worse. The
  fitted max error (−2.4 pp) is at ~2.4% progress (work below
  `leaves_start`). It is an underestimate, so the estimate stays below
  10% and only the status display is affected.
- S2_hard: (b′) was better in every cross-validated number, while (a)
  extrapolated worse than the old constants. The largest error
  (+5.6 pp at ~22% progress at 1e22) sits where the true work density
  roughly doubles at `low = y`, which the model can't represent. The
  largest alphas (factor 4) have ~10–11 pp errors with any constants.
  1e23 doesn't fit the 10-minute limit.
- LMO: (b′) falls apart when extrapolated, because the MT data covers
  only 1e16–1e18. LMO only supports x up to ~9e18. The 1e17
  default-alpha runs were replaced by unpinned reruns because of the
  NUMA issue ([section 5](#5-numa-cpu-pinning)). Sieving is a much
  larger share of LMO's work than the old model assumed.
- SMT: in a first run with SMT on, neither fit of D beat the old
  constants at 1e23 and the ST/MT check failed (2.44 pp). Without SMT,
  both fits beat them clearly.
- Short MT runs with many threads (< 2 s) deviated 3–47 pp from the
  single-threaded truth, hence the 5 s target with fewer threads.
- A single slow run per configuration is harmless (median shift
  <= 0.02 pp); half the runs slow is not (2.5 pp).
- Lag at the largest sizes: D +2.7 pp, S2_hard +3.1 pp, LMO +5.7 pp.

## Appendix A: collect.py

<!-- file: collect.py -->
```python
#!/usr/bin/env python3
"""
Collect hard special leaves work-distribution data for tuning the
StatusS2 progress model (WorkModelS2) of one algorithm.

Requires a primecount binary built with primecount-s2-tuning-log.patch.
Every run writes the chunks of LoadBalancerS2 to a log file, see
fit.py for the analysis.

Phases (in priority order). The budget is checked before each
repetition, ST run and non-default-alpha run; the calibration runs and
the first run of each reduced-thread and ST size always run (bounded by
the per-run time limits), since the later phases depend on them. A run
that fails (nonzero exit status) aborts the collection with status 1.
 1. Calibration: one multi-threaded (MT) run with all threads and the
    default alpha for x = 10^n, n = min_n.., until the next run would
    exceed the time limit.
 2. Default alpha, MT, repeated (more repetitions for short runs).
    Sizes whose all-thread runs take less than --target-secs (5 s)
    use fewer threads so that each run takes about --target-secs:
    very short runs with many threads are dominated by thread startup
    and timing noise. These runs are not pinned (with --no-smt: pinned
    to the physical cores of one socket first). Pinning a few threads
    to cores on both sockets made runs bimodal (NUMA): half of them
    used up to 27% more CPU time in the first third of the work.
 3. Sizes too small even for 2 threads: single-threaded (ST) with
    fixed 2% chunks, many repetitions, run in parallel pinned to
    distinct physical cores. The smallest MT size is also run ST to
    check ST/MT consistency.
 4. Non-default alphas, MT (same thread counts), large x first.

Usage example:
  ./collect.py --algo D --binary ./build/primecount --out data_D --budget-min 110 --rep-scale 1.5
"""
import argparse, json, math, os, subprocess, sys, time

ALGO_OPT = {"D": "--D", "S2-hard": "--S2-hard", "lmo": "--lmo"}
# LMO supports only 64-bit x
MAX_N = {"D": 23, "S2-hard": 23, "lmo": 18}
# Non-default alphas as multiples of the default alpha. They cover
# the branches of the leaves_start formulas, e.g. D uses its
# x/y^(3/2) term only for factors < ~0.21 at 1e20 (< ~0.15 at 1e23)
# and its x/y^(4/3) term for factors > ~1.4, S2_hard uses its y term
# for factors < ~0.75 at 1e20.
ALPHA_FACTORS = {"D": [0.15, 0.3, 0.6, 1.6, 3, 6],
                 "S2-hard": [0.25, 0.5, 0.75, 1.5, 2.5, 4],
                 "lmo": [0.45, 0.7, 1.4, 2, 3]}

def iroot(x, n):
    r = int(round(x ** (1.0 / n)))
    while r ** n > x: r -= 1
    while (r + 1) ** n <= x: r += 1
    return r

def trunc3(a):
    return math.floor(a * 1000) / 1000

def default_alpha(algo, x):
    """Replica of get_alpha_gourdon/deleglise_rivat/lmo(x) in util.cpp"""
    L = math.log(x)
    if algo == "lmo":
        a = 0.001103 * L * L - 0.00896211 * L + 1.00404
    elif algo == "S2-hard":
        a = (0.078173 * L + 1 if x <= 1e9 else
             0.00148918 * L ** 3 - 0.0691909 * L ** 2 + 1.00165 * L + 0.372253)
    else:
        if x <= 3e11:
            a = 0.156346 * L + 1
        else:
            a = (0.021078652642134526 * L ** 3 - 2.0334362096119443 * L ** 2
                 + 67.97860096709101 * L - 758.0125793471067)
            if int(iroot(x, 3) * a) > 386083 * 386083 - 1:
                a *= 0.7
    x16 = iroot(x, 6)
    return trunc3(min(max(a, 1), x16))

def physical_cpus():
    """One logical CPU per physical core (avoid SMT siblings), only
    CPUs this process may use, as one list per socket."""
    allowed = os.sched_getaffinity(0)
    try:
        out = subprocess.run(["lscpu", "-p=CPU,CORE,SOCKET"], capture_output=True, text=True).stdout
        seen, by_sock = set(), {}
        for line in out.splitlines():
            if line.startswith("#"): continue
            cpu, core, sock = line.split(",")[:3]
            if int(cpu) in allowed and (core, sock) not in seen:
                seen.add((core, sock)); by_sock.setdefault(int(sock), []).append(int(cpu))
        socks = [by_sock[k] for k in sorted(by_sock)]
        if socks: return socks
    except Exception:
        pass
    return [sorted(allowed)]

class Collector:
    def __init__(self, args):
        self.a = args
        self.t0 = time.time()
        self.budget = args.budget_min * 60
        os.makedirs(args.out, exist_ok=True)
        self.index = os.path.join(args.out, "runs.jsonl")
        socks = physical_cpus()
        # Interleaved across sockets, a prefix spreads over all
        # sockets: used for the concurrent single-threaded runs
        self.cpus = [lst[i] for i in range(max(map(len, socks))) for lst in socks if i < len(lst)]
        # Socket by socket, a prefix fills one socket first: used for
        # reduced-thread MT runs with --no-smt (avoids NUMA effects)
        self.socket_cpus = [c for lst in socks for c in lst]
        # --no-smt: pin MT runs to one hardware thread per physical
        # core, so that no two threads share a core (as on CPUs
        # without SMT). Not needed if SMT is turned off in the kernel.
        self.mt_cpus = ",".join(map(str, sorted(self.cpus))) if args.no_smt else None
        self.mt_threads = args.threads or (len(self.cpus) if args.no_smt else 0)
        # thread count of all-thread runs, used to size reduced runs
        self.full_threads = self.mt_threads or len(os.sched_getaffinity(0))
        self.threads = {}   # n -> MT thread count used for that size
        try:
            smt = open("/sys/devices/system/cpu/smt/active").read().strip()
        except OSError:
            smt = "unknown"
        print("SMT active: %s, physical cores: %d (%s per socket), MT runs: %s, reduced-thread runs: %s" % (
              smt, len(self.cpus), "/".join(str(len(s)) for s in socks),
              "pinned, %d threads" % self.mt_threads if args.no_smt else
              "%d threads" % self.mt_threads if self.mt_threads else "primecount default threads",
              "pinned to one socket first" if args.no_smt else "not pinned"), flush=True)

    def remaining(self):
        return self.budget - (time.time() - self.t0)

    def tag(self, n, alpha, mode, rep):
        return "%s_1e%d_%s_%s_r%d" % (self.a.algo, n, alpha, mode, rep)

    def cmd(self, n, alpha, mode, threads=None):
        c = [self.a.binary, "1e%d" % n, ALGO_OPT[self.a.algo]]
        if alpha != "def":
            c += ["-a", str(alpha)]
        if mode == "st":
            c += ["-t1"]
        elif threads:
            c += ["-t%d" % threads]
        elif self.mt_threads:
            c += ["-t%d" % self.mt_threads]
        return c

    def taskset(self, mode, threads=None, cpu=None):
        """CPU pinning prefix of a run's command (empty: not pinned)"""
        if cpu is not None:
            return ["taskset", "-c", str(cpu)]
        if threads and threads < self.full_threads:
            # Never split a few threads over both sockets (NUMA, see
            # above). Unpinned runs were verified to be clean.
            if self.a.no_smt:
                return ["taskset", "-c", ",".join(map(str, sorted(self.socket_cpus[:threads])))]
            return []
        if mode in ("mt", "cal") and self.mt_cpus:
            return ["taskset", "-c", self.mt_cpus]
        return []

    def run(self, n, alpha, mode, rep, limit, cpu=None, threads=None):
        """Run one computation, returns wall seconds or None on timeout.
        A failed run (nonzero exit status) aborts the collection.
        mode: "cal" (all-thread calibration), "mt" or "st". MT runs with
        a reduced thread count are not pinned, with --no-smt they are
        pinned to physical cores of one socket first."""
        tag = self.tag(n, alpha, mode, rep)
        log = os.path.join(self.a.out, tag + ".log")
        if os.path.exists(log):
            os.remove(log)
        env = dict(os.environ, PRIMECOUNT_S2_LOG=log)
        if mode == "st":
            env["PRIMECOUNT_S2_CHUNK_FRAC"] = "0.02"
        c = self.taskset(mode, threads, cpu) + self.cmd(n, alpha, mode, threads)
        t = time.time()
        try:
            p = subprocess.run(c, env=env, capture_output=True, text=True, timeout=limit)
            wall, status, res, err = time.time() - t, p.returncode, p.stdout.strip(), p.stderr.strip()
        except subprocess.TimeoutExpired:
            wall, status, res, err = None, "timeout", "", ""
        if status != 0 and os.path.exists(log):
            os.remove(log)   # incomplete or failed runs are useless
        rec = dict(tag=tag, algo=self.a.algo, n=n, alpha=alpha, mode=mode, rep=rep,
                   threads=1 if mode == "st" else (threads or self.full_threads),
                   wall=wall if status == 0 else None, status=status, result=res,
                   stderr=err[-2000:], cmd=c)
        self.last_rec = rec
        with open(self.index, "a") as f:
            f.write(json.dumps(rec) + "\n")
        if status == "timeout":
            print("%-32s TIMEOUT" % tag, flush=True)
            return None
        if status != 0:
            print("%-32s FAILED (exit status %s): %s" % (tag, status, err[-500:]), flush=True)
            sys.exit("collection aborted: %s failed, see runs.jsonl" % tag)
        print("%-32s %.2fs" % (tag, wall), flush=True)
        return wall

    def run_parallel(self, jobs, limit):
        """ST jobs (n, alpha, rep), pinned to distinct physical cores.
        After a failed run (nonzero exit status) no new jobs are
        started; the running ones finish, then the collection aborts."""
        procs, pending, failed = {}, list(jobs), []
        cpus = self.cpus[:self.a.st_parallel]
        free = list(cpus)
        while pending or procs:
            while pending and free and not failed and self.remaining() > limit:
                n, alpha, rep = pending.pop(0)
                cpu = free.pop(0)
                tag = self.tag(n, alpha, "st", rep)
                log = os.path.join(self.a.out, tag + ".log")
                if os.path.exists(log): os.remove(log)
                env = dict(os.environ, PRIMECOUNT_S2_LOG=log, PRIMECOUNT_S2_CHUNK_FRAC="0.02")
                c = ["taskset", "-c", str(cpu)] + self.cmd(n, alpha, "st")
                p = subprocess.Popen(c, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                procs[p] = (cpu, time.time(), n, alpha, rep, tag, log, c)
            if not procs: break   # done, out of budget or failed
            time.sleep(0.05)
            for p in list(procs):
                cpu, t, n, alpha, rep, tag, log, c = procs[p]
                killed = False
                if p.poll() is None and time.time() - t > limit:
                    p.kill(); p.wait(); killed = True
                if p.poll() is not None:
                    wall = time.time() - t
                    status = "timeout" if killed else p.returncode
                    out, err = p.communicate()
                    if status != 0 and os.path.exists(log): os.remove(log)
                    rec = dict(tag=tag, algo=self.a.algo, n=n, alpha=alpha, mode="st", rep=rep,
                               threads=1, wall=wall if status == 0 else None, status=status,
                               result=out.strip(), stderr=err.strip()[-2000:], cmd=c)
                    with open(self.index, "a") as f: f.write(json.dumps(rec) + "\n")
                    if status == "timeout":
                        print("%-32s TIMEOUT" % tag, flush=True)
                    elif status != 0:
                        print("%-32s FAILED (exit status %s): %s" % (tag, status, err.strip()[-500:]), flush=True)
                        failed.append(tag)
                    del procs[p]; free.append(cpu)
        if failed:
            sys.exit("collection aborted: %d ST run(s) failed, e.g. %s, see runs.jsonl" % (len(failed), failed[0]))

    def mt_reps(self, wall):
        return max(1, round(self.a.rep_scale * self._mt_reps(wall)))

    @staticmethod
    def _mt_reps(wall):
        # Short runs are noisy: repeat them more often
        if wall < 2: return 10
        if wall < 10: return 8
        if wall < 30: return 6
        if wall < 120: return 4
        if wall < 300: return 3
        return 2

    def st_reps(self, cpu_secs):
        return max(1, round(self.a.rep_scale * self._st_reps(cpu_secs)))

    @staticmethod
    def _st_reps(cpu_secs):
        if cpu_secs < 0.1: return 20
        if cpu_secs < 1: return 12
        if cpu_secs < 10: return 8
        return 5

    def main(self):
        a = self.a
        lim_def, lim_alpha = a.limit_default, a.limit_alpha
        target = a.target_secs
        # 1. Calibration with the default alpha, all threads
        wall = {}
        for n in range(a.min_n, MAX_N[a.algo] + 1):
            if n - 1 in wall and n - 2 in wall:
                growth = max(wall[n - 1] / max(wall[n - 2], 1e-3), 3.0)
                if wall[n - 1] * growth > 0.9 * lim_def:
                    break
            w = self.run(n, "def", "cal", 1, lim_def)
            if w is None: break
            wall[n] = w
            if w >= target:
                # long enough: the calibration run is regular MT data.
                # Record it under its new tag too, e.g. for replays.
                src = os.path.join(a.out, self.tag(n, "def", "cal", 1) + ".log")
                if os.path.exists(src):
                    os.replace(src, os.path.join(a.out, self.tag(n, "def", "mt", 1) + ".log"))
                    rec = dict(self.last_rec, tag=self.tag(n, "def", "mt", 1), mode="mt",
                               promoted_from=self.last_rec["tag"])
                    with open(self.index, "a") as f:
                        f.write(json.dumps(rec) + "\n")
        print("calibration (default alpha, all threads):", {n: round(w, 2) for n, w in wall.items()}, flush=True)

        # Thread count per size: all threads if the run takes at least
        # target seconds, else fewer threads so that a run takes about
        # target seconds; ST if even 2 threads would be too fast.
        full = [n for n, w in wall.items() if w >= target]
        small = sorted(n for n, w in wall.items() if w < target)
        for n in full:
            self.threads[n] = self.full_threads
        st_ns = []
        for n in small:
            k = int(self.full_threads * wall[n] / target)
            if k >= 2: self.threads[n] = min(k, self.full_threads)
            else: st_ns.append(n)
        def_wall = {n: wall[n] for n in full}

        # 2. Default alpha, MT repetitions: all-thread sizes largest
        # first, then the reduced-thread sizes
        for n in sorted(full, reverse=True):
            for rep in range(2, self.mt_reps(wall[n]) + 1):
                if self.remaining() < wall[n] * 1.5: break
                self.run(n, "def", "mt", rep, lim_def)
        for n in sorted(set(self.threads) - set(full), reverse=True):
            k = self.threads[n]
            w = self.run(n, "def", "mt", 1, lim_def, threads=k)
            if w is not None and w < 0.6 * target and k > 2:
                # thread startup made the calibration run look longer
                k = max(2, int(k * w / target))
                self.threads[n] = k
                w = self.run(n, "def", "mt", 1, lim_def, threads=k)
            if w is None: continue
            def_wall[n] = w
            print("1e%d: %d threads, %.2fs per run" % (n, k, w), flush=True)
            for rep in range(2, self.mt_reps(w) + 1):
                if self.remaining() < w * 1.5: break
                self.run(n, "def", "mt", rep, lim_def, threads=k)
        mt_ns = sorted(def_wall)

        # 3. ST sizes (and the smallest MT size for the ST/MT check).
        # Calibrate the ST runtime with one pinned run per size, then
        # run the repetitions and the non-default alphas in parallel
        # on distinct physical cores.
        check = [min(mt_ns)] if mt_ns else []
        st_secs = {}
        for n in st_ns + check:
            w = self.run(n, "def", "st", 1, a.st_max_secs, cpu=self.cpus[0])
            if w is None or w > a.st_max_secs: break
            st_secs[n] = w
        jobs = []
        for n, w in st_secs.items():
            for rep in range(2, self.st_reps(w) + 1):
                jobs.append((n, "def", rep))
            if n in check: continue
            for f in ALPHA_FACTORS[a.algo]:
                alpha = self.alpha(n, f)
                if alpha is None: continue
                for rep in range(1, max(2, self.st_reps(w) // 3) + 1):
                    jobs.append((n, alpha, rep))
        self.run_parallel(jobs, a.st_max_secs * 5)

        # 4. Non-default alphas, MT, with the same thread count as the
        # default alpha runs of that size. Measure each alpha's runtime
        # ratio at the smallest MT size, then run the remaining sizes
        # largest first, skipping runs predicted to exceed the limit.
        if mt_ns:
            n0 = mt_ns[0]
            ratio = {}
            for f in ALPHA_FACTORS[a.algo]:
                alpha = self.alpha(n0, f)
                if alpha is None: continue
                w = self.run(n0, alpha, "mt", 1, lim_alpha, threads=self.threads[n0])
                if w is not None: ratio[f] = w / def_wall[n0]
            for n in reversed(mt_ns):
                for f, r in ratio.items():
                    alpha = self.alpha(n, f)
                    if alpha is None: continue
                    pred = def_wall[n] * r * 1.2
                    if pred > 0.9 * lim_alpha: continue
                    reps = max(1, self.mt_reps(pred) // 2)
                    for rep in range(2 if n == n0 else 1, reps + 1):
                        if self.remaining() < pred * 1.5: break
                        self.run(n, alpha, "mt", rep, lim_alpha, threads=self.threads[n])
        print("threads per size:", self.threads, "ST sizes:", sorted(st_ns), flush=True)
        print("done, elapsed %.1f min" % ((time.time() - self.t0) / 60), flush=True)

    def alpha(self, n, factor):
        x = 10 ** n
        lo, hi = 1.0, float(iroot(x, 6))
        alpha = min(max(default_alpha(self.a.algo, x) * factor, lo), hi)
        alpha = trunc3(alpha)
        if abs(alpha / default_alpha(self.a.algo, x) - 1) < 0.05:
            return None   # too close to the default alpha
        return alpha

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--algo", required=True, choices=list(ALGO_OPT))
    p.add_argument("--binary", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--budget-min", type=float, default=110)
    p.add_argument("--threads", type=int, default=0, help="MT threads (default: primecount's default)")
    p.add_argument("--no-smt", action="store_true",
                   help="pin MT runs to one hardware thread per physical core and use that many threads")
    p.add_argument("--min-n", type=int, default=12)
    p.add_argument("--limit-default", type=float, default=600)
    p.add_argument("--limit-alpha", type=float, default=300)
    p.add_argument("--target-secs", "--mt-min-secs", type=float, default=5.0,
                   help="minimum MT run time; shorter sizes use fewer threads")
    p.add_argument("--st-max-secs", type=float, default=120)
    p.add_argument("--st-parallel", type=int, default=48)
    p.add_argument("--rep-scale", type=float, default=1.0, help="multiply all repetition counts")
    Collector(p.parse_args()).main()
```

## Appendix B: fit.py

<!-- file: fit.py -->
```python
#!/usr/bin/env python3
"""
Fit the WorkModelS2 progress model of one algorithm (D, S2-hard or
lmo) to the logs written by collect.py.

Ground truth: for each run the CPU time of every LoadBalancerS2 chunk
[low, high[ is known, hence F(pos) = CPU time below pos / total CPU
time is the true progress at position pos, independent of the thread
count (verified: 1 vs 4 threads differ less than run-to-run noise).
Repetitions are combined using the per-position median.

Model (StatusS2 with WorkModelS2, positions relative to leaves_start s):
  work(low)  = sieve * low / L + leaves(low)
  sieve      = sieve_work * (default_alpha / alpha)^alpha_exponent * (L / s)^sieve_exponent
  leaves     = 1 - (1 + (low - s) / (leaves_scale * s))^-leaves_shape   for low > s
  percent    = 100 * work(low) / work(L)
where L = sieve_limit and s = in_between(1, leaves_start, L) comes from
the algorithm's leaves_start formula with fitted multipliers. In C++:
  WorkModelS2 model(leaves_start, sieve_work * (default_alpha / alpha)^alpha_exponent,
                    leaves_scale, leaves_shape, sieve_exponent);

Objective per configuration: work-weighted mean |model - truth| in
percentage points + rem_weight * work-weighted mean |log2(remaining
time ratio)| for 10% <= truth <= 99%, as used by the load balancer.
Configurations are weighted by x (larger x more important) and the
default alpha gets a higher weight.

Checks printed before the fit: incomplete logs (chunks must cover
[0, sieve_limit[ without gaps, match the logged chunk count and have
valid CPU times), replica checks (the logged WorkModelS2 values vs
CURRENT, the logged getPercent() vs StatusS2's formula), the bimodality
check (runs with > 5% more CPU time than the fastest run of their
configuration) and the ST vs MT curves at the overlap size.
"""
import argparse, glob, json, math, os, sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect import default_alpha, iroot

PARAMS = {
    "D": ["start_x^(2/5)", "start_x/y^(4/3)", "start_x/y^(3/2)",
          "leaves_scale", "leaves_shape", "sieve_work", "sieve_exponent", "alpha_exponent"],
    "S2-hard": ["start_x/y^(3/2)", "start_y",
                "leaves_scale", "leaves_shape", "sieve_work", "sieve_exponent", "alpha_exponent"],
    "lmo": ["start_x/y^2",
            "leaves_scale", "leaves_shape", "sieve_work", "sieve_exponent", "alpha_exponent"],
}
# The constants in primecount's source code (D.cpp, S2_hard.cpp,
# pi_lmo_parallel.cpp), keep in sync. The replica checks compare them
# with the values the instrumented binary logged.
CURRENT = {
    "D": [1.6473, 1.1504, 2.0473, 1.9307, 0.41271, 0.0149, 0.34966, 1.082],
    "S2-hard": [1.4064, 0.45799, 4.3158, 0.61257, 0.076243, 0.1403, 0.9068],
    "lmo": [1.3352, 1.1764, 0.36355, 0.47309, 0.030408, 2.5849],
}

def leaves_start(algo, p, x, y):
    if algo == "D":
        s = min(p[0] * x ** 0.4, p[1] * x / y ** (4 / 3))
        return max(s, p[2] * x / y ** 1.5)
    if algo == "S2-hard":
        return min(p[0] * x / y ** 1.5, p[1] * y)
    return p[0] * x / (y * y)

def read_blocks(path):
    blocks, cur = [], None
    for line in open(path):
        if line.startswith("#"):
            t = line[1:].split()
            cur = {"meta": {t[i]: float(t[i + 1]) for i in range(0, len(t) - 1, 2)}, "rows": []}
            blocks.append(cur)
        elif cur is not None and line.strip():
            cur["rows"].append([float(v) for v in line.split()])
    return blocks

class Config:
    """All repetitions of one (n, alpha, mode)."""
    def __init__(self, algo, key, runs, ngrid):
        self.algo, self.key = algo, key
        n, alpha, mode = key
        self.n, self.mode, self.is_default = n, mode, alpha == "def"
        m = runs[0]["meta"]
        self.meta = m
        self.x, self.y, self.L = m["x"], m["y"], m["sieve_limit"]
        self.alpha = self.y / iroot(10 ** n, 3)
        self.def_alpha = default_alpha(algo, 10 ** n)
        self.G = np.geomspace(max(self.L * 1e-7, 1.0), self.L, ngrid)
        curves = []
        for r in runs:
            rows = np.array(sorted(r["rows"]))
            hi = np.minimum(rows[:, 1], self.L)
            w = rows[:, 3]
            xs = np.concatenate([[rows[0, 0]], hi])
            ys = np.concatenate([[0.0], np.cumsum(w)]) / w.sum()
            curves.append(np.interp(self.G, xs, ys))
        C = 100 * np.array(curves)
        self.C = C
        # total CPU time of each run (bimodality check)
        self.cpu = np.array([sum(row[3] for row in r["rows"]) for r in runs])
        self.P = np.median(C, axis=0)
        self.reps = len(runs)
        # noise: median absolute deviation of single runs (pp, max over grid)
        self.noise = float(np.max(np.median(np.abs(C - self.P), axis=0))) if len(runs) > 1 else float("nan")
        self.dP = np.diff(np.concatenate([[0.0], self.P]))
        self.rem_mask = (self.P >= 10) & (self.P <= 99)
        self.runs = runs

    def percent(self, p):
        algo = self.algo
        k = {"D": 3, "S2-hard": 2, "lmo": 1}[algo]
        scale, shape, sw, se, ae = p[k:k + 5]
        s = min(max(leaves_start(algo, p, self.x, self.y), 1.0), self.L)
        sieve = sw * (self.def_alpha / self.alpha) ** ae * (self.L / s) ** se
        def work(low):
            low = np.asarray(low, float)
            u = np.maximum(low - s, 0) / (scale * s)
            return sieve * low / self.L + 1 - (1 + u) ** (-shape)
        return np.clip(100 * work(self.G) / work(self.L), 0, 100), work

    def errors(self, p):
        Pm, _ = self.percent(p)
        d = np.abs(Pm - self.P)
        mean = float((d * self.dP).sum() / self.dP.sum())
        m = self.rem_mask
        if m.any():
            pm = np.clip(Pm[m], 10, 99.99)
            r = np.abs(np.log2((100 / pm - 1) / (100 / self.P[m] - 1)))
            rem = float((r * self.dP[m]).sum() / self.dP[m].sum())
        else:
            rem = 0.0
        return float(d.max()), mean, rem

    def wall_error(self, p):
        """Wall-clock check for MT runs: model percent at each assigned
        chunk's low vs elapsed / total wall time (includes load
        balancer lag and timing noise)."""
        if self.mode != "mt": return float("nan")
        _, work = self.percent(p)
        tot = work(self.L)
        errs = []
        for r in self.runs:
            rows = np.array(sorted(r["rows"], key=lambda v: v[4]))
            T = r["meta"]["wall"]
            t = rows[:, 4]
            Pm = np.clip(100 * work(rows[:, 0]) / tot, 0, 100)
            dt = np.diff(np.concatenate([t, [T]]))
            errs.append(float((np.abs(Pm - 100 * t / T) * dt).sum() / dt.sum()))
        return float(np.median(errs))

    def wall_lag(self):
        """Signed lag for MT runs: true progress at each assigned chunk's
        low minus elapsed / total wall time, time-weighted mean (pp).
        Independent of the model. Positive = the assigned position is
        ahead of the finished work (in-flight work). Returns the median
        and the min/max over the repetitions."""
        if self.mode != "mt": return float("nan"), float("nan"), float("nan")
        lags = []
        for r in self.runs:
            rows = np.array(sorted(r["rows"]))
            hi = np.minimum(rows[:, 1], self.L)
            xs = np.concatenate([[rows[0, 0]], hi])
            ys = np.concatenate([[0.0], np.cumsum(rows[:, 3])]) / rows[:, 3].sum()
            rows = rows[np.argsort(rows[:, 4])]
            T = r["meta"]["wall"]
            t = rows[:, 4]
            Pf = 100 * np.interp(rows[:, 0], xs, ys)
            dt = np.diff(np.concatenate([t, [T]]))
            lags.append(float(((Pf - 100 * t / T) * dt).sum() / dt.sum()))
        return float(np.median(lags)), min(lags), max(lags)

def check_block(b):
    """None if the block's chunks cover [0, sieve_limit[ completely,
    else the reason. Replica checks can't detect incomplete logs."""
    m, rows = b["meta"], sorted(b["rows"])
    if "chunks" in m and len(rows) != int(m["chunks"]):
        return "%d of %d chunks" % (len(rows), int(m["chunks"]))
    if rows[0][0] != 0:
        return "first chunk starts at %d" % rows[0][0]
    for prev, row in zip(rows, rows[1:]):
        if row[0] != prev[1]:
            return "gap or overlap at %d" % row[0]
    if rows[-1][1] < m["sieve_limit"]:
        return "chunks end at %d < sieve_limit %d" % (rows[-1][1], m["sieve_limit"])
    times = [v for row in rows for v in row[2:4]]
    if not all(math.isfinite(v) and v >= 0 for v in times):
        return "invalid CPU time"
    if sum(row[3] for row in rows) <= 0:
        return "total CPU time 0"
    return None

def load(algo, paths, ngrid, min_n, mt_min_wall):
    groups = {}
    files = [f for p in paths for f in sorted(glob.glob(os.path.join(p, "%s_1e*_r*.log" % algo)))]
    for f in files:
        name = os.path.basename(f)[:-4]
        _, ns, alpha, mode, rep = name.split("_")
        n = int(ns[2:])
        # all-thread calibration runs shorter than collect.py's target
        if n < min_n or mode == "cal": continue
        x = 10.0 ** n
        # LMO computes nested pi(x) using other LoadBalancerS2
        # instances, select the block of the requested x
        bl = [b for b in read_blocks(f) if abs(b["meta"]["x"] / x - 1) < 1e-9 and b["rows"]]
        if not bl: continue
        reason = check_block(bl[-1])
        if reason:
            print("warning: incomplete log %s: %s, not used" % (f, reason))
            continue
        groups.setdefault((n, alpha, mode), []).append(bl[-1])
    configs = [Config(algo, k, v, ngrid) for k, v in sorted(groups.items())]
    # MT runs shorter than mt_min_wall are dominated by thread
    # startup and timing noise, they are only used for the ST/MT check
    for c in configs:
        c.fit = c.mode == "st" or np.median([r["meta"]["wall"] for r in c.runs]) >= mt_min_wall
    return configs

def model_work(L, start, sieve_work, scale, shape, exponent):
    """StatusS2::work() for the WorkModelS2 values passed by an algorithm"""
    s = min(max(start, 1.0), L)
    sieve = sieve_work * (L / s) ** exponent
    def work(low):
        low = np.asarray(low, float)
        u = np.maximum(low - s, 0) / (scale * s)
        return sieve * low / L + 1 - (1 + u) ** (-shape)
    return work

def replica_checks(algo, configs, cur):
    """Compare the WorkModelS2 values logged by the instrumented binary
    with this file's replica (CURRENT, default_alpha) and the logged
    getPercent() of every chunk with StatusS2's formula."""
    k = {"D": 3, "S2-hard": 2, "lmo": 1}[algo]
    scale, shape, sw, se, ae = cur[k:k + 5]
    keys = ["leaves_start", "sieve_work", "leaves_scale", "leaves_shape", "sieve_exponent"]
    bad, worst, checked = {}, 0.0, 0
    old = [c for c in configs if "sieve_work" not in c.meta]
    if old:
        print("note: %d configs were logged by an older logging patch, only leaves_start is checked" % len(old))
    for c in configs:
        m = c.meta
        exp = dict(zip(keys, [leaves_start(algo, cur, c.x, c.y),
                              sw * (c.def_alpha / c.alpha) ** ae, scale, shape, se]))
        for key in keys:
            if key in m and abs(exp[key] / m[key] - 1) > 1e-6:
                bad.setdefault(key, []).append("1e%d %s %s: replica %.10g, logged %.10g"
                                               % (c.key + (exp[key], m[key])))
        if all(key in m for key in keys):
            checked += 1
            work = model_work(c.L, *[m[key] for key in keys])
            tot = work(c.L)
            for r in c.runs:
                rows = np.array(r["rows"])
                pct = np.clip(100 * work(rows[:, 0]) / tot, 0, 100)
                worst = max(worst, float(np.abs(pct - rows[:, 6]).max()))
    for key in keys:
        if key in bad:
            print("warning: %s replica mismatch in %d configs, e.g. %s" % (key, len(bad[key]), bad[key][0]))
    if worst > 1e-6:
        print("warning: getPercent() replica mismatch, max %.3g pp" % worst)
    if not bad and worst <= 1e-6:
        print("replica checks OK (%d configs), getPercent() %s" % (len(configs),
              "max diff %.2g pp" % worst if checked else "not checked (older logging patch)"))

def bimodality_checks(configs):
    """Runs with > 5% more CPU time than the fastest run of their
    configuration were disturbed (e.g. NUMA). Report how much the
    median curve moves without them: a single outlier is harmless."""
    found = False
    for c in configs:
        slow = c.cpu > 1.05 * c.cpu.min()
        if c.reps < 2 or not slow.any(): continue
        found = True
        shift = float(np.abs(np.median(c.C[~slow], axis=0) - c.P).max())
        print("%s: 1e%d %s %s: %d/%d runs use > 5%% more CPU time than the fastest run (max +%.0f%%), "
              "without them the median curve moves %.2f pp"
              % (("warning" if shift >= 0.2 else "note",) + c.key
                 + (slow.sum(), c.reps, 100 * (c.cpu.max() / c.cpu.min() - 1), shift)))
    if not found:
        print("bimodality check OK: no run uses > 5% more CPU time than the fastest run of its configuration")

def max_error(c, p):
    """Largest model error (signed, pp) and the true progress (%) where
    it occurs, overall and for true progress >= 10%: below 10% the
    error only affects the status, LoadBalancerS2::remaining_secs()
    clamps the percentage to >= 10."""
    Pm, _ = c.percent(p)
    d = Pm - c.P
    i = int(np.argmax(np.abs(d)))
    j = int(np.argmax(np.where(c.P >= 10, np.abs(d), -1)))
    return d[i], c.P[i], d[j], c.P[j]

def cpp_code(algo, p):
    """The fitted constants in the form used by primecount's source"""
    f = lambda v: "%.5g" % v
    k = {"D": 3, "S2-hard": 2, "lmo": 1}[algo]
    scale, shape, sw, se, ae = p[k:k + 5]
    if algo == "D":
        start = "max(min(%s * x^(2/5), %s * x / y^(4/3)), %s * x / y^(3/2))" % (f(p[0]), f(p[1]), f(p[2]))
    elif algo == "S2-hard":
        start = "min(%s * x / y^(3/2), %s * y)" % (f(p[0]), f(p[1]))
    else:
        start = "%s * x / y^2" % f(p[0])
    return ["leaves_start = " + start,
            "sieve_work = %s * pow(default_alpha / alpha, %s)" % (f(sw), f(ae)),
            "WorkModelS2 model(leaves_start, sieve_work, %s, %s, %s);" % (f(scale), f(shape), f(se))]

def weights(configs, xpow, def_weight):
    # each (n, default?) group shares the weight (n - 10)^xpow
    grp = {}
    for c in configs:
        grp.setdefault((c.n, c.is_default), []).append(c)
    w = {}
    for (n, d), cs in grp.items():
        for c in cs:
            w[c.key] = max(n - 10, 1) ** xpow * (def_weight if d else 1) / len(cs)
    return w

def objective(p, configs, w, rem_weight):
    if np.any(p <= 0) or not np.all(np.isfinite(p)): return 1e9
    tot = 0.0
    # Extreme trial parameters may overflow, treat them as very bad
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for c in configs:
            _, mean, rem = c.errors(p)
            tot += w[c.key] * (mean + rem_weight * rem)
    tot /= sum(w[c.key] for c in configs)
    return tot if np.isfinite(tot) else 1e9

def fit(configs, w, rem_weight, x0, starts, seed=1):
    rng = np.random.default_rng(seed)
    f = lambda q: objective(np.exp(q), configs, w, rem_weight)
    best = None
    inits = [np.log(x0)] + [np.log(x0) + rng.normal(0, 0.4, len(x0)) for _ in range(starts)]
    for q0 in inits:
        r = minimize(f, q0, method="Nelder-Mead",
                     options=dict(maxiter=6000, maxfev=12000, xatol=1e-4, fatol=1e-5, adaptive=True))
        r = minimize(f, r.x, method="Powell", options=dict(maxiter=20000, xtol=1e-4, ftol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return np.exp(best.x), best.fun

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--algo", required=True, choices=list(PARAMS))
    ap.add_argument("--data", required=True, nargs="+", help="one or more collect.py output directories")
    ap.add_argument("--min-n", type=int, default=12)
    ap.add_argument("--grid", type=int, default=400)
    ap.add_argument("--mt-min-wall", type=float, default=2.0)
    ap.add_argument("--x-weight-power", type=float, default=1.0)
    ap.add_argument("--default-weight", type=float, default=3.0)
    ap.add_argument("--rem-weight", type=float, default=10.0)
    ap.add_argument("--starts", type=int, default=6)
    ap.add_argument("--no-cv", action="store_true")
    ap.add_argument("--check-only", action="store_true", help="print the data checks and exit")
    ap.add_argument("--current", default=None, metavar="FILE[:KEY]",
                    help="use the constants of a fit.py output JSON instead of CURRENT "
                         "(KEY: current or params, default current), e.g. for logs of an older binary")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    configs = load(a.algo, a.data, a.grid, a.min_n, a.mt_min_wall)
    if not configs: sys.exit("no logs found")
    cur = np.array(CURRENT[a.algo])
    if a.current:
        path, _, key = a.current.partition(":")
        d = json.load(open(path))[key or "current"]
        cur = np.array([d[k] for k in PARAMS[a.algo]])
    # data summary and consistency checks against the logged values
    for c in configs:
        if c.is_default:
            print("1e%d %s %s: %d reps, noise %.2f pp, wall %.2f s%s"
                  % (c.key + (c.reps, c.noise, np.median([r["meta"]["wall"] for r in c.runs]),
                              "" if c.fit else " (too short, not fitted)")))
    replica_checks(a.algo, configs, cur)
    bimodality_checks(configs)
    # ST/MT consistency: both must measure the same work distribution
    by = {c.key: c for c in configs}
    for c in configs:
        if c.mode == "st" and (c.n, c.key[1], "mt") in by:
            m = by[(c.n, c.key[1], "mt")]
            print("ST vs MT median curves at 1e%d %s: max diff %.2f pp (noise ST %.2f, MT %.2f)"
                  % (c.n, c.key[1], float(np.abs(c.P - m.P).max()), c.noise, m.noise))
    if a.check_only:
        return
    configs = [c for c in configs if c.fit]
    w = weights(configs, a.x_weight_power, a.default_weight)
    p, fun = fit(configs, w, a.rem_weight, cur, a.starts)
    print("objective: current %.4f, fitted %.4f" % (objective(cur, configs, w, a.rem_weight), fun))
    for name, v0, v in zip(PARAMS[a.algo], cur, p):
        print("  %-18s current %9.5g   fitted %9.5g" % (name, v0, v))

    # Cross-validation: leave one decade out, and extrapolation
    # to the two largest decades from the smaller ones
    ns = sorted(set(c.n for c in configs))
    lono = {}
    if not a.no_cv and len(ns) > 2:
        for n in ns:
            train = [c for c in configs if c.n != n]
            # Start from CURRENT, not from the full fit p: the held-out
            # decade must not influence the fold through the start point
            q, _ = fit(train, weights(train, a.x_weight_power, a.default_weight), a.rem_weight, cur, a.starts)
            for c in configs:
                if c.n == n: lono[c.key] = q
        top = ns[-2:]
        train = [c for c in configs if c.n not in top]
        qx, _ = fit(train, weights(train, a.x_weight_power, a.default_weight), a.rem_weight, cur, a.starts)
    else:
        qx = None

    print()
    print("errors: max/mean pp, remaining-time log2 error; wall = MT wall-clock check (mean pp)")
    print("%-22s %4s %6s | %-19s | %-19s | %-19s | %-19s | %-11s"
          % ("config", "reps", "noise", "current", "fitted", "leave-decade-out", "extrapolated", "wall cur/fit")
          + " | lag med [min, max]")
    fmt = lambda e: "%5.2f %5.2f %6.3f" % e
    for c in configs:
        lo = fmt(c.errors(lono[c.key])) if c.key in lono else "-"
        ex = fmt(c.errors(qx)) if qx is not None and c.n in ns[-2:] else "-"
        print("%-22s %4d %6.2f | %-19s | %-19s | %-19s | %-19s | %5.2f %5.2f"
              % ("1e%d %s %s" % c.key, c.reps, c.noise, fmt(c.errors(cur)), fmt(c.errors(p)), lo, ex,
                 c.wall_error(cur), c.wall_error(p))
              + " | %5.2f [%5.2f, %5.2f]" % c.wall_lag())

    print()
    print("max error at the default alpha (signed pp at true progress %, and for true progress >= 10%)")
    sets = [("current", cur), ("fitted", p)] + ([("extrapolated", qx)] if qx is not None else [])
    for c in configs:
        if c.is_default and c.n in ns[-2:]:
            print("  %-12s %s" % ("1e%d %s %s" % c.key, "  ".join(
                  "%s %+5.2f at %4.1f%% (>= 10%%: %+5.2f at %4.1f%%)" % ((name,) + tuple(max_error(c, q)))
                  for name, q in sets)))

    print()
    print("C++ (5 significant digits):")
    for line in cpp_code(a.algo, p):
        print("  " + line)

    summary = {}
    for c in configs:
        if c.is_default:
            e = dict(current=c.errors(cur), fitted=c.errors(p))
            if c.key in lono: e["leave_decade_out"] = c.errors(lono[c.key])
            if qx is not None and c.n in ns[-2:]: e["extrapolated"] = c.errors(qx)
            summary["1e%d %s" % (c.n, c.mode)] = {k: [round(v, 4) for v in e[k]] for k in e}
    res = dict(algo=a.algo, params=dict(zip(PARAMS[a.algo], map(float, p))),
               current=dict(zip(PARAMS[a.algo], map(float, cur))),
               objective_current=objective(cur, configs, w, a.rem_weight), objective_fitted=fun,
               extrapolation_params=None if qx is None else dict(zip(PARAMS[a.algo], map(float, qx))),
               default_alpha_errors_max_mean_rem=summary, cpp=cpp_code(a.algo, p))
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)
        print("\nwrote", a.out)

if __name__ == "__main__":
    main()
```

## Appendix C: logging patch

Applies to primecount with `WorkModelS2` (commit "Tune StatusS2 work
model per algorithm" and later), from the repository root with
`git apply primecount-s2-tuning-log.patch`. Never commit it.

<!-- file: primecount-s2-tuning-log.patch -->
```diff
diff --git a/src/LoadBalancerS2.cpp b/src/LoadBalancerS2.cpp
index 01a72326..055ec0ef 100644
--- a/src/LoadBalancerS2.cpp
+++ b/src/LoadBalancerS2.cpp
@@ -44,6 +44,12 @@
 #include <stdint.h>
 #include <atomic>
 #include <cmath>
+#include <cstdio>
+#include <cstdlib>
+
+#if defined(_OPENMP)
+  #include <omp.h>
+#endif

 namespace {

@@ -107,6 +113,44 @@ LoadBalancerS2::LoadBalancerS2(maxint_t x,

   if (is_print_)
     print_S2_status(start_time_);
+
+  log_path_ = std::getenv("PRIMECOUNT_S2_LOG");
+  log_x_ = double(x);
+  log_model_ = model;
+  x15_ = max(int64_t(std::pow(double(x), 1.0 / 5.0)), 1 << 10);
+
+  // Fixed chunks of chunk_frac * low (single-threaded measurements only)
+  if (const char* frac = std::getenv("PRIMECOUNT_S2_CHUNK_FRAC"))
+    if (threads == 1)
+      chunk_frac_ = std::atof(frac);
+}
+
+/// Append one block per load balancer to PRIMECOUNT_S2_LOG.
+/// Nested pi(x) computations (e.g. in LMO) write their own,
+/// smaller-x blocks, select the block by its x value.
+///
+LoadBalancerS2::~LoadBalancerS2()
+{
+  if (!log_path_)
+    return;
+
+  FILE* f = std::fopen(log_path_, "a");
+  if (!f)
+    return;
+
+  // WorkModelS2 as passed by the algorithm, sieve_work includes the alpha factor
+  const WorkModelS2& m = log_model_;
+  std::fprintf(f, "# x %.17g y %lld sieve_limit %lld leaves_start %.17g sieve_work %.17g leaves_scale %.17g "
+                  "leaves_shape %.17g sieve_exponent %.17g threads %d wall %.9f chunks %zu\n",
+               log_x_, (long long) y_, (long long) sieve_limit_, m.leaves_start, m.sieve_work, m.leaves_scale,
+               m.leaves_shape, m.sieve_exponent, threads_, get_time() - start_time_, log_chunks_.size());
+
+  // low high init_secs secs start stop percent thread
+  for (const LogChunk& c : log_chunks_)
+    std::fprintf(f, "%lld %lld %.9f %.9f %.9f %.9f %.12g %d\n", (long long) c.low, (long long) c.high,
+                 c.init, c.secs, c.start, c.stop, c.pct, c.tid);
+
+  std::fclose(f);
 }

 /// Remaining seconds till finished
@@ -138,6 +182,35 @@ void LoadBalancerS2::store_packed(uint64_t segment_size, uint64_t segments)
 ///
 bool LoadBalancerS2::get_work(ThreadData& thread)
 {
+  if (log_path_ &&
+      thread.segments != 0)
+  {
+    int tid = 0;
+    #if defined(_OPENMP)
+      tid = omp_get_thread_num();
+    #endif
+    int64_t high = min(thread.low + thread.segment_size * thread.segments, sieve_limit_);
+    LogChunk c = { thread.low, high, thread.init_secs(), thread.secs(),
+                   thread.start_time - start_time_, thread.stop_time - start_time_,
+                   status_.getPercent(thread.low), tid };
+    #pragma omp critical (s2_tuning_log)
+    log_chunks_.push_back(c);
+  }
+
+  if (chunk_frac_ > 0)
+  {
+    int64_t low = low_.load(std::memory_order_relaxed);
+    int64_t segment_size = in_between(x15_, low / 16, L1_segment_size);
+    segment_size = Sieve::align_segment_size(segment_size);
+    int64_t segments = (int64_t) std::llround(chunk_frac_ * double(low) / double(segment_size));
+    segments = max(segments, 1);
+    thread.low = low_.fetch_add(segment_size * segments, std::memory_order_relaxed);
+    thread.segment_size = segment_size;
+    thread.segments = segments;
+    thread.init_time = 0;
+    return thread.low < sieve_limit_;
+  }
+
   if (is_print_ &&
       thread.segments != 0)
     print_S2_status(thread.stop_time);
diff --git a/src/LoadBalancerS2.hpp b/src/LoadBalancerS2.hpp
index ab9372e4..a4ddb164 100644
--- a/src/LoadBalancerS2.hpp
+++ b/src/LoadBalancerS2.hpp
@@ -18,6 +18,7 @@
 #include <stdint.h>
 #include <algorithm>
 #include <atomic>
+#include <vector>

 namespace primecount {

@@ -49,6 +50,7 @@ class LoadBalancerS2
 public:
   LoadBalancerS2(maxint_t x, int64_t y, int64_t sieve_limit, const WorkModelS2& model, int threads, bool is_print);
   bool get_work(ThreadData& thread);
+  ~LoadBalancerS2();

 private:
   double remaining_secs(const ThreadData& thread, int64_t low) const;
@@ -65,6 +67,15 @@ private:
   bool is_print_ = false;
   StatusS2 status_;

+  // Tuning instrumentation, enabled by PRIMECOUNT_S2_LOG=<file>
+  struct LogChunk { int64_t low, high; double init, secs, start, stop, pct; int tid; };
+  std::vector<LogChunk> log_chunks_;
+  const char* log_path_ = nullptr;
+  double log_x_ = 0;
+  WorkModelS2 log_model_{0, 0, 0, 0, 0};
+  double chunk_frac_ = 0;
+  int64_t x15_ = 0;
+
   MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
   std::atomic<int64_t> low_{0};
   MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
```
