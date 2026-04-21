# Reproducibility sweep — 2026-04-16

## What this folder is

An end-to-end sanity check that the Comment 1.9 paired-run results can be
reproduced from nothing but a fresh `git clone` of the public repository on
a clean build tree. No local state, no previously-built `.mod` files, no
reuse of the development worktrees.

Both the ion-on run (`newNetwork-transport-fixes`, commit
`73e8011704e0dda4fe0bb0f0e0075de475b3e80c`) and the ion-off control
(`newNetwork-reviewer-noion`, commit `a0bedbb643c4d08ef3d05b422871927403b8f9f6`)
reproduce **bit-exactly** against the reference artifacts in
`../comment_1.9_paired/`.

## Result

| Quantity                | Paired-run reference | Cold-clone replay |
|-------------------------|---------------------:|------------------:|
| Ion-on  unweighted RMSD |             2.011529 |      **2.011529** |
| Ion-on  weighted RMSD   |             2.429285 |      **2.429285** |
| Ion-on  peak bO3 %      |                26.80 |         **26.80** |
| Ion-off unweighted RMSD |            13.863138 |     **13.863138** |
| Ion-off weighted RMSD   |            16.460337 |     **16.460337** |
| Ion-off peak bO3 %      |                 3.66 |          **3.66** |

Additionally, `bO3.csv` from each cold-clone run is byte-identical to the
corresponding file archived under `../comment_1.9_paired/` (verified with
`diff`).

## What's in this folder

```
ion-on/
    bO3.csv                  Full bO3 trajectory vs fluence (99 rows).
    model_evaluation.png     Plot produced by singleRMSD.py.
    singleRMSD_stdout.txt    Full stdout of singleRMSD.py, per-point table + summary.
ion-off/
    bO3.csv                  As above, ion-off side.
    model_evaluation.png
    singleRMSD_stdout.txt
provenance.txt               Commit SHAs, toolchain, FLAGS, uname, date.
```

## How to reproduce (copy/paste)

```bash
# fresh clone, ion-on
git clone git@github.com:cnshingledecker/Dday.git ion-on
cd ion-on
git checkout newNetwork-transport-fixes
make clean
make FLAGS="-g -O2 -march=native -ffree-line-length-512 -J ."
bash run.sh
python3 singleRMSD.py      # -> Unweighted RMSD 2.011529, Weighted 2.429285

# fresh clone, ion-off
cd ..
git clone git@github.com:cnshingledecker/Dday.git ion-off
cd ion-off
git checkout newNetwork-reviewer-noion
make clean
make FLAGS="-g -O2 -march=native -ffree-line-length-512 -J ."
bash run.sh
python3 singleRMSD.py      # -> Unweighted RMSD 13.863138, Weighted 16.460337
```

The `FLAGS` override drops `-pg`, which on current macOS pulls in `gcrt1.o`
(not installed by default). No numerical effect. This is a pre-existing
Makefile / platform interaction, not introduced by any of the revision-era
commits.

## Interpretation

The Comment 1.9 quantitative claim — factor-of-~7 RMSD degradation and
peak O3 yield dropping to ~14% of the ion-on value when the twelve
charged-product photoprocess rows are muted — rests on numbers that are
now verifiably reproducible from the public remote, not just from the
development machine state. No hidden local edits, no unstaged files, no
orphaned executables are needed to get these numbers. The revisions can
be written with this provenance on the record.
