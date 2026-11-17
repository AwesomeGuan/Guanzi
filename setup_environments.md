# Guanzi environment setup

Run the installer from the repository root on Linux with Conda available. It
creates the six environments required by the pipeline and validates every
existing environment without updating it:

```bash
bash setup_environments.sh
```

Missing environments are created from the specifications below. Existing
environments are validated in place. Successful creation or validation writes
an environment completion marker; failure removes the corresponding marker and
stops setup with manual recovery instructions. Models and reference datasets
are installed separately through `setup_resources.md`.

Stage wrappers enter their required environment through
`scripts/run_in_env.sh`; manual activation is not required for normal pipeline
runs.

## 1. `defense_system`

This is the general CPU analysis environment. It supports Step1 table and
sequence processing in both modes, Step2 table filtering and import operations,
Step3 assembly analysis, the CPU portions of Step4 quality control, and Step5
information integration.

Main packages:

- Python 3.10
- pandas, Biopython, requests and tqdm
- DIAMOND
- CD-HIT
- TM-align

The installer validates the Python imports and executable availability. CD-HIT
is checked by resolving its executable because `cd-hit -h` can print valid help
text while returning a nonzero status.

## 2. `prodigal`

This environment is used by Custom Step1 to predict CDS features and translated
proteins from user-supplied genome or MAG assemblies.

Main packages:

- Python 3.10
- Prodigal

Python is required because the stage wrapper uses the shared Python logging
utilities while Prodigal is active. Validation therefore requires Python 3.10,
`time.monotonic()` and a working Prodigal executable.

## 3. `esmfold`

This GPU environment predicts protein structures in Custom Step1 and predicts
retained GTDB Mining-component structures in Replication Step4.

Main packages:

- Python 3.9
- PyTorch 1.12.1, torchvision 0.13.1 and torchaudio 0.12.1
- CUDA Toolkit 11.6
- Transformers 4.27.0, Accelerate, ModelScope and safetensors 0.3.3
- NumPy 1.24.3, pandas, SciPy, scikit-learn and Biopython

The host NVIDIA driver must support the installed CUDA runtime. ESMFold model
weights are not installed by this environment script; follow
`setup_resources.md` to deploy them.

## 4. `foldseek`

This environment builds and searches the Custom Step1 structure database and
runs the Replication Step1 structural search against AFDB50.

Main packages:

- Python 3.10
- FoldSeek
- pandas, requests and Biopython

AFDB50 and other structure databases are resources rather than environment
packages and are described in `setup_resources.md`.

## 5. `defensefinder`

This environment runs DefenseFinder annotation for Custom Step2 and supplies
the DefenseFinder model-update command used during resource setup. Replication
Step2 imports release-matched precomputed DefenseFinder results instead of
rerunning annotation across GTDB.

Main packages:

- Python 3.12 and pip
- mdmparis-defense-finder
- HMMER
- pandas and NumPy

DefenseFinder model definitions are downloaded separately as described in
`setup_resources.md`.

## 6. `eggnog`

This CPU environment runs the local eggNOG annotation operation in Step2 for
both modes. It is separate from `defense_system` and from database deployment.

Main packages:

- Python 3.10
- eggNOG-mapper 2.1.13
- DIAMOND 2.2.8
- Biopython, psutil and XlsxWriter, installed as tool dependencies
- SQLite support through Python's standard library

The installer verifies Python and package versions, `emapper.py --help`, and
DIAMOND. Software validation does not require a downloaded database. Follow
`setup_resources.md` to deploy eggNOG 5.0.2; the annotation wrapper verifies its
installed version before running. No MMseqs2, HMMER, novel-family or de novo
Pfam database is required for this DIAMOND annotation route.

## 7. Environment completion markers and recovery

Successful creation or validation writes these empty markers:

```text
output/.done/setup_envs/
├── 01_defense_system.done
├── 02_prodigal.done
├── 03_esmfold.done
├── 04_foldseek.done
├── 05_defensefinder.done
└── 06_eggnog.done
```

The markers record the latest successful environment check; they do not bypass
validation on later installer runs. Every invocation still validates all six
existing environments in order. A failed validation removes that environment's
marker and stops. A stale marker for a missing environment is removed before
the environment is recreated.

The installer never repairs, updates or removes an existing environment. To
replace an environment that fails validation, remove only that environment and
run the installer again:

```bash
conda env remove -n <environment> -y
bash setup_environments.sh
```

Do not remove an environment solely because another environment failed its
independent validation.
