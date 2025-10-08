# SVAR-SVSR Experimental Suite

This repository hosts the MATLAB and Python code used to compare several
structural VAR identification strategies under sign, narrative, and
zero-restriction schemes. The workflow produces impulse-response
functions (IRFs) from three benchmark data sets and then visualises the
results in the accompanying Jupyter notebook.

## Repository Layout

- `run_svar_compare_*.m` – entry-point MATLAB scripts for each data set /
  identification scheme (Rubio-Ramírez narrative restrictions, Kilian
  oil-market model, and sign/zero restrictions).
- `functions/` – MATLAB helpers called by the driver scripts. Only the
  routines still referenced by the project remain in this folder.
- `output/` – **initially empty** destination for `.mat` files produced by the MATLAB runs.
  Sub-folders (`killian/`, `SignZero/`, …) mirror the experiment type. You must run the main scripts with appropriate parameter settings to generate these files (see "Parameter Configuration" below).
- `draw_irfs.ipynb` – Python notebook that loads the `.mat` files,
  computes summary statistics, and plots the IRFs with Plotly.

The repository assumes the raw data live in `data/` relative to the
project root:

| Script                              | Expected data file                            |
| ----------------------------------- | --------------------------------------------- |
| `run_svar_compare_rubioramirez.m`   | `data/Uhlig_Data_Updated.mat`           |
| `run_svar_compare_kilian.m`         | `data/Kilian_Data_Updated.mat`          |
| `run_svar_compare_signzero.m`       | `data/dataset_RubioRamirez_SignZero.csv`|

## Prerequisites

### MATLAB

- MATLAB R2021a or later (tested on recent releases with Parallel
  Computing Toolbox optional).
- Ability to run scripts from the repository root (`svar-svsr`) so
  relative paths resolve correctly.
- Write permissions for the `output/` directory.

> **Note:** The scripts call `addpath(genpath('/functions'));`. When
> running from MATLAB on Windows or macOS you may prefer to replace this
> with `addpath(genpath(fullfile(pwd, 'functions')));` to avoid issues
> with absolute paths.

### Python

- Python 3.9+ with `pip`.
- Recommended packages: `numpy`, `scipy`, `h5py`, `plotly`, `nbformat`,
  and `jupyter` (or `jupyterlab`) for notebook execution.

You can install the Python dependencies into a virtual environment:

```bash
python -m venv .venv
. .venv/Scripts/activate            # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install numpy scipy h5py plotly jupyter nbformat
```

## Producing the IRFs with MATLAB

### Parameter Configuration for Different Methods

**Important:** The `output/` folder is initially empty. Each main script must be run **multiple times** with different parameter settings to generate the various output files. Open the script in MATLAB, modify the parameters at the top of the file, and run it.

#### For `run_svar_compare_rubioramirez.m`:

To generate different identification methods, modify these four key parameters while leaving all others unchanged:

| Output File | `regime_on` | `prior_on` | `svar_svsr_on` | `penalty_offdiagonal_on` |
|------------|-------------|------------|----------------|--------------------------|
| `output/irfs_narrative_rubioramirez.mat` | 0 | 0 | 0 | 0 |
| `output/irfs_narrative_dawis_regime.mat` | 1 | 0 | 1 | 0 |
| `output/irfs_narrative_dawis_regime_penalty.mat` | 1 | 0 | 1 | 1 |
| `output/irfs_narrative_carriero_regime.mat` | 1 | 0 | 0 | 0 |

#### For `run_svar_compare_kilian.m` and `run_svar_compare_signzero.m`:

Similar parameter patterns apply. Modify the same parameters to generate corresponding outputs in `output/killian/` and `output/SignZero/` subfolders.

### Running the Scripts

Each experiment follows this pattern:

1. Launch MATLAB and set the working directory to the repository root
   (`svar-svsr`).
2. Confirm the corresponding data file is available in `../../data/`.
3. Open the appropriate script in the MATLAB editor:
   - `run_svar_compare_rubioramirez.m`
   - `run_svar_compare_kilian.m`
   - `run_svar_compare_signzero.m`
4. **Modify the parameter settings** at the top of the script according to the configuration table above.
5. Run the script. It will save a `.mat` file under `output/…`.
6. **Repeat steps 4-5** with different parameter settings to generate all required output files.

> **Heads-up:** Two of the scripts still contain `system('shutdown …')`
> commands at the end. Comment these lines out if you do not want the
> machine to schedule a shutdown after completion.

If you prefer to run MATLAB non-interactively from the command line (after editing parameters):

```bash
matlab -batch "run('run_svar_compare_rubioramirez.m')"
```

## Recreating the Figures in `draw_irfs.ipynb`

Once the `.mat` files are available:

1. Activate your Python environment and launch Jupyter:

   ```bash
   jupyter notebook draw_irfs.ipynb
   ```

2. In the notebook, run all cells. The helper functions inside load the
   MATLAB outputs via `scipy.io.loadmat`/`h5py`, normalise the IRFs, and
   produce interactive Plotly figures for each identification comparison
   (narrative, elasticity, and sign-zero cases).

3. The notebook reads files from the paths listed above; ensure the cell
   inputs match the filenames produced by your MATLAB runs (e.g.
   `irfs_dawis_regime_penalty.mat`).

The resulting figures are displayed inline and can be exported from the
Plotly UI as static images or HTML for sharing.

## Suggested Workflow

1. Run the MATLAB script that matches the experiment you wish to
   reproduce **with each required parameter configuration**. Verify `.mat` files appear under `output/`.
2. Optionally run the other scripts if you plan to compare across all
   identification schemes.
3. Launch `draw_irfs.ipynb` and execute the relevant sections to create
   publication-ready plots.

Remember: The `output/` folder starts empty. You must generate all output files by running the main scripts multiple times with different parameter settings as documented in the "Parameter Configuration" section above.

Feel free to adapt the scripts or notebook functions—only the utilities
still referenced by the pipeline remain in `functions/`, so the project
is ready for custom extensions.
