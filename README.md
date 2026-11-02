# Piuma / MTS Data Analyzer

**Full-loading Hertz contact fitting with interactive review and optional R²-based Excel exclusion**

Companion script: `piuma_mts_analyzer_r2_cutoff.py`

## Overview

This Python program imports tab-delimited force–displacement `.txt` files, fits each measurement using a spherical Hertz contact model, and presents the results in a graphical interface. It supports manual review and exclusion of measurements, optional automatic exclusion based on fit R², and export of accepted results to an Excel spreadsheet.

**Important:** The R² cutoff controls **which results are exported**; it does **not** restrict which curves are fitted or change the fitted parameter estimates. Excluded measurements remain visible in the application's results table.

## 1. Requirements

- Python 3, with Tkinter GUI support (Tkinter is commonly included with standard Python installations on Windows).
- Python packages: `numpy`, `pandas`, `matplotlib`, `scipy`, and `openpyxl`.

Install the required packages from a terminal or Windows Command Prompt:

```bash
python -m pip install numpy pandas matplotlib scipy openpyxl
```

If Tkinter is unavailable on Linux, install your distribution's `python3-tk` package. On Windows, installing Python from python.org with Tcl/Tk support generally provides Tkinter.

## 2. Prepare the input files

Select the **top-level folder containing your measurement folders**, not an individual `.txt` file. The program scans each immediate subfolder of the selected folder and reads `.txt` files there. If a measurement folder itself contains subfolders, the program reads `.txt` files from its **first** subfolder instead; it does not search all nested folders recursively.

Example — measurements directly inside sample folders:

```text
My_Experiment/
├── Sample_01/
│   ├── test_01.txt
│   ├── test_02.txt
│   └── calibration.txt      # Ignored
├── Sample_02/
│   ├── test_03.txt
│   └── test_04.txt
└── Calib/                    # Ignored as a top-level folder
```

Example — one nested measurement folder:

```text
My_Experiment/
└── Sample_01/
    └── Data/
        ├── test_01.txt
        └── test_02.txt
```

Each measurement file must contain a tab-separated heading line **beginning with `Time`** (and containing tabs), followed by numeric data in this column order:

```text
Time (s)    Load (uN)    Crosshead (nm)
0.00        0.000        0.0
0.10        0.005        30.0
0.20        0.012        60.0
```

*The example above is shown with spaces for readability; actual files must use **tabs** as separators.* Additional text before the `Time` header is allowed. The parser reads the first three tab-separated data fields as time (s), load (µN), and crosshead (nm), respectively. Rows with nonnumeric values in these fields are skipped.

### Calibration exclusion

The current script ignores, **case-insensitively**:

- A measurement file with the exact base name `calibration`, such as `calibration.txt` or `Calibration.TXT`.
- A folder named exactly `calibration` when selecting the top-level measurement folders or a nested data folder.
- Top-level folders named exactly `Figure`, `default`, or `Calib` (with the specified capitalization).

**A filename that merely contains `calibration` is not ignored.** For example, `sample_calibration_02.txt` will still be processed. Rename or remove it if it should be excluded.

## 3. Launch the analyzer

Place the script in a convenient location, then run:

```bash
python piuma_mts_analyzer_r2_cutoff.py
```

The program opens a graphical window. A console/terminal may remain open while it runs.

## 4. Typical workflow

1. **Select directory:** Click **Browse**, and choose the top-level experiment folder containing your measurement subfolders.
2. **Enter probe/sample radius (nm):** Enter the radius in **nanometers**, not micrometers. For example, a spherical probe with radius 110 µm should be entered as `110000` nm. This value applies to every file in the run.
3. **Optional R² threshold:** Enter a value between `0` and `1` in **Auto-exclude fits with R² below** (for example, `0.90`). Leave the field **blank** to disable automatic exclusion.
4. **Click Run:** The program reads the data, performs Hertz fits, generates plots, and populates the results table. Fitting can take time because it uses multiple optimization starting points.
5. **Inspect results:** Use the plot-selector buttons **Crosshead VS Load**, **Time VS Load**, and **Time VS Crosshead**. Use **Previous** and **Next** to navigate measurements, or double-click a table row to jump to that result.
6. **Refine exclusions:** Change the R² threshold and click **Apply R²** (or press Enter) to reapply the threshold **without refitting**. Select one or more table rows and click **Exclude/Include Selected** to manually override the current automatic decision.
7. **Export to Excel:** Click **Export to Excel**. The program uses the currently entered R² threshold and manual overrides when deciding what to include in the workbook.

**Table shading:** Excluded rows are shown in gray. They remain available for review but are omitted from Excel export.

### Automatic R² exclusion rules

| R² result | Cutoff entered: `0.90` | Export behavior |
|---|---|---|
| `0.97` | Pass | Included |
| `0.90` | Pass (equal to cutoff) | Included |
| `0.85` | Fail | Excluded |
| Missing / `NaN` | Fail | Excluded |

- With a **blank** cutoff, no rows are automatically excluded based on R², including failed fits with unavailable R²; you can still exclude these manually.
- A cutoff must be a finite number from **0 to 1**, inclusive.
- **Manual choices take precedence** over the automatic R² rule. For example, a manually included poor fit remains included even if the cutoff is reapplied, until that manual choice is changed or a new analysis is run.
- Clicking **Run** again clears previous manual overrides and recomputes the results.
- R² is used to evaluate fit quality; the program does **not** optimize fits by choosing the highest-R² segment.

## 5. How data are processed and fitted

### 5.1 Unit conversion and force adjustment

The script converts recorded crosshead positions from nm to m, loads from µN to N, and the entered probe radius from nm to m before fitting.

For each measurement, the force series is adjusted as:

```text
F_adjusted[i] = F_raw[i] + abs(F_raw[0])
```

This is the exact behavior of the script. **It is not equivalent to conventional baseline subtraction** (`F_raw[i] - F_raw[0]`) when the first force value is positive.

### 5.2 Select the loading branch

The script identifies the **first prominent load peak** using SciPy's peak-finding routine, with a minimum prominence of 5% of the measured force range and a minimum separation of 1% of the series length. If no qualifying peak exists, it uses the global maximum. The fit includes the data from the **beginning of the recording through that peak**. It does not choose a smaller, highest-R² sliding window.

### 5.3 Hertz model and optimization

The adjusted force is modeled as:

```text
F(x) = C + (4/3) E* sqrt(R) [max(0, s(x - x0))]^(3/2)
```

where:

- `F`: modeled force (N)
- `E*`: fitted effective elastic modulus (Pa)
- `R`: spherical probe radius (m)
- `x`: recorded crosshead position (m)
- `x0`: fitted contact-point offset (m)
- `C`: fitted force offset (N)
- `s`: loading direction, `+1` or `−1`, inferred from the force–crosshead relationship

The algorithm fits `E*`, `x0`, and `C` by bounded nonlinear least squares using multiple initial parameter guesses. It retains the successful solution with the **lowest RMSE** among those attempts. There is no mandatory minimum R² required for a fit to be reported.

The program then calculates its reported modulus using Poisson's ratio **ν = 0.42**:

```text
E_corrected = E* (1 - ν²) = E* (1 - 0.42²)
```

The fitted modulus is **in pascals (Pa)**. This model treats recorded crosshead displacement as the relevant displacement coordinate. Interpretation as an intrinsic Young's modulus requires that the Hertz contact assumptions and the measurement's displacement calibration be appropriate.

### 5.4 Fit quality metrics

- **R²:** Coefficient of determination, calculated over the selected loading data.
- **RMSE (µN):** Root mean square of the force-fitting residuals.
- **NRMSE (%):** `100 × RMSE / (maximum fitted-segment force − minimum fitted-segment force)`.

If there is insufficient usable data or fitting fails, the affected fit outputs are set to `NaN`, while the corresponding measurement can still appear in the results table.

### 5.5 Plasticity-related calculation

The script also computes a displacement-based value named `Plasticity_nm`. It estimates a baseline force over **10–50 s**, adds **0.008 µN** to create a threshold, and searches for upward threshold crossings within **60–65 s**. When crossings are found, it subtracts the first crosshead position from the position at the **last** qualifying crossing. `Plasticity_percent` is that displacement divided by the maximum recorded crosshead position, multiplied by 100, when the denominator is positive.

This is a **script-specific threshold/displacement metric**, rather than an independently validated measure of irreversible plastic deformation. If the required crossing or valid denominator is unavailable, the output may be `NaN`.

### 5.6 Plots

The force–time plot uses a bilateral filter **for display only**; the Hertz fitting uses the adjusted, unfiltered force values. The crosshead–load figure displays the measured curve, highlights the loading interval used for fitting, and overlays the fitted model when available.

## 6. Output files and locations

The program writes results under the **selected experiment folder**:

```text
My_Experiment/
├── output.xlsx
└── Figure/
    ├── crosshead VS load/
    │   └── test_01.pkl
    ├── time VS load/
    │   └── test_01.pkl
    └── time VS crosshead/
        └── test_01.pkl
```

### Excel workbook

- Filename: **`output.xlsx`**
- Worksheet: **`Results`**
- Only rows currently **included** are exported.
- The GUI-only column `Exclude_from_Excel` is **not** written to the workbook.

| Column | Meaning |
|---|---|
| `Test_Name` | Measurement name, derived from the `.txt` filename without its extension |
| `E_corrected_Pa` | Poisson-corrected modulus (Pa) |
| `Hertz_fit_R2` | Coefficient of determination of the Hertz fit |
| `Hertz_fit_RMSE_uN` | Force-fitting RMSE (µN) |
| `Hertz_fit_NRMSE_percent` | Force-fitting normalized RMSE (%) |
| `Fit_Start_Index` | Starting index of the selected loading segment (Python zero-based indexing) |
| `Fit_End_Index` | End boundary of the selected loading segment (Python exclusive indexing) |
| `Fit_Direction` | Inferred direction of loading, `+1` or `−1` |
| `Plasticity_nm` | Script-defined displacement metric (nm) |
| `Plasticity_percent` | Script-defined displacement metric normalized to maximum crosshead position (%) |

### Figure files

The figure outputs are **pickled Matplotlib figures (`.pkl`)**, not PNG or PDF image files. The GUI reloads these files when you select a plot. They are intended for reuse with a compatible Python/Matplotlib environment; avoid loading pickle files from untrusted sources.

**Note:** Running again or exporting to the same experiment folder can overwrite files with the same output names. When different folders contain measurement files with identical base filenames (for example, two different `test_01.txt` files), their plot filenames may collide, and the manual selection behavior can affect all rows with that same `Test_Name`. Use unique measurement filenames across the experiment for predictable behavior.

## 7. Troubleshooting

| Problem | What to check |
|---|---|
| No measurements found | Select the parent folder containing sample subfolders; confirm `.txt` files are within those folders or their first subfolder. |
| Data file not read | Verify the file has a tab-separated line beginning with `Time`, followed by at least three numeric tab-separated fields. |
| A calibration-like file was processed | Only a file whose base name is exactly `calibration` is ignored; names containing the word are not automatically ignored. |
| Modulus or R² is `NaN` | There may be too few valid data points, too little force/displacement variation, or no successful numerical fit. Inspect the source curve. |
| An R²-excluded row still appears | This is expected: the row stays in the GUI but is omitted from Excel. |
| A low-R² row is exported | Check whether it was manually included (manual overrides have priority). |
| Changing the cutoff did not change visible exclusions | Click **Apply R²** or press Enter. **Export to Excel** also applies the current cutoff. |
| Excel export fails | Close `output.xlsx` if it is open in Excel; confirm the selected folder is writable. |
| A plot is missing or points to the wrong measurement | Check for duplicate `.txt` filenames across folders; `.pkl` plots are named by the filename stem only. |
| GUI window cannot open | Confirm Tkinter is installed and that Python is running in a desktop environment. |

## 8. Interpretation and reporting notes

- The results depend on the choice of spherical probe radius and **ν = 0.42**; document both in reports and manuscripts.
- The fit uses the entire start-to-first-peak interval, which may include pre-contact data and other departures from ideal Hertz assumptions. Examine individual fit plots rather than relying solely on R².
- Report the R² threshold used for export and any manual overrides when describing data exclusions.
- The Excel file contains accepted results, not a complete exclusion audit trail. The application table displays exclusions during the current session.
- A reproducible analysis should retain raw `.txt` files, the script version, the selected radius, the R² cutoff, and notes of manual exclusions.

---

**Related document:** `Hertz_Contact_Model_Methods_Editable_Equations.docx` provides manuscript-oriented fitting descriptions and editable Microsoft Word equations.
