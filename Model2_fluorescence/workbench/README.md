# Fluorescence Unmixing Workbench

A browser interface for applying the M3 two-channel fluorescence calibration in the Model section of ShanghaiTech-China's iGEM 2026 project. It estimates the GFP-positive fraction from GFP and mCherry plate-reader readings.

## Files

- `index.html`: workbench interface and prediction code.
- `README.md`: usage and model scope.

## Use

Download the HTML file and open it in a browser. GitHub's file view displays source code, rather than running the workbench.

1. Set the blank and test wells on the 96-well layout, or click **Load example layout**.
2. Download the layout template and enter the measurements, or prepare an Excel/CSV file with the columns below.
3. Upload the file, check the matched wells, and click **Start analysis**.
4. Inspect the plate map and per-well results. Click **Download results (.csv)** to export them.

| Column | Meaning |
| --- | --- |
| `Well` | Well identifier, such as A1 or H12 |
| `GFP_raw` | Raw GFP-channel fluorescence |
| `mCherry_raw` | Raw mCherry-channel fluorescence |

Include measurements for the blank wells selected in the layout. The workbench subtracts their mean reading separately for each channel, then applies:

```text
fraction = 2.6363458858e-7 * GFPc
         - 5.5041279846e-6 * mCherryc
         + 0.1932729405
```

## Requirements

No Python installation or local server is required. This HTML version loads SheetJS Community Edition 0.18.5 from cdnjs for spreadsheet import/export, and fonts from Google Fonts. 

Measurement processing runs in the browser; the workbench code does not upload the measurement file to a server. External resource loading still makes requests to the dependency providers.

## Model scope

Use samples with comparable strains, expression states, cell concentrations and instrument settings. Specified volume fractions approximate cell fractions only when cell concentrations are comparable.

Out-of-range signals are flagged. Plate colors and the mean summary use clipped values, while the per-well table retains raw predictions and the CSV includes both raw and clipped predictions. 
