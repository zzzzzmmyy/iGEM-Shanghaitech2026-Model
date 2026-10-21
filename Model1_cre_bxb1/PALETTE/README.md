# PALETTE

A browser-based companion to the MOSAIC Cre–Bxb1 model for comparing circuit designs, visualizing predicted within-well cell fractions, and generating draft GenBank files.

## Open the tool

Download this folder and open `index.html` in a modern desktop browser such as Chrome, Edge, or Firefox. No Python installation, package installation, or backend server is required. GitHub's file viewer shows the source; download the HTML to use the interface.

## Inputs

| Input | Meaning | Default |
| --- | --- | --- |
| Fraction min / max | Desired A-type (amilGFP) fraction range among differentiated A/B cells; enter percentages | 20–60% |
| Colony picks | Number of independently picked colonies, one colony per well | 96 |
| Left arm (bp) | Left-arm length supplied to the prediction model | 773 bp |
| Fragment name | Label used for the payload and exported files | YourFavouriteGene |
| Readout time | Time after Bxb1 induction | 24 h |
| Payload sequence | Optional DNA or FASTA input for sequence construction | Empty |
| Insertion arm | Payload insertion on the amilGFP or mCherry side | Right / mCherry |

Changes update the results automatically. Click a row in Ranked Designs to select a candidate.

The sequence parser removes FASTA headers and retains only A, C, G, T, and N. Other ambiguity symbols are discarded; check the displayed sequence length before export.

## Results

- **Ranked Designs:** the top 12 candidates from a scan of 1–10 lox71 sites and 150–2,000 bp NS spacers in 50 bp increments. Ranking prioritizes the number of usable fraction bins inside the target range, then the overlap of the predicted fraction span with that range, then shorter array length, then fewer sites.
- **Predicted Ratio Ladder:** the A-type fraction associated with each displayed Cre architecture. The `w` labels are expected well counts, not observed integer counts.
- **Plate Simulation:** a fixed-seed illustrative plate layout. Only the first 96 picks are displayed. Colors illustrate composition and are not measured fluorescence colors.
- **Ratio types:** displayed fractions grouped into fixed 5-percentage-point bins, with expected counts and the probability of observing at least one well of that type.
- **P(hit):** probability of obtaining at least one displayed usable architecture within the requested fraction range, rather than the probability of covering every fraction type.
- **Plasmid Map:** a schematic of the generated sequence with selected annotations.


## Export

| Button | Download |
| --- | --- |
| Plasmid .gb | Draft circular plasmid sequence with feature annotations |
| Synthesis fragment .gb | Draft linear fragment from the start of amilGFP through the final lox71 site |
| Show text | An abbreviated preview of the plasmid GenBank text; it is not the complete export |

Use the download buttons for complete files. A compatible sequence viewer can open the `.gb` files.

## Interpretation and limitations

Predictions describe colony-derived wells after Bxb1 induction under the model assumptions.

Sequences are draft designs based on the embedded backbone and spacer template. Verify sequence coordinates, recombination-site orientations, reading frames, junctions, actual arm lengths, and assembly feasibility before synthesis or cloning. The tool does not verify payload expression or function. Model parameters and predicted fractions require experimental calibration and validation.

## Files

- `index.html`: interface, JavaScript model, embedded sequence template, SVG visualizations, and GenBank export.
- `README.md`: usage and interpretation.


