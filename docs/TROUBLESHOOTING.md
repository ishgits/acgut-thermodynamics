# Troubleshooting

**Repository not found:** run from the checkout or pass `atcgu --root /absolute/path/to/checkout COMMAND ...`. The package requires the accompanying data/configuration.

**Output directory already populated:** choose a new run name. This protects previous runs and the released data.

**"Inside the checkout, write outputs under …/outputs":** every command refuses output inside the checkout except below `outputs/` (which git ignores). Use `outputs/<name>` or any directory outside the checkout.

**Missing package/RDKit:** activate the environment and install the appropriate requirements file. RDKit is optional for data-based result reproduction, required for building molecules. Gaussian/CREST are external programs, not Python dependencies.

**Integrity or record checksum mismatch:** check for unintended file changes. Do not update the checksum just to suppress a failure. Recover the released bytes, or create and version a new dataset intentionally.

**Missing species, duplicate fixed product, wrong method/state, a family panel that no longer matches `family_count`, or incomplete primary panel:** inspect configuration and calculation manifests. The workflow stops rather than silently shrinking the comparison or selecting a diagnostic structure.

**New Gaussian extraction fails:** inspect the extraction audit and `record.json` failure reasons. Unsupported elements, missing/imaginary modes, unfinished jobs or ambiguous concatenated jobs need correction/manual review. Failed candidates can be retained but cannot supply selected thermochemistry.

**CREST parse/refinement failure:** a multi-XYZ ensemble must have complete atom blocks, a finite numeric energy in every comment, matching atom order/composition and the energy minimum first. Source conformer numbers are one-based.

**Plot/font differences:** numerical table/point comparisons are authoritative. Font rendering can vary by platform. The Figure 1 canvas is 17.8 × 13.2 cm at 1200 dpi; PDF/SVG are vector exports. Use `MPLBACKEND=Agg` in unattended environments. `--dpi 200` is a quick preview; it changes only PNG resolution.

**New dataset does not match the published results:** `reproduce` verifies this study's published results. For an intentional extension use a separate registry/analysis ID and `analyze`, then establish its own reference results and plots.
