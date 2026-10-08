# Troubleshooting

**Repository not found:** run from the checkout or pass `acgut --root /absolute/path/to/checkout COMMAND ...`. The package requires the accompanying data/configuration.

**Output directory already populated:** choose a new run name. This protects previous runs and the released data.

**"Inside the checkout, write outputs under …/outputs":** every command refuses output inside the checkout except below `outputs/` (which git ignores). Use `outputs/<name>` or any directory outside the checkout.

**Missing package/RDKit:** activate the environment and install the appropriate requirements file. RDKit is optional for data-based result reproduction, required for building molecules. Gaussian/CREST are external programs, not Python dependencies.

**Record checksum mismatch:** check for unintended file changes. Recover the released bytes rather than editing manifests to suppress a failure.

**Missing species, duplicate fixed product, wrong method/state, a family panel that no longer matches `family_count`, or incomplete primary panel:** inspect configuration and calculation manifests. The workflow stops rather than silently shrinking the comparison or selecting a diagnostic structure.

**A comparison fails:** `reproduce` verifies this study's published results against the bundled records. If a check fails, something changed — inspect the run output and the records before anything else.

**Plot/font differences:** numerical table/point comparisons are authoritative. Font rendering can vary by platform. The Figure 1 canvas is 17.8 × 13.2 cm at 1200 dpi; PDF/SVG are vector exports. Use `MPLBACKEND=Agg` in unattended environments. `--dpi 200` is a quick preview; it changes only PNG resolution.
