# Attribution and software

Structure acquisition and Gaussian input generation were adapted from Ishaan
Madan's **PubChem Gaussian Pipeline**, v2.1.0, DOI
[10.5281/zenodo.21365151](https://zenodo.org/records/21365151).
Upstream code and citation are available at
[ishgits/pubchem-gaussian-pipeline](https://github.com/ishgits/pubchem-gaussian-pipeline).
The study's explicit CID choices, neutral linkage builders, source atom indices,
sampling inputs and calculation routes define the adaptations here. The study
is not represented as having executed the upstream release unchanged.

Frozen starting structures acquired from PubChem retain PubChem metadata and
compound identifiers in their SDF files. `data/source/manifest.csv` identifies
the source type and file hash for each structure. A manuscript link and its
source/CID table will be added when public. Embedded CIDs do not establish
a new retrieval date or guarantee that current PubChem responses match frozen bytes.

The runtime uses NumPy, pandas, Matplotlib and PyYAML. Molecule construction uses
optional RDKit. External calculations use Gaussian 16 and CREST/GFN2-xTB; none
of their executables is distributed. Each dependency retains its own license.
See the per-record Gaussian revision and pinned Python environment for tested
versions.
