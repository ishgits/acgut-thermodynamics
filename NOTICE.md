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
compound identifiers. Source/CID mappings and file hashes are provided under
`data/source/`. Some source URLs are reconstructed API locators from embedded
PubChem CIDs where the original acquisition CSV lacked a row; they do not establish
a new retrieval date or guarantee that the current response matches frozen bytes.

The runtime uses NumPy, pandas, Matplotlib and PyYAML. Molecule construction uses
optional RDKit. External calculations use Gaussian 16 and CREST/GFN2-xTB; none
of their executables is distributed. Each dependency retains its own license.
See the per-record Gaussian revision and pinned Python environment for tested
versions.
