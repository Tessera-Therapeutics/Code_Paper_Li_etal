# Public Database Files

This directory contains public cancer gene databases used for annotation.

## Included Files

### MSK Cancer Gene List

**Filename**: `msk.cancerGeneList.tsv`

**Source**: Memorial Sloan Kettering Cancer Center  
**URL**: https://www.cancergenecensus.org/ or https://www.oncokb.org/

**Description**: Curated list of cancer-associated genes classified as:
- Oncogenes
- Tumor suppressors (TSG)
- Both (ONCOGENE_AND_TSG)

**Citation**: Chakravarty et al., JCO Precision Oncology 2017

**Required columns**:
- `Hugo Symbol` - Gene symbol
- `Gene Type` - Classification (ONCOGENE, TSG, ONCOGENE_AND_TSG)
- `# of occurrence within resources (Column J-P)` - Confidence score

**Download instructions** if not included:
```bash
# Visit https://www.oncokb.org/cancerGenes and download the TSV file
# Or use the Cancer Gene Census from COSMIC (requires registration)
```

---

### CancerMine Database

**Filename**: `cancermine_collated.tsv`

**Source**: CancerMine text-mining database  
**URL**: http://bionlp.bcgsc.ca/cancermine/

**Description**: Automatically extracted gene-cancer associations from literature.

**Citation**: Lever et al., Nature Methods 2019

**Required columns**:
- `gene_normalized` - Gene symbol
- `role` - Gene role (Oncogene, Tumor_Suppressor, Driver)
- `cancer_normalized` - Associated cancer type
- `citation_count` - Supporting citations

**Download instructions**:
```bash
cd metadata/public_databases/
wget http://bionlp.bcgsc.ca/cancermine/cancermine_collated.tsv
```

---

### UniProt Oncogene/Tumor Suppressor Lists

**Filenames**: 
- `uniprot-Proto-oncogene.tsv`
- `uniprot-Tumor-suppressor.tsv`

**Source**: UniProt database  
**URL**: https://www.uniprot.org/

**Description**: Manually curated protein annotations from UniProt.

**Download instructions**:
```bash
# Visit https://www.uniprot.org/
# Search: keyword:"Proto-oncogene" OR keyword:"Tumor suppressor"
# Filter: Reviewed (Swiss-Prot)
# Download: TSV format with columns:
#   - Gene names (primary)
#   - Gene names
#   - Organism
#   - Status

# Or use the API:
cd metadata/public_databases/
wget "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=keyword:Proto-oncogene%20AND%20reviewed:true" \
  -O uniprot-Proto-oncogene.tsv

wget "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=keyword:Tumor%20suppressor%20AND%20reviewed:true" \
  -O uniprot-Tumor-suppressor.tsv
```

---

## Usage in Analysis

These files are used by the `isa_utils.R` functions:

```r
# Load oncogenes and tumor suppressors
oncots_db <- loadOncoTSgenes(
  onco_db_file = "metadata/public_databases/uniprot-Proto-oncogene.tsv",
  tumsup_db_file = "metadata/public_databases/uniprot-Tumor-suppressor.tsv",
  species = "human"
)

# Or use CancerMine
cancermine_db <- loadOncoTSgenes_fromCancerMine(
  cancermine_collated_file = "metadata/public_databases/cancermine_collated.tsv"
)
```

## File Formats

All files should be tab-delimited (TSV) with headers. Ensure your files match the expected column names in the analysis scripts.

## License & Citation

These are **public databases**. When using them in publications, please cite the original sources listed above.

## Updates

Database versions used in the paper:
- MSK Cancer Gene List: 2026-01 version
- CancerMine: 2019-10 version
- UniProt: 2018-06 version

For reproducibility, we include the specific versions used. For new analyses, you may want to download updated versions.
