# Sample Metadata Files

This directory contains experimental metadata and sample annotations.

## Files Needed

### 1. Master Sample File

**Filename**: `MasterFile_DISTALseq.xlsx` or similar

**Purpose**: Maps sample IDs to experimental conditions, replicates, and file locations.

**Required columns**:
- `SampleID` - Unique sample identifier
- `SampleName` - Descriptive sample name
- `FileName` - Name of the input data file
- `FileFolder` - Path to data files (relative to project root)
- `ToProcess` - Boolean (TRUE/FALSE) indicating if sample should be analyzed
- `Replica` - Replicate number
- `Vector` - Vector construct used
- `PoolID` - Pool/batch identifier
- `GroupName` - Experimental group
- `Transgene` - Transgene name
- `SampleType` - Sample type (e.g., "gDNA", "edited", "control")
- `Separator` - File delimiter ("," or "\\t")
- `Promoter_start`, `Promoter_end` - Vector promoter coordinates
- `pA_start`, `pA_end` - PolyA signal coordinates
- `Functional_min`, `Functional_max` - Functional region boundaries
- `VectorOrientation` - Vector orientation
- `MinReadLen` - Minimum read length threshold
- `MinTargetAlmLen` - Minimum target genome alignment length
- `CommonPrefixAllFiles` - Common prefix for output files
- `DNA_ID` - DNA sample identifier

**Example structure**:

| SampleID | SampleName | Vector | Transgene | SampleType | ToProcess | MinReadLen |
|----------|------------|--------|-----------|------------|-----------|------------|
| RTE25_B1 | RTE-25_Rep1 | RTE-25 | GeneX | gDNA | TRUE | 1000 |
| RTE3_B1 | RTE-3_Rep1 | RTE-3 | GeneY | gDNA | TRUE | 1000 |

### 2. Clone/IS Annotation Files

**Filename**: `Clones.knownIS.xlsx` or similar

**Purpose**: Annotates known integration sites in clonal populations.

**Required columns**:
- Clone identifiers
- Chromosome
- Integration locus
- Strand
- Gene name
- ToProcess flag

### 3. Vector Annotations

If using separate vector file:

**Required columns**:
- `VectorName_BAM` - Vector name as it appears in BAM files
- `UTR5_start`, `UTR5_end` - 5' UTR coordinates
- `UTR3_start`, `UTR3_end` - 3' UTR coordinates
- `Promoter_start`, `Promoter_end`
- `pA_start`, `pA_end`

### 4. Flanking Sequence Data

If applicable:

**Columns**:
- `CommonPrefixAllFiles` - Links to IS sequence data
- ReadID
- Sequence information

## Creating Your Own Metadata Files

1. **Start with templates** in this directory
2. **Fill in your experimental data**
3. **Update the paths** in `config.R` to point to your files
4. **Ensure file formats match** the expected column names

## Example Templates

See `template_master_file.xlsx` and `template_clone_annotations.xlsx` for examples.

## Path Configuration

Once you've added your metadata files here, update the paths in the main analysis scripts or create a `config.R` file:

```r
# In config.R or at the top of your analysis scripts
MASTER_FILE <- "metadata/sample_info/MasterFile_DISTALseq.xlsx"
CLONE_ANNOTATIONS <- "metadata/sample_info/Clones.knownIS.xlsx"
```
