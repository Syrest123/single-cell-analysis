# **Single-Cell RNA-Seq Preprocessing Pipeline**  
*Demultiplexing with `bcl2fastq` and counts generation with `cellranger multi`*

## **Overview**
1. **Demultiplexing**: Convert BCL to FASTQ (`bcl2fastq`)  
2. **Counts generation**: Multi-sample processing (`cellranger multi`)  
3. **Quality control**: 10X Genomics QC metrics  

## **Pipeline Steps**
### 1. Demultiplexing (BCL → FASTQ)
```bash
bcl2fastq --runfolder-dir BCL_DIR \
          --output-dir FASTQ_OUTPUT \
          --sample-sheet SampleSheet.csv \
          --no-lane-splitting
```
*Output*:  
- `FASTQ_OUTPUT/` (per-sample FASTQs)

### 2. Counts Generation (`cellranger multi`)
- Configure `multi_config.csv`:
  ```csv
  [gene-expression]
  reference=/path/to/refdata-gex-GRCh38-2020-A

  [libraries]
  fastq_id,fastqs,feature_types
  sample1,/path/to/fastq,Gene Expression

  [samples]
  sample_id,description
  sample1,control
  ```
- Run:
  ```bash
  cellranger multi --id=OUTPUT_DIR \
                   --csv=multi_config.csv
  ```

### 3. Quality Control
Check in `OUTPUT_DIR/outs/per_sample_outs/`:  
- `web_summary.html` (key metrics)  
- `metrics_summary.csv` (tabular QC)  

## **Outputs**
- `filtered_feature_bc_matrix.h5` (counts per sample)  
- `cloupe.cloupe` (Loupe Browser visualization)  

## **References**
- [Illumina bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)  
- [10X cellranger multi](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi)  

---
