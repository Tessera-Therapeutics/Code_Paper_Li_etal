# GM25256 CRAM Locations (S3)

All CRAMs are name-sorted and require coordinate sorting before use with DRAGEN.
Experiment: EXP21001214 | Reference: hg38

## B samples (experiment)

| Sample | S3 Path |
|--------|---------|
| GM25256-B1 | `s3://fvl58-compbio/archive/baltshuler/baltshuler-wgs-01/scratch/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B1/GM25256-B1.name-sorted.cram` |
| GM25256-B2 | `s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-02/from-s3/2021-08_rt-wgs_EXP21001214/GM25256-B2/GM25256-B2.name-sorted.cram` |
| GM25256-B3 | `s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-03/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B3/GM25256-B3.name-sorted.cram` |
| GM25256-B4 | `s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-04/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B4/GM25256-B4.name-sorted.cram` |
| GM25256-B5 | `s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-05/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B5/GM25256-B5.name-sorted.cram` |
| GM25256-B6 | `s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-06/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-B6/GM25256-B6.name-sorted.cram` |

## C samples (control)

| Sample | S3 Path |
|--------|---------|
| GM25256-C8 | `s3://fvl58-compbio/analyses/baltshuler/wgs-ec2-storage-to-sort/wgs-06/from-s3/rt-platform/rt-wgs/2021-08_rt-wgs_EXP21001214/GM25256-C8/GM25256-C8.name-sorted.cram` |

## Notes

- B1 is in the **archive** bucket (`fvl58-compbio/archive/`) — may have slower download speeds (Glacier retrieval may be needed).
- B2 and C8 BAMs are already processed and ready on the DRAGEN instance.
- Pipeline per sample: `aws s3 cp` → `samtools sort` → `samtools index` → `samtools view -b` (CRAM→BAM) → `samtools index`
- Process one sample at a time due to disk space constraints; move or delete intermediate files before starting the next.
