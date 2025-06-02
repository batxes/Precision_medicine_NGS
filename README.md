# 🧬 Targeted NGS Analysis: RNA-Seq & DNA Amplicon Data

This repository contains two separate next-generation sequencing (NGS) analyses, focusing on targeted datasets:

- **Task 1**: Analysis of a **targeted RNA-seq dataset**, including gene enrichment, variant calling, and interpretation of soft-clipped reads.
- **Task 2**: Analysis of **amplicon-based DNA sequencing data**, covering read mapping, mutation detection, and coverage evaluation across specific gene regions.

Each task is structured with code, outputs, and a final report summarizing insights and results.

---

## 📂 Tasks Overview

### 🧪 Task 1 — Targeted RNA-Seq Analysis

- Quality control and exploration of aligned reads (sorted BAM format)
- Identification of enriched genes and barcoded read quantification
- Variant analysis (SNVs and structural variations)
- Investigation of soft-clipped sequences (>15 bp) and their possible biological implications

### 🔬 Task 2 — DNA Amplicon Sequencing Analysis

- Mapping analysis of FASTQ reads and interpretation of BAM alignment
- Variant calling and prioritization using technical and biological filters
- Evaluation of read constitution (adapters, primers, amplified region)
- Amplicon-level coverage profiling and hypotheses for coverage variation

---

## 🛠️ Technologies & Tools Used

| Category                  | Tools / Libraries                                                 |
|---------------------------|-------------------------------------------------------------------|
| **NGS Tools**             | `samtools`, `bcftools`, `FreeBayes`, `GATK`, `FastQC`             |
| **Python Libraries**      | `pysam`, `pandas`, `numpy`, `matplotlib`, `seaborn`               |
| **Coverage & Composition**| `BEDTools`, `pysamstats`                                          |
| **Visualization**         | `IGV`, UCSC Genome Browser, Python plots                          |
| **Variant Annotation**    | `VEP`, `SnpEff`, or `ANNOVAR` *(optional)*                        |
| **Reporting**             | `Jupyter Notebooks`, Markdown, LaTeX                              |

---

## 📁 Repository Structure
```
targeted-ngs-analysis/
├── rna_seq/
│ ├── bam_files/
│ ├── scripts/
│ ├── figures/
│ └── Task1_Report.pdf
├── amplicon/
│ ├── fastq/
│ ├── bam/
│ ├── scripts/
│ ├── figures/
│ └── Task2_Report.pdf
├── README.md
└── requirements.txt
```

---

## 📄 How to Use

1. Clone the repository and install dependencies:
    ```
    git clone https://github.com/batxes/Precision_medicine_NGS.git
    cd Precision_medicine_NGS/RNA-seq #or Amplicon
    pip install -r requirements.txt
    ```
2. Run the analysis notebooks or scripts found in each task folder.
3. Open the reports in PDF format for a detailed description of the analysis, figures, and interpretations.


