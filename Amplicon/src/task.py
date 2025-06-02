import pandas as pd
import pybedtools as pb
import matplotlib.pyplot as plt
import numpy as np
import sh
import os
import sys
import re
import time
from matplotlib.colors import LogNorm
from collections import Counter

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

## set commands ##
samtools_cmd = sh.Command("samtools")
less_cmd = sh.Command("less")
head_cmd = sh.Command("head")
annotatePeaks_cmd = sh.Command("annotatePeaks.pl")
gatk_cmd = sh.Command("~/Downloads/gatk-4.6.1.0/gatk")

def bam_QC(bam_file_path) -> str:
    """Function that applies extra filters to the bam file.
    
    Input the path to the bam file. Output the name of the filtered bam file
    """
    start = time.time()
    bam_name = bam_file_path.split("/")[-1]
    name = ".".join(bam_name.split(".")[:-1])
    samtools_cmd("view","-b","-q","30","-f","2","-F","4","../sequenceFile/{}".format(bam_name), _out= "../sequenceFile/{}.filtered.bam".format(name))
    samtools_cmd("stats","../sequenceFile/{}.filtered.bam".format(name),_out="../output/{}_filtered_stats.txt".format(name))
    samtools_cmd("index","../sequenceFile/{}.filtered.bam".format(name))
    end = time.time()
    print (f"bam QC -> Succesfull [elapsed {end - start:.4f} seconds]")
    return "{}.filtered.bam".format(name)

def remove_under_covered_regions(bam_file_name,min_reads) -> str:
    """Function that removes covered regions with a specific depth.
    
    Input bam file name and the minimum number of reads in the covered region. Output the name of the filtered bam file
    """
    start = time.time()
    name = ".".join(bam_file_name.split(".")[:-1])
    bamfile_bed = pb.BedTool("../sequenceFile/{}".format(bam_file_name)) 
    covered = bamfile_bed.genome_coverage(bg=True)
    df = covered.merge(c=4,o="median").to_dataframe()
    high_coverage_df = df[df['name'] >= min_reads]
    high_coverage_df[['chrom', 'start', 'end']].to_csv("../output/{}.high_coverage.bed".format(name), sep="\t", header=False, index=False)
    samtools_cmd("view", "-b", "-L", "../output/{}.high_coverage.bed".format(name), "../sequenceFile/{}".format(bam_file_name),_out= "../sequenceFile/{}.high_coverage.bam".format(name))
    end = time.time()
    print (f"remove under covered regions -> Succesfull [elapsed {end - start:.4f} seconds]")
    return "{}.high_coverage.bam".format(name)

def get_genes_from_bam(bam_file_name) -> str:
    """Function that using the bam alignment file, gets the coverage of the genome and outputs the genes where the reads overlap.
    
    It also generates a plot showing the position of the genes.
    In the process, it creates a bed file with the coverage (using bedtools) and a text file with the list of genes
    Input bam file name . Output the coverage bed file
    """
    start = time.time()
    name = ".".join(bam_file_name.split(".")[:-1])
    bamfile_bed = pb.BedTool("../sequenceFile/{}".format(bam_file_name)) 
    covered = bamfile_bed.genome_coverage(bg=True)
    df = covered.to_dataframe()
    #df["chrom"] = "chr"+df["chrom"]
    covered = pb.BedTool.from_dataframe(df)

    gtf = pb.BedTool("../annotations/hg19.refGene.tsv") 
    overlapped_genes = covered.intersect(gtf,wb=True)
    overlapped_genes_df = overlapped_genes.to_dataframe(names=["chrom1", "start1", "end1","number_reads","chrom2", "start2", "end2", "name"])
    covered_genes = list(set(overlapped_genes_df["name"]))
    with open('../output/genes_overlapping_with_{}.txt'.format(name), 'w') as f:
        for line in covered_genes:
            f.write(f"{line}\n")
    
    # do the figure
    fig,axs = plt.subplots(24,1,figsize=(20,10),sharex=True,sharey=True)
    all_chr = ["chr"+str(i) for i in range(1,23)]
    all_chr.append("chrX")
    all_chr.append("chrY")
    n_reads = list(overlapped_genes_df["number_reads"])
    vmax = max(overlapped_genes_df["number_reads"])
    vmin = min(overlapped_genes_df["number_reads"])
    for i, chr_ in enumerate(all_chr):
        subset_df = overlapped_genes_df.query("chrom1 == '{}'".format(chr_))
        scatter = axs[i].scatter(subset_df["start1"],np.zeros(len(subset_df)),cmap="viridis_r",c=subset_df["number_reads"],s=20,norm=LogNorm(vmin=vmin, vmax=vmax))
        axs[i].set_ylabel(chr_)
    axs[i].set_xlabel("Genomic position")
    cax = fig.add_axes([0.65, 0.3, 0.25, 0.03]) 
    cbar = plt.colorbar(scatter,cax=cax,orientation='horizontal')
    cbar.set_label('Number of reads (Log10)')
    plt.savefig("../output/genes_position_{}.pdf".format(name))

    fig,axs = plt.subplots(1,24,figsize=(20,10),sharex=False,sharey=True)
    for i, chr_ in enumerate(all_chr):
        subset_df = overlapped_genes_df.query("chrom1 == '{}'".format(chr_))
        scatter = axs[i].scatter(subset_df["start1"],subset_df["number_reads"],cmap="viridis_r",c=subset_df["number_reads"],s=20,norm=LogNorm(vmin=vmin, vmax=vmax))
        axs[i].set_xlabel(chr_)
        axs[i].spines[['right', 'top']].set_visible(False)
    axs[0].set_ylabel("Number of reads")
    cax = fig.add_axes([0.65, 0.3, 0.25, 0.03]) 
    cbar = plt.colorbar(scatter,cax=cax,orientation='horizontal')
    cbar.set_label('Number of reads (Log10)')
    plt.savefig("../output/genes_position_manhattan_{}.pdf".format(name))

    # prepare this step for second function, where we annotatePeaks

    #bedtools merge -i sorted.bed -c 4 -o median > merged.bed

    df = covered.merge(c=4,o="median").to_dataframe()
    df.to_csv("../output/{}.genome_coverage_merged.tsv".format(name),sep="\t")
    df["strand"] = ["+"]*(len(df))
    df=df[["chrom","start","end","strand"]] #output for homer
    df.to_csv("../output/bam_coverage_{}.bed".format(name),header=False,sep="\t")
    end = time.time()
    print (f"get genes from coverage -> Succesfull [elapsed {end - start:.4f} seconds]")
    return "bam_coverage_{}.bed".format(name)
    
def get_gene_regions_from_bam_coverage(bam_coverage_file_name):
    """Function that characterizes the location of the aligned reads in respect to the genes.

        Input, the bed file with the coverage of the bam. It has no output. 
    """    
    start = time.time()
    name = ".".join(bam_coverage_file_name.split(".")[:-1])
    #"../output/bam_coverage_{}.bed".format(bam_file_name)
    #annotatePeaks.pl bam_coverage.bed hg19 > bam_coverage_annotated.txt
    print ("../output/{}".format(bam_coverage_file_name))
    annotatePeaks_cmd("../output/{}".format(bam_coverage_file_name),"hg19",_out="../output/bam_coverage_annotated_{}.tsv".format(name))

    annotated_df = pd.read_csv("../output/bam_coverage_annotated_{}.tsv".format(name),sep="\t")
    #annot_labels = ["3'", "5'", 'Intergenic', 'TTS', 'exon', 'intron', 'non-coding', 'promoter-TSS']

    annotations = annotated_df.Annotation
    annotations_list = []
    for annot  in annotations:
        try:
            annotations_list.append(annot.split(" ")[0])
        except:
            print (f"{annot} can not be split")

    annotations_set = Counter(annotations_list)

    labels = list(annotations_set.keys())
    values = list(annotations_set.values())

    # bar plot of gene regions counts
    fig, ax = plt.subplots()
    ax.bar(labels, values, color='darkblue')
    ax.set_xlabel('Labels')
    ax.set_ylabel('Count')
    ax.set_title('Annotation Counts')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("../output/gene_regions_{}.pdf".format(name))

    # get gene names of reads falling in promoters, TSS and exons
    annotated_df = annotated_df[annotated_df['Annotation'].notna()]
    subset_labels = ["promoter","TSS","exon"]
    subset_df = annotated_df[annotated_df['Annotation'].str.startswith(tuple(subset_labels))]
    subset_df.to_csv("../output/gene_regions_overlapping_with_{}.tsv".format(name),sep="\t")
    print ("Genes found in exons/Promoters/TSS")
    print (sorted(list(set(subset_df["Gene Name"]))))
    end = time.time()
    print (f"get gene regions from coverage -> Succesfull [elapsed {end - start:.4f} seconds]")

def get_variants_with_mutect2(bam_file_name):
    """ Function that given a bam file, it computes the somatic variants using mutect2 and it also add extra filtering process.

        Input is the name of the bam file. There is no output.
    """
    start = time.time()
    name = ".".join(bam_file_name.split(".")[:-1])

    # Not needed with GATK tools
    #samtools_cmd("view","-H","../sequenceFile/{}".format(bam_file_name),_out="../sequenceFile/header.txt")
    #samtools_cmd("reheader","../sequenceFile/modified_header.txt","../sequenceFile/{}".format(bam_file_name),_out="../sequenceFile/{}.newheader.bam".format(name))
    #samtools_cmd("index","../sequenceFile/{}.newheader.bam".format(name))
    try:
        gatk_cmd("CreateSequenceDictionary","-R","../annotations/hg19.fa","-O","../annotations/hg19.dict")
    except:
        #print ("hg19.dict already exists")
        pass
    #gatk_cmd("Mutect2","-R","../annotations/hg19.fa","-I","../sequenceFile/{}.newheader.bam".format(name),"-O","../output/{}.mutect2_calling.vcf.gz".format(name))
    # add read group tp @RG tag and reindex the bam
    samtools_cmd("addreplacerg", "-r", "@RG\tID:samplename\tSM:samplename", "../sequenceFile/{}.bam".format(name), "-o" "../sequenceFile/{}.renamed.bam".format(name))
    samtools_cmd("index", "../sequenceFile/{}.renamed.bam".format(name))
    gatk_cmd("Mutect2","-R","../annotations/hg19.fa","-I","../sequenceFile/{}.renamed.bam".format(name),"-O","../output/mutect2_calling_{}.renamed.vcf.gz".format(name))
    gatk_cmd("FilterMutectCalls","-V","../output/mutect2_calling_{}.vcf.gz".format(name),"-R","../annotations/hg19.fa","-O", "../output/mutect2_calling_{}.passed.vcf.gz".format(name))
    # This step I skipped, already done, takes too long
    # ~/Downloads/gatk-4.6.1.0/gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
    end = time.time()
    print (f"get variants with mutect2 -> Succesfull [elapsed {end - start:.4f} seconds]")

def get_variants_relevance_funconator(bam_file_name):
    """ Function uses the vcf file generated in the mutect2 calling and annotates teh variants against clinical databases.

        Input is the name of the bam file. There is no output.
    """
    start = time.time()

        #"--java-options", "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true",
   
    name = ".".join(bam_file_name.split(".")[:-1])

    gatk_cmd(
        "Funcotator",
        "--variant", "../output/mutect2_calling_{}.passed.vcf.gz".format(name),
        "--reference", "../annotations/hg19.fa",
        "--ref-version", "hg19",
        "--data-sources-path", "../annotations/funcotator_dataSources.v1.8.hg19.20230908s",
        "--output", "../output/variants.funconator_{}.passed.vcf".format(name),
        "--output-file-format", "VCF"
    )
    end = time.time()
    print (f"get variants relevance with funconator -> Succesfull [elapsed {end - start:.4f} seconds]")


def filter_variants_for_diseases(bam_file_name):
    """ Function that postprocesses the vcf file with the clinical annotations.

        Input is the name of the bam file. There is no output.
    """
    start = time.time()

    name = ".".join(bam_file_name.split(".")[:-1])

    print (name)
    funcotator_fields = "Gencode_43_hugoSymbol|Gencode_43_ncbiBuild|Gencode_43_chromosome|Gencode_43_start|Gencode_43_end|Gencode_43_variantClassification|Gencode_43_secondaryVariantClassification|Gencode_43_variantType|Gencode_43_refAllele|Gencode_43_tumorSeqAllele1|Gencode_43_tumorSeqAllele2|Gencode_43_genomeChange|Gencode_43_annotationTranscript|Gencode_43_transcriptStrand|Gencode_43_transcriptExon|Gencode_43_transcriptPos|Gencode_43_cDnaChange|Gencode_43_codonChange|Gencode_43_proteinChange|Gencode_43_gcContent|Gencode_43_referenceContext|Gencode_43_otherTranscripts|Achilles_Top_Genes|CGC_Name|CGC_GeneID|CGC_Chr|CGC_Chr_Band|CGC_Cancer_Somatic_Mut|CGC_Cancer_Germline_Mut|CGC_Tumour_Types__(Somatic_Mutations)|CGC_Tumour_Types_(Germline_Mutations)|CGC_Cancer_Syndrome|CGC_Tissue_Type|CGC_Cancer_Molecular_Genetics|CGC_Mutation_Type|CGC_Translocation_Partner|CGC_Other_Germline_Mut|CGC_Other_Syndrome/Disease|ClinVar_HGMD_ID|ClinVar_SYM|ClinVar_TYPE|ClinVar_ASSEMBLY|ClinVar_rs|ClinVar_VCF_AF_ESP|ClinVar_VCF_AF_EXAC|ClinVar_VCF_AF_TGP|ClinVar_VCF_ALLELEID|ClinVar_VCF_CLNDISDB|ClinVar_VCF_CLNDISDBINCL|ClinVar_VCF_CLNDN|ClinVar_VCF_CLNDNINCL|ClinVar_VCF_CLNHGVS|ClinVar_VCF_CLNREVSTAT|ClinVar_VCF_CLNSIG|ClinVar_VCF_CLNSIGCONF|ClinVar_VCF_CLNSIGINCL|ClinVar_VCF_CLNVC|ClinVar_VCF_CLNVCSO|ClinVar_VCF_CLNVI|ClinVar_VCF_DBVARID|ClinVar_VCF_GENEINFO|ClinVar_VCF_MC|ClinVar_VCF_ORIGIN|ClinVar_VCF_RS|ClinVar_VCF_ID|ClinVar_VCF_FILTER|Cosmic_overlapping_mutations|CosmicFusion_fusion_genes|CosmicFusion_fusion_id|CosmicTissue_total_alterations_in_gene|CosmicTissue_tissue_types_affected|DNARepairGenes_Activity_linked_to_OMIM|DNARepairGenes_Chromosome_location_linked_to_Genome_Data_Viewer|DNARepairGenes_Accession_number_linked_to_NCBI_Entrez|Familial_Cancer_Genes_Syndrome|Familial_Cancer_Genes_Synonym|Familial_Cancer_Genes_Reference|Gencode_XHGNC_hgnc_id|Gencode_XRefSeq_mRNA_id|Gencode_XRefSeq_prot_acc|HGNC_HGNC_ID|HGNC_Approved_name|HGNC_Status|HGNC_Locus_type|HGNC_Locus_group|HGNC_Previous_symbols|HGNC_Previous_name|HGNC_Alias_symbols|HGNC_Alias_names|HGNC_Chromosome|HGNC_Date_modified|HGNC_Date_symbol_changed|HGNC_Date_name_changed|HGNC_Accession_numbers|HGNC_Enzyme_IDs|HGNC_NCBI_Gene_ID|HGNC_Ensembl_gene_ID|HGNC_Pubmed_IDs|HGNC_RefSeq_IDs|HGNC_Gene_group_ID|HGNC_Gene_group_name|HGNC_CCDS_IDs|HGNC_Vega_IDs|HGNC_NCBI_Gene_ID(supplied_by_NCBI)|HGNC_OMIM_ID(supplied_by_OMIM)|HGNC_RefSeq(supplied_by_NCBI)|HGNC_UniProt_ID(supplied_by_UniProt)|HGNC_Ensembl_ID(supplied_by_Ensembl)|HGNC_UCSC_ID(supplied_by_UCSC)|Oreganno_Build|Oreganno_ID|Oreganno_Values|Simple_Uniprot_uniprot_entry_name|Simple_Uniprot_DrugBank|Simple_Uniprot_alt_uniprot_accessions|Simple_Uniprot_uniprot_accession|Simple_Uniprot_GO_Biological_Process|Simple_Uniprot_GO_Cellular_Component|Simple_Uniprot_GO_Molecular_Function|dbSNP_ASP|dbSNP_ASS|dbSNP_CAF|dbSNP_CDA|dbSNP_CFL|dbSNP_COMMON|dbSNP_DSS|dbSNP_G5|dbSNP_G5A|dbSNP_GENEINFO|dbSNP_GNO|dbSNP_HD|dbSNP_INT|dbSNP_KGPhase1|dbSNP_KGPhase3|dbSNP_LSD|dbSNP_MTP|dbSNP_MUT|dbSNP_NOC|dbSNP_NOV|dbSNP_NSF|dbSNP_NSM|dbSNP_NSN|dbSNP_OM|dbSNP_OTH|dbSNP_PM|dbSNP_PMC|dbSNP_R3|dbSNP_R5|dbSNP_REF|dbSNP_RS|dbSNP_RSPOS|dbSNP_RV|dbSNP_S3D|dbSNP_SAO|dbSNP_SLO|dbSNP_SSR|dbSNP_SYN|dbSNP_TOPMED|dbSNP_TPA|dbSNP_U3|dbSNP_U5|dbSNP_VC|dbSNP_VLD|dbSNP_VP|dbSNP_WGT|dbSNP_WTD|dbSNP_dbSNPBuildID|dbSNP_ID|dbSNP_FILTER|gnomAD_exome_AF|gnomAD_exome_AF_afr|gnomAD_exome_AF_afr_female|gnomAD_exome_AF_afr_male|gnomAD_exome_AF_amr|gnomAD_exome_AF_amr_female|gnomAD_exome_AF_amr_male|gnomAD_exome_AF_asj|gnomAD_exome_AF_asj_female|gnomAD_exome_AF_asj_male|gnomAD_exome_AF_eas|gnomAD_exome_AF_eas_female|gnomAD_exome_AF_eas_jpn|gnomAD_exome_AF_eas_kor|gnomAD_exome_AF_eas_male|gnomAD_exome_AF_eas_oea|gnomAD_exome_AF_female|gnomAD_exome_AF_fin|gnomAD_exome_AF_fin_female|gnomAD_exome_AF_fin_male|gnomAD_exome_AF_male|gnomAD_exome_AF_nfe|gnomAD_exome_AF_nfe_bgr|gnomAD_exome_AF_nfe_est|gnomAD_exome_AF_nfe_female|gnomAD_exome_AF_nfe_male|gnomAD_exome_AF_nfe_nwe|gnomAD_exome_AF_nfe_onf|gnomAD_exome_AF_nfe_seu|gnomAD_exome_AF_nfe_swe|gnomAD_exome_AF_oth|gnomAD_exome_AF_oth_female|gnomAD_exome_AF_oth_male|gnomAD_exome_AF_popmax|gnomAD_exome_AF_raw|gnomAD_exome_AF_sas|gnomAD_exome_AF_sas_female|gnomAD_exome_AF_sas_male|gnomAD_exome_ID|gnomAD_exome_FILTER|gnomAD_genome_AF|gnomAD_genome_AF_afr|gnomAD_genome_AF_afr_female|gnomAD_genome_AF_afr_male|gnomAD_genome_AF_amr|gnomAD_genome_AF_amr_female|gnomAD_genome_AF_amr_male|gnomAD_genome_AF_asj|gnomAD_genome_AF_asj_female|gnomAD_genome_AF_asj_male|gnomAD_genome_AF_eas|gnomAD_genome_AF_eas_female|gnomAD_genome_AF_eas_male|gnomAD_genome_AF_female|gnomAD_genome_AF_fin|gnomAD_genome_AF_fin_female|gnomAD_genome_AF_fin_male|gnomAD_genome_AF_male|gnomAD_genome_AF_nfe|gnomAD_genome_AF_nfe_est|gnomAD_genome_AF_nfe_female|gnomAD_genome_AF_nfe_male|gnomAD_genome_AF_nfe_nwe|gnomAD_genome_AF_nfe_onf|gnomAD_genome_AF_nfe_seu|gnomAD_genome_AF_oth|gnomAD_genome_AF_oth_female|gnomAD_genome_AF_oth_male|gnomAD_genome_AF_popmax|gnomAD_genome_AF_raw|gnomAD_genome_ID|gnomAD_genome_FILTER"
    header = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","S1"]
    skip_rows = 372 # there should be a better way to clean this
    func_df = pd.read_csv("../output/{}.variants.funconator.vcf".format(name),sep="\t",skiprows=skip_rows,names=header)
    print (func_df)
    subset_df = func_df[header[:-3]].copy()

    subset_df['FUNCOTATOR'] = func_df['INFO'].str.extract(r'FUNCOTATION=\[(.*)\];MBQ')
    subset_df['FUNC_AUX'] = subset_df['FUNCOTATOR'].str.split(",")
    subset_df['FUNC_PRIMARY'] = subset_df['FUNC_AUX'].str[0]
    func_cols_df = subset_df["FUNC_PRIMARY"].str.split(r'\|', expand=True)
    final_df = pd.concat([subset_df[header[:-3]], func_cols_df], axis=1)
    header2 = funcotator_fields.split("|")
    header1 = header[:-3]
    final_header = header1+header2
    final_df.columns = final_header
    final_df.to_csv("../output/{}.final_mutation_list.tsv".format(name),sep="\t")

    final_df["gnomAD_exome_AF"] = pd.to_numeric(final_df["gnomAD_exome_AF"],errors="coerce").astype(float)

    filtered_mutations_df = final_df.query(
        "`CGC_Tumour_Types__(Somatic_Mutations)` != '' or \
        `ClinVar_rs` != '' or \
        `gnomAD_exome_AF` <= 1 or \
        `Cosmic_overlapping_mutations` != ''"
    )

    filtered_mutations_df.to_csv("../output/{}.final_filtered_mutation_list.tsv".format(name),sep="\t")

    end = time.time()
    print (f"variants filtered for biological relevance -> Succesfull [elapsed {end - start:.4f} seconds]")

def main():

    # TASK 0
    # QC

    # TASK 1
    print ("task1")
    bam_file_name = "my_aln_sorted.bam"
    filtered_file_name = bam_QC("../sequenceFile/{}".format(bam_file_name))
    filtered_file_name = remove_under_covered_regions(filtered_file_name,min_reads=5000)
    coverage_file_name = get_genes_from_bam(filtered_file_name)
    get_gene_regions_from_bam_coverage(coverage_file_name)

    # TASK 2
    print ("task2")
    get_variants_with_mutect2(filtered_file_name)
    get_variants_relevance_funconator(filtered_file_name)
    filter_variants_for_diseases(filtered_file_name)

if __name__ == "__main__":
    main()