# Origins-of-adaptive-genomic-variation-in-a-wild-rodent  

 Data  structure, scripts and input files



Set of VCF (Variant Call Format) files and txt files with filenames that were used as input for allele frequency, gene genealogy, allele age, Dxy and correlation analysis.

The output of the gene trees out of RAxML are in newick format files. The user can import those files into figtree or R to visualize the trees.

Description of the data and file structure

The files with filenames are required to run most the scripts. They are just a list of file names inside a file.

VCF files are also required to run most of the scripts. Although they are VCF files I used a flags on bcftool that allows the generation of VCF files that also contain invariant sites because I found it appropriate to include variant and invariant sites in the gene gene genealogy analysis.

A description of how VCF files are structured and formatted can be found in the following link ([https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf). This description was written by the samtools management team

##

## Code/Software.  &#x20;

## *# Scripts for gene genealogy analysis* &#x20;

All the bash scripts can be executed using the command bash scriptname.sh. The description of each script (in italic) are in the order I executed them to run the analysis.

***extract_regions_bcftools_v4.sh:*** extracts 50 kb upstream and downstream each of the 1074 outlier loci. Requires bcftools software, a file (Outliers_QvalPcadapt_Qval_pRDA_95Fst.txt) with the indexes for each outlier locus, and vcf.gz file. Generates 1074 vcf.gz files each containing 140 haplotypes of 101kb.

***phasing_script.sh:*** phases each 1074 vcf.gz file containing 140 haplotypes of 101Kb.

***bcftools_extracting_one_site_v2.sh:*** takes 1074 phased vcf.gz files with 101kb sites and extract just only the outlier locus and it stores in another vcf.gz file that contain only the alleles in outlier locus.Requires bcftools software, a file with a list of names of the phased vcf.gz file (phased_filenames_sorted_vcfgz.txt), a file (phased_outliers_coordinates_sorted_list.txt) with the indexes for each outlier locus and vcf.gz files. Generates 1074 vcf.gz files each containing only the alleles in the outlier locus.

***genotype_from_phase_onesite_v2.sh:*** Extract  the genotype information and sample ID information of a vcf file that contain a single outlier locus. Each sample ID is associated with its respective genotype for a given individual. Requires a file (onesite_vcfgz_filenames.txt) with the names of the vcf.gz files that contain only the alleles in the outlier loci.

***editing_ids_genotype_files.sh:*** remove all the characters that are not of interest from each of 1074 vcf.gz file containing only the alleles in the outlier locus. Requires a file (genotypes_and_ids_vcfgz_filenames_sorted.txt) with a list of names of vcf.gz files containing only the alleles in the outlier locus. Edits 1074 vcf.gz files in such a way that the  data is laid out following the structure ID:genotype. Those files will be used later in the pipeline.

***samtool_faidx_reference_subset_loop.sh:*** extract the regions of the reference of the genome of the bank vole that corresponds to each of the 1074 haplotypes with 101kb that contain one of the outlier locus to use as a reference to convert the phased vcf.gz with the 101kb haplotypes into fasta files in a later step of the pipeline. Requires samtools software, a file with a list of index of where each of 1074 outliers loci is within 101kb stretchs of the reference genome, and the reference genome of the bank vole in a fasta file. Generates 1074 reference fasta files that will be used to convert pashed vcf.gz files into phased fasta files.

***editing_the_header_of_referencemakere_reference_file.sh:*** remove unwanted characters from the fasta headers of the reference files. This need to be done to approprietaly convert vcf.gz files into fasta files in the later steps of this pipeline.

***samtools_index_loop.sh:*** indexes each of the 1074 reference fasta files. Requires, samtools, bwa, picard, a file (reference_fasta_files_list.txt) with a list of names of the reference fasta files.

***vcf2fasta_loop3.sh:*** converts 1074 phased vcf.gz files each containing 140 haplotypes with length of 101kb into phased fasta files. Requires vcf2fasta software, a file with a list of names of the phased vcf.gz file (phased_filenames_sorted_vcfgz.txt) and a file with a list of names of the reference fasta files (reference_fasta_files_list_sorted.txt). Generates 150,360 phased fasta files each containing a single haplotype with length of 101kb.

***merge_phased_files_loop.sh:*** merges 150,360 phased fasta files into 1074 fasta files each containing 140 phased haplotypes based on the position of the genome they are supposed to be. Requires a file (positions_scaffold_names_file_sorted.txt) with the indexes of the position of the genome the outlier loci is supposed to be.

***put_localities_names_into_fastaheader.sh:*** Transfer the locality information in a pop file (france_netherlands_slovakia_czech_heslington_lizard_maud_pop.txt) into each of the 140 fasta headers in the 1074 fasta files. Requires a pop file and a file (fasta_filenames_phased_sorted.txt) with the list of names of the phased fasta files.

***merging_genotypes_files.sh:*** transfer the ID:genotype information in 1074 vcf.gz files into a single text file (all_genotypes_ids_sorted_outliers.txt) that will be used later in the pipeline to include the genotype information into the fasta headers of the 1074 phased fasta files. Requires a file (genotypes_and_ids_vcfgz_filenames_sorted.txt) with the names of the vcf.gz files that contain only ID:genotype information.

***add_genotype_info_in_fastaheader_v14.sh:*** transfer the genotype information of the outlier loci available in the file "all_genotypes_ids_sorted_outliers.txt" to the fasta headers of the phased fasta files. The genotype are coded as 0|0 or 0|1 or 1|1 or 1|0 depending on whether the haplotypes pairs - haplotype 1 left side of the bar or haplotype 2 right side of the bar - has the reference (0) or the alternative allele (1) in the outlier loci. In the fasta headers the haplotype 1 is labeled with -0 and the haplotype 2 is labeled with -1. In the vcf.gz files the reference allele is coded as 0 and the alternative allele is coded as 1. Requires a file (all_genotypes_ids_sorted_outliers.txt) with ID:genotype information and a file with the names of the phased fasta files

***identifying_alternative_reference_in_samples_v3.sh:*** label the haplotype with the /R character if  it contains the reference allele in one of the outlier loci. Label the the haplotype with the /A character if it contains the alternative allele in one of the outlier loci. If the haplotype contains two different alternative alleles in the outlier loci, then label it with /A to represent one of the alternative alleles and label it with /E to represent the other alternative allele. Requires a file (fasta_filenames_phased_sorted.txt) with a list of names of the phased fasta files.

***edit_fasta_files.sh:*** remove characters from the fasta files that might cause problems in the folowing steps of the pipeline. Requires a list of phased fasta files.

***raxml_analysis_loop_v6.sh:*** Execute raxml in a loop with each iteration of the loop using one of the 1074 fasta files as input. Requires a raxml and a file (phased_fasta_renamed_sorted.txt) with the names of the edited phased fasta files.

***plotting_trees_rooting.R:*** Read all the 1074 newick files generated by RAxML, set the outgroup, removed unwanted charaters from tip label, give different colors for tip labels representing haplotypes that contain the reference or the alternative allele in the outlier locus.

## *# Scripts for the dxy analysis  based only on the haplotypes that contain the reference allele in the outlier loci.*&#x20;

All the bash scripts can be executed using the command bash scriptname.sh. The description of each script (in italic) are in the order I executed them to run the analysis.

***reference_separate_vcf_samples.sh:*** take 1074 vcf.gz each containing 140 phased haplotypes and generate 75,180 vcf files each containing 2 phased haplotypes. Requires bcftools, a file (vcffilenames.txt) with the names of the phased  vcf.gz files, a file with sample ids (sample_ids.txt).

***reference_creating_filenames_file.sh:*** create a file (reference_vcf_by_sample_filenames.txt) with the names of the 75,180 vcf files that contains a pair of haplotypes. Requires a file (bank_vole_loci_under_selection_chrom_positions.txt) with the indexes for the position of the genome where the outlier loci are supposed to be.

***reference_filenames_indexed_by_genotypes.sh:*** create a file (reference_genotypes_and_separate_vcffilename.txt) with the file names of the 75,180 vcf files each containing 2 phased haplotypes and index those files names with the genotype that is found in the outlier locus of the haplotypes. Requires (reference_vcf_by_sample_filenames.txt) a file with the names of the 75,180 vcf files

***get_reference_haplotypes_from_vcf.sh:*** exclude the haplotypes that have the altenartive allele in the outlier locus resulting in vcf files that contain only the haplotypes that have the reference allele in the outlier loci. Requires a file (reference_genotypes_and_separate_vcffilename.txt) with the file names of the 75,180 vcf files each containing 2 phased haplotypes and the genotype in the outlier loci of those haplotypes.

***reference_vcf_compression_loop.sh:*** compress the 75,180 vcf files into 75,180 vcf.gz files so we can merge them into 1074 vcf.gz files later in the pipeline. Requires a file (reference_genotypes_and_separate_vcffilename.txt) with the file names of the 75,180 vcf files each containing 2 phased haplotypes and the genotype in the outlier loci of those haplotypes. It also requires bcftools software.

***reference_index_loop_v2.sh:*** indexes 75,180 vcf.gz files so we can merge them into 1074 vcf.gz files later in the pipeline. Requires a file (reference_genotypes_and_separate_vcffilename.txt) with the file names of the 75,180 vcf files each containing 2 phased haplotypes and the genotype in the outlier loci of those haplotypes.

***reference_creating_filenames_by_position.sh:*** create txt files containing the names that will be necessary to merge the 75,180 vcf.gz files into 1074 .vcf.gz files later in the pipeline. Requires a file (vcfgz_filenames.txt) with the names of 1074 vcf.gz files

***reference_merging_vcfiles_by_position_and_samples.sh:*** merges the 75,180 vcf.gz files each containing only the haplotypes that have the reference allele in the outlier loci into 1074 vcf.gz files according to the position of the genome where the outlier loci is supposed to be. Requires a file (reference_concatenated_filenames.txt) that contains the information about where each haplotype is supposed to be in the genome.

***reference_vcf_uncompression_loop.sh:*** after the last merging step, characters that will cause problems later down the pipeline are introduced  in the vcf.gz files. This script uncompress 1074 vcf.gz files into 1074 vcf files so it is possible to replace the problematic characters with the appropriate characters in the the next step of the pipeline. Requires a file (reference_concatenated_filenames.txt) with the names of the vcf.gz files.

***reference_replacing_characters_loop.sh:*** replace the unwanted characters with the appropriate characters in 1074 vcf files that will make it possible to run the pixy software using a compressed and indexed version of these VCF files in the later steps of the pipeline. Requires a file (reference_concatenated_filenames.txt) with the names of the vcf files.

***reference_vcf_compression_loop.sh:*** compress 1074 edited vcf files into 1074 vcf.gz files. Requires a file (reference_concatenated_filenames.txt) with the names of the vcf files.

***reference_index_loop_v2.sh:*** indexes 1074 vcf.gz files. Requires a file (reference_concatenated_filenames.txt) with the names of the vcf files.

***tyding_pxy_analysis_reference_v2.sh:*** run the pixy software on 1074 vcf.gz files that contain only the haplotypes with the reference allele in the outlier loci. Requires pixy, a file  (reference_concatenated_filenames.txt) with the names of the vcf.gz files.

\**# To run the pixy software using only the haplotypes containing the alternative allele in the outlier loci repeat all the steps above using the files with the label alternative in them.    **

## *#Scripts to execute the allele frequency analysis*&#x20;

All the scripts can be executed as bash scriptname.sh. The description of each script (in italic) are in the order I executed them to run the analysis.

***bcftools_extracting_one_site_for_plink.sh:*** extract only the outlier loci from the whole genome and ouputs each loci on a separate vcf.gz file. Requires bcftools and a  file (Outliers_QvalPcadapt_Qval_pRDA_95Fst.txt) with the indexes of the positions where the outliers are in the genome.

***concat_vcf_bcf_files.sh:*** merge the separate vcf.gz files containing the outlier locus into a single vcf.gz that contain all the outlier loci. Requires bcftools and a file (onesite_plink_list.txt) containing the names of the separate vcf.gz files each containg a single outlier loci.

***plink2_analysis_v2.sh:*** executes plink to calculate the allele frequencies in the outlier loci. Requires a pop file (france_netherlands_slovakia_czech_heslington_lizard_maud_pop_for_plink.txt), a vcf.gz file with the outlier loci.
