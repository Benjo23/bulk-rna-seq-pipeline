# bulk-rna-seq-pipeline
For Nextflow pipeline

Steps to set up a new run

Create working directory and retrieve files from Github using this command:
git clone https://github.com/Biomedical-Genetics/bulk-rna-seq-pipeline/tree/main/nextflow

Generate params.infile in TSV format using create_params.infile.R and the example of TSV file can be found in params.infile_example.tsv.

Edit RNA_Seq.config file:

params.infile: full path to the tsv file
params.output_dir: full path to working directory created in step 1
params.prefix: prefix use to label output files
params.genome: set genomic parameters which include species, ucsc, assembly, set and ensemble
params.read_length: read length used
params.paired_end: true for paired end; false for unpaired end
params.stranded: true for stranded while false for not stranded
params.modules: change the module versions used if it is necessary
Start the pipelines using submit_RNA_Seq.qsub.
This will run Nextflow script RNA_Seq.nf using params.infile_example.tsv and RNA_Seq.config.

Description of RNA_Seq.nf

This Nextflow pipeline contains the processes as shown as https://github.com/compbiomed/RNA_Seq except that generateFASTA , runRSeQCsexcheck and runOutliersdetection are included as additional steps.

Note:

Add // in front of the line to skip running the lines.
Adjust the memory and time as needed in RNA_Seq.config file.

Future direction

Create singularity container to improve reproducibility and minimize the errors.
