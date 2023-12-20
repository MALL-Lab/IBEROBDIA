# Study of the faecal microbiome of patients with metabolic syndrome and type 2 diabetes

## Citation

Article is open access here.

DOI: link DOI

## Cohort description
The participants of this research belong to the project IBEROBDIA: "Obesity and Diabetes in Iberoamerica: Risk factors and new pathogenic and predictive biomarkers", funded by the Iberoamerican programme of Science and Technology for Development (CyTED) (918PTE0540) and by the Spanish State Research Agency (PCI2018-093245, and PCI2018-093284). This project has the approval of the ethics committee of the Galician Health System (SERGAS, Xunta de Galicia), code 2018/270. All data were processed in accordance with the organic law 3/2018, of 5 December, on the protection of personal data and the guarantee of digital rights. A total of 79 volunteers from the autonomous community of Galicia met the inclusion/exclusion criteria of the project and participated in this study. This cohort presents 59 healthy individuals, 15 pre-diabetic, and 5 patients with T2D. Of these 59 healthy patients 36 suffer from obesity. In addition, this cohort has 24 (only 2 patients with normal weight) patients classified as having MS and 55 who do not. Of these 55 patients who do not suffer from MS 34 are obese. These volunteers underwent an 8-hour fasting glucose tolerance test using the 75 g oral glucose solution (200 mL, Lambra Glucomedics®, Lambra, Madrid, Spain) and a body composition assessment using the InBody 127 digital scale (Microcaya S. L., Bilbao, Spain). Participants were provided with a sterile anaclin for the collection of faecal samples at home. Once collected, the sample had to be transported to the laboratory within two hours, where it was kept frozen until analysis.

## Project workflow

### 00.Preprocess Metadata

[00_preprocess_metadata.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/00_preprocess_metadata/code/00_preprocess_metadata.r) This script extracts the info for FASTQ files and merge with the clinical and saves the metadata for each the cohort. This script saves a ```.rds``` with metadata.

### 01.Sequencing Data
[00_get_fastqc.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/01_sequencing_data/code/00_get_fastqc.r) This script passes a quality filter to the samples using the FASTQC program. It generates an ```.html``` file for each sample.

[01_pipeline_16S.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/01_sequencing_data/code/01_pipeline_16S.r) This script process FASTQ files following the pipeline established in the DADA2 package. Return the ASV table and the taxonomic table for each cohort. Both in ```.rds``` format.

[02_make_phy.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/01_sequencing_data/code/02_make_phy.r) In this script, the ```phyloseq``` object is built from the ```ASV table```, the ```taxonomic table``` and the ```metadata```.

### 02.Preprocess
[00_preprocess_phy](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/02_preprocess/code/00_preprocess_phy.r) This script agglomerates the phyloseq from the previous script by a given taxonomic level. It also filters out taxa and samples that do not meet minimum abundance criteria. The Phyloseq is saved in ```.rds``` format.

[01_exploratory_analysis.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/02_preprocess/code/01_exploratory_analysis.R) In this script we perform an exploratory analysis of the cohort, studying the alpha and beta diversity of the cohort.

[02_linda_analysis.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/blob/main/02_preprocess/code/02_linda_analysis.r) In this script the LinDA method  was used to measure the differences in abundance between the different groups.

### Figures
[figures](https://github.com/MALL-Machine-Learning-in-Live-Sciences/IBEROBDIA/tree/main/figures/code) In this folder are all the scripts needed to generate the figures of the paper. In the scripts 00_prepare_data.r, 01_extract_top_FamGen.r and 02_make_custom_palette.r intermediate data are generated to generate the figures.
  - Twitter: [@8JoseLinares](https://twitter.com/8JoseLinares)
  - GitHub: [jlinaresb](https://github.com/jlinaresb)
- Carlos Fernández Lozano
  - Email: <carlos.fernandez@udc.es>
  - Twitter: [@cafernandezlo](https://twitter.com/cafernandezlo)
  - GitHub: [cafernandezlo](https://github.com/cafernandezlo)
