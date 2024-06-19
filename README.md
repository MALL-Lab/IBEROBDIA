# Study of the faecal microbiome of patients with metabolic syndrome and type 2 diabetes

## Authors

Laura Isabel Sinisterra Loaiza<sup>1†</sup>, Diego Fernandez-Edreira<sup>2†</sup>, Jose Liñares-Blanco<sup>2</sup>, Alberto Cepeda<sup>1</sup>, Alejandra Cardelle-Cobas<sup>1</sup>, Carlos Fernandez-Lozano<sup>2*</sup>

<sup>1</sup>Machine Learning in Life Sciences Lab. Dept. of Computer Science and Information Technologies, Universidade da Coruña (CITIC), A Coruña, Spain.

<sup>2</sup>Department of Analytical Chemistry, Nutrition and Bromatology. Faculty of Veterinary Medicine. Universidade de Santiago de Compostela, Campus de Lugo. Lugo, Spain.

<sup>*</sup>Corresponding author(s). E-mail(s): carlos.fernandez@udc.es;

Contributing authors: laura.sinisterra@usc.es; diego.fedreira@udc.es; j.linares@udc.es; alberto.cepeda@usc.es; alejandra.cardelle@usc.es;

<sup>†</sup>These authors contributed equally to this work.

## Citation

Article is open access here.

DOI: link DOI

## Recruitment and Participants
The participants of this research belong to the project IBEROBDIA: “Obesity and Diabetes in Iberoamerica: Risk factors and new pathogenic and predictive biomarkers”, funded by the Iberoamerican Program of Science and Technology (CyTED) (918PTE0540) and by the Spanish State Research Agency (PCI2018-093284). This project has the approval of the ethics committee of the Galician Health System (SERGAS, Xunta de Galicia), code 2018/270. All the participants included in the study were adequately informed about the process and signed an informed consent. Data from volunteers were processed in accordance with the organic law 3/2018, of 5 December, on the protection of personal data and the guarantee of digital rights.
The participants in this study were recruited through the dissemination of the project in different media, press, radio, posters, etc. in different locations of from the autonomous community of Galicia, Spain.
Inclusion criteria for participation in the project included: i) Spanish adults aged 40 to 70 years, ii) Body Mass Index (BMI) of $18.5-24.9$ (normal weight), or $IMC \geq 27 kg/m2$, iii) FINDRISC SCORE $\geq 12$, indicating a moderate of T2D development.
Exclusion criteria were: i) antibiotic and probiotic treatment in the past two months, ii) diagnostic of T2D, iii) presence of other chronic diseases, iii) pregnancy iv) antibiotic consumption in the last to 2 months prior the study, v) chronic medication (hypertension, cholesterol, contraceptives, proton pump inhibitors, etc., v) drug consumption and vi) low risk of developing diabetes according the FINDRISC test.

## Funding
The authors thank the CyTED, Spain and each National Organism for Science and Technology for funding the IBEROBDIA project (P918PTE0409). In this regard, Spain specifically thanks the Ministry of Economy and Competitiveness for the financial support for this project through the State Program of I+D+I Oriented to the Challenges of Society 2017–2020 (International Joint Programming 2018), project (PCI2018-093245).

CITIC is funded by the Xunta de Galicia through the collaboration agreement between the Ministry of Culture, Education, Vocational Training, and Universities and the Galician universities for the strengthening of research centers in the University System of Galicia (CIGUS). The authors acknowledge the support of CESGA (Centro de Supercomputación de Galicia) for providing computing resources and related technical support that contributed to the research results reported in this paper. 

JLB work was financed by the Spanish Ministry of Universities by means of the Margarita Salas (RSUC.UDC.MS06) linked to the European Union through the NextGenerationEU program.

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
