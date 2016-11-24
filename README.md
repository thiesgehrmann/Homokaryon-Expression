# Homokaryon expression
The discovery and analysis of homokaryon specific expression
![The General Overview of Homokaryon Specific Expression](figure0.png)
## Usage

### Configuration
  In the /pipeline.sh/ file, change the variable ARALOC, to match your own setup.

### Running the example dataset
  By running ./run_example.sh, the example dataset will be run. Output will be produced in output_example

### Real data

  Additionally, if you are willing dow download your own data, there are two configuration files provided for two publicly available datasets:
  * agabi_compost_nodup.tsv: https://www.ncbi.nlm.nih.gov/bioproject/309475
  * agabi_compost_data.tsv: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA275107

These can be run with the script "pipeline.sh"


## Dependencies
This pipeline makes use of
  * Ibidas (https://github.com/thiesgehrmann/ara)
  * Numpy
  * Matplotlib
  * 
  * ara (https://github.com/thiesgehrmann/ara)
  * R
  * DESeq

