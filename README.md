# resmarkers_utilities


## format_resmarkers.R

This R script converts mad4hatter's resmarkers_table.txt (long format) into the old "report" format. It uses the `tidyr` and `optparse` libraries for data manipulation and command-line argument parsing, respectively.

### Usage

```shell
Rscript format_resmarkers.R --input resmarkers_table.txt --output resmarkers_table_old_format.csv
