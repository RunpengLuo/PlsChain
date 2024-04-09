## Getting Started
```sh
git clone https://github.com/RunpengLuo/PlsChain.git
cd PlsChain && make
# create an index for the plasmid library with k=15
./plschain -i -k 15 -o lib_idx/ backbone.fa promotor.fa peptide.fa gene.fa terminal.fa terminator.fa
# classify the reads against the indexed library
./plschain -q lib_idx/ -o qry_res/ query.fastq.gz
# perform fuzzy match and group the classification
python scripts/plschain_postprocess.py qry_res/ lib_idx/
```

## About PlsChain
PlsChain is an algorithm to classify Oxford Nanopore noisy reads (~5% error rate) sequenced from the plasmid mixtures, it solves the cyclic co-linear chaining problem in the cyclic manner.


## Installation
The program is designated for Unix-like system (Linux & MacOS), C compiler, GNU make and zlib development files are required to compile the program.

Run the python script `scripts/plschain_postprocess.py` for grouping the results with a Python3 environment with no additional library been required.

## Program Usage
```sh
Usage: plschain -i -k INT -o DIRECTORY FILE1 FILE2 FILE3 ...
       plschain -q DIRECTORY -o DIRECTORY <query.fa>
Options:
    -i            Indexing mode
    -q DIRECTORY  Query mode, index directory
    -k INT        k-mer size [15,32]
    -o DIRECTORY  output directory
    -h            show this message
```
* `FILE1 FILE2 ...` consists the plasmid components, the order should follow the plasmid structure, cyclic order is allowed, e.g., `backbone.fa promotor.fa peptide.fa gene.fa terminal.fa terminator.fa`.

```sh
$python scripts/plschain_postprocess.py
scripts/plschain_postprocess.py <query_dir> <index_dir>
```
* `index_dir`refers to the output directory after running PlsChain with `-i` indexing mode, and `query_dir` refers to the output directory after running PlsChain with `-q` query mode.

### Output
* `<out_dir>/qry_total.csv` and `<out_dir>/qry_total.fuzzy.csv` stores the classification result per read with and without fuzzy match opertaions. Each row consists read name, followed by the ordered list of classified components. `*` indicates the corresponding component is not decided by PlsChain. `fail` indicates unclassified record. `contamination` indicates the filtered unclassified record as contamination based on read length.

* `<out_dir>/qry_total.group.csv` and `<out_dir>/qry_total.group.fuzzy.csv` stores the grouped results based on `<out_dir>/qry_total.csv` and `<out_dir>/qry_total.fuzzy.csv`, respectively.