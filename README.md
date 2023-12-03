## Getting Started
```sh
git clone https://github.com/RunpengLuo/PlsChain.git
cd PlsChain && make
# create an index for the plasmid library with k=21
./plschain -i -k 21 -o lib_idx/ backbone.fa promotor.fa peptide.fa gene.fa terminal.fa terminator.fa
# classify the reads against the indexed library
./plschain -q lib_idx/ -o qry_res/ query.fastq.gz
# summarize the classification via python script
python groupby.py qry_res/
```

## About PlsChain
PlsChain is an algorithm to classify noisy reads (~5% error rate) sequenced from the plasmid mixtures, it solves the co-linear chaining problem in cyclic manner.


## Installation
The program is designated for Unix-like system (Linux & MacOS), C compiler, GNU make and zlib development files are required to compile the program.

To run the python script `groupby.py` for grouping the results, python module pandas is required.

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
* `FILE1 FILE2 ...` the plasmid components, the order should follow the plasmid structure, cyclic order is allowed.

### Output
* `<out_dir>/qry_total.csv` stores the classification result per read
* `<out_dir>/qry_total.group.csv` stores the grouped results based on `<out_dir>/qry_total.csv`.