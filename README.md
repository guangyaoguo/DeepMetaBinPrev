
---

<div align="center">

#  A graph-based Gaussian Mixture Variational Autoencoder for Metagenome Binning

<a href="https://pytorch.org/get-started/locally/"><img alt="PyTorch" src="https://img.shields.io/badge/PyTorch-ee4c2c?logo=pytorch&logoColor=white"></a>
<a href="https://pytorchlightning.ai/"><img alt="Lightning" src="https://img.shields.io/badge/-Lightning-792ee5?logo=pytorchlightning&logoColor=white"></a>
<a href="https://hydra.cc/"><img alt="Config: Hydra" src="https://img.shields.io/badge/Config-Hydra-89b8cd"></a>
<a href="https://github.com/ashleve/lightning-hydra-template"><img alt="Template" src="https://img.shields.io/badge/-Lightning--Hydra--Template-017F2F?style=flat&logo=github&labelColor=gray"></a><br>
</div>

## Description

This repo contains first implementations of deep learning method of metagenomic binning so-called DeepMetaBin,

The code repository is organized into the following components:
| Component | Description |
| --- | --- |
| [datamodules](https://github.com/mx-ethan-rao/deepmetabin/tree/multi_sample_finished/src/datamodules) | Contains torch dataset objects and pl.LightningDataModules for Graph embeded GMVAE in DeepMetaBin|
| [models](https://github.com/mx-ethan-rao/deepmetabin/tree/multi_sample_finished/src/models) | Contains torch module objects and pl.LightningModules for GMVAE, deepbin |
| [utils](https://github.com/mx-ethan-rao/deepmetabin/tree/multi_sample_finished/src/utils) | Contains util functions in the project, shared visualization and evaluation functions across different model backbones. |
| [configs](https://github.com/mx-ethan-rao/deepmetabin/tree/multi_sample_finished/configs) | Contains hydra based config files to control the experiments across differernt models. |

<p align="center"><img src="./pictures/Poster_Recomb_2023.png" width=90% height=50%></p>
<p align="center"><em>Figure.</em> Poster Recomb 2023.</p>



## To Do List:
- :black_square_button: Update the Project Page of DeepMetaBin.
- :black_square_button: Released the arXiv version of DeepMetaBin (Under review).
- :white_check_mark: Released the early version of sample code.

## DeepMetaBin was tested using the following software
* [Gurobi==10.0.1](https://anaconda.org/Gurobi/gurobi) (need to request for license)
* [bwa==0.7.17](https://github.com/lh3/bwa)
* [samtools==1.9](https://github.com/samtools/samtools)
* [checkM==1.1.2](https://github.com/Ecogenomics/CheckM)
* [HMMER==3.3.2](http://hmmer.org/)
* [Prodigal==2.6.3 ](https://github.com/hyattpd/Prodigal)
* [FragGeneScan==1.31 (optional)](https://sourceforge.net/projects/fraggenescan/)

## Examples and experiments run on a Linux server with the following specifications:
* Dell PowerEdge R6525
* CPU: Dual 64-core AMD EPYC 7742 2.25GHz 256MB L3 cache
* Memory: 1T
* Operating System: Oracle Linux 8.7 (64-bit)

## Installation (estimated time: 10min)

Install dependencies

```bash
# clone project
git clone https://github.com/mx-ethan-rao/deepmetabin.git
cd deepmetabin

# [OPTIONAL] create conda environment
conda create -n myenv python=3.9
conda activate myenv

# install pytorch according to instructions
# https://pytorch.org/get-started/

# install requirements
pip install -r requirements.txt
```
## Single sample
### Preprocess Dataset From Scratch
```bash
bwa index my_contigs.fna 
bwa mem contigs.fasta reads_file_1 reads_file_2 -t 100 | samtools sort -@ 40 -n -o contigs.map.sorted.bam
python preprocessing.py --outdir /path/to/out --fasta my_contigs.fna --bamfiles *.bam -m 1000
```

### Preprocess Dataset from sample data (estimated running time for samle data: 45min)
Please download the sample dataset [sample_data](https://drive.google.com/drive/folders/1G8Wlws3HT4BtrBWG_KZk8w755MICpQXI?usp=sharing) and put it under deepmetabin folder
```bash
python preprocessing.py --outdir sample_data --fasta ./sample_data/contigs.fasta --bamfiles ./sample_data/contigs.map.sorted.bam -m 1000
```

### Run for DeepMetaBin

```bash
python run.py datamodule.zarr_dataset_path=sample_data/data.zarr datamodule.output=deepmetabin_out model.contignames_path=sample_data/contignames.npz model.contig_path=sample_data/contigs.fasta
```
You can specify the trainning epochs, e.g trainer.max_epochs=500
```bash
cd ./deepmetabin_out/results/pre_bins && checkm lineage_wf -t 100 -x fasta --tab_table -f checkm.tsv ./ ./
```
```bash
python secondary_clustering.py --primary_out ./deepmetabin_out --contigname_path ./sample_data/contignames.npz --contig_path ./sample_data/contigs.fasta --output_path  ./deepmetabin_out/results --binned_length 1000 
```

The binning result is under ./deepmetabin_out

## Multiple sample
### Preprocess Dataset From Scratch
```bash
python concatenate.py /path/to/catalogue.fna.gz /path/to/assemblies/sample1/contigs.fasta /path/to/assemblies/sample2/contigs.fasta  [ ... ]
bwa index /path/to/catalogue.fna.gz
bwa mem -p /path/to/reads_file_1 /path/to/catalogue.fna.gz -t 100 | samtools sort -@ 40 -n -o /path/to/output_bamfile
```

```bash
python preprocessing.py --outdir out --fasta /path/to/catalogue.fna.gz --bamfiles /path/to/bamfiles/*.bam --outdir /path/to/out -m 1000
python ./src/utils/split_samples.py \
        *contigs.fast \
        --concat_tnf ./out/tnf.npz \
        --concat_rkpm ./out/rpkm.npz \
        --concat_contignames ./out/contignames.npz \
        --out ./out/multisample_out
```

### Sample run for DeepMetaBin
```bash
bash ./bash/run_multisample.sh multisample_name ./out/multisample_out/tnf.npz ./out/multisample_out/rkpm.npz ./out/multisample_out/ontignames.npz ./out/multisample_out
```
The binning result is under ./out/multisample_out


## Reference Code

- https://github.com/ashleve/lightning-hydra-template
- https://github.com/samtools/samtools
- https://github.com/lh3/bwa

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
