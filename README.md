# MetaG: A Comprehensive Metagenomic Analysis Tool

This application is a metagenomic visualization platform that is created using Dash Framework.

MetaG, a standalone metagenomic analysis tool helps users in converting higher dimensional sequence data into lower dimensional visualizations through deep learning techniques. MetaG is uniquely designed to produce rich lower dimensional representations of oligonucleotide frequency vectors of the population-level genomic diversity of the microbial organisms in a given metagenome sample.

## Getting Started

### Using MetaG

To get started, upload a metagenomic reads fasta file you want to obtain visualizations of. Scatter plots and other related charts appears on each tab of MetaG one you insert the relevant parameters.



### Running the app locally

First create a virtual environment with conda or venv inside a temp folder, then activate it.

```
# With conda:
conda create -n dash-ma-venv python=3.x anaconda

activate dash-ma-venv

# With venv:
virtualenv dash-ma-venv

# Windows
dash-ma-venv\Scripts\activate
# Or Linux
source dash-ma-venv/bin/activate
```

Clone this git repo, then install the requirements with pip

```
git clone https://github.com/kaymers/MetaG
cd MetaG

pip install -r requirements.txt
```

Configure the Backend
```
Clone https://github.com/kaymers/MetaG_backend into THIS directory (root directory of this repository). (.i.e resulting file structure would be MetaG[this directory]/MetaG_backend[backend])

Follow the instructions in the https://github.com/kaymers/MetaG_backend repository.
```


Run the app

```
python app.py
```


