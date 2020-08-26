# METAGENOMICS ANALYZER

This application is a visualization platform that is created using Dash Framework.

This application will display the the plots the given metagenomic sequences with various options for visualizations. For the dimensionality reduction Autoencoders will be used when processing the data.

## Getting Started

### Using METAGENOMICS ANALYZER

To get started, upload a metagenomic .fasta file you want to visualize. When the scatter plots and other charts appears on the first tab.



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
git clone 
cd MetaG

pip install -r requirements.txt
```

Configure the Backend
```
Clone https://github.com/kaymers/MetaG_backend into THIS directory. (.i.e resulting file structure would be MetaG[this directory]/MetaG_backend[backend])

Follow the instructions in the https://github.com/kaymers/MetaG_backend repository.
```


Run the app

```
python app.py
```


