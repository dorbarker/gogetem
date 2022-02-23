# gogetem

Get nucleotide data from EMBL for Gene Ontology terms

## Install

### pip
```commandline
git clone https://github.com/dorbarker/gogetem.git

python -m pip install gogetem
```

### conda
```commandline
conda install -c conda-forge -c dorbarker gogetem 
```

## Typical invocation
```commandline
gogetem --go-terms 0003677 --limit 2000 --download-path results/ --include-amino-acids
```