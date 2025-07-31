#!/bin/bash
#

BASE_URL="https://opm-assets.storage.googleapis.com/pdb"
out="$HOME/Downloads/pdbs"

while IFS= read -r line; do
    echo "Downloading from OPM database: $line.pdb"
    url="$BASE_URL/$line.pdb"
    wget $url -P $out -q
done < pdbid
