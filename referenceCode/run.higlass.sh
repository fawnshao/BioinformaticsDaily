#!/bin/bash
for f in *_allValidPairs
do
	~/Documents/GitHub/HiC-Pro/bin/utils/hicpro2higlass.sh -i $f -r 5000 -c ~/Documents/GitHub/HiC-Pro/annotation/chrom_hg19.sizes
done

docker pull gehlenborglab/higlass # Ensure that you have the latest.
docker run --detach \
           --publish 8888:80 \
           --volume ~/hg-data:/data \
           --volume ~/hg-tmp:/tmp \
           --name higlass-container \
           gehlenborglab/higlass

# COOLER=dixon2012-h1hesc-hindiii-allreps-filtered.1000kb.multires.cool 
# wget -P ~/hg-tmp https://s3.amazonaws.com/pkerp/public/$COOLER

# Confirm that the file is visible inside the container:
docker exec higlass-container ls /tmp

# Ingest:
for COOLER in SRR*.mcool
do
	docker exec higlass-container python higlass-server/manage.py ingest_tileset --filename /tmp/$COOLER --filetype cooler --datatype matrix
done
