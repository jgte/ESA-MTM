#!/bin/bash -u

URL="ftp://ig2-dmz.gfz-potsdam.de/ESAESM/mtmshc/"

FLAGS="--continue --recursive -A *.tar.gz,index.html --level=2"

wget $FLAGS $URL

#removing index.html files
find . -name \*index.html -delete

#pruning empty dirs
find . -type d -empty -delete