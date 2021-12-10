#!/bin/bash -u

mkdir -p gather

cd gather

for i in `find .. -type f -name \*.tar.gz`
do
  ln -sf $i .
done

cd - > /dev/null

for i in `find gather -type l`
do
  filename=$(basename $i)
  YEAR=${filename%_*}
  YEAR=${YEAR#*_}
  MODEL=${filename##*_}
  MODEL=${MODEL/.tar.gz/}
  mkdir -p data/$YEAR/$MODEL
  tar -zxvkf $i --directory data/$YEAR/$MODEL
done