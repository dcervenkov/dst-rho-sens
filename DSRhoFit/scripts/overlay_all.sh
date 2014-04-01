#!/bin/bash

min=$1
max=$2

let 'count = max - min + 1'
let 'i = min - 1'
counter=0

while [ $i -lt $max ]; do
    let 'i = i + 1'
    let 'counter = counter +1'
    echo "Processing $counter/$count ..."
    ./overlay plots/projections_evt$i.root plots/projections_gen$i.root plots/projections_b_$i.root
done
