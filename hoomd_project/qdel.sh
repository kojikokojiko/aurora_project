#!/bin/bash

# Loop from 1 to 10
for i in {746736..746747}
do
    # Run qdel command with the current number
    qdel $i
done