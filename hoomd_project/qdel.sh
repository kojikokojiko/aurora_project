#!/bin/bash

# Loop from 1 to 10
for i in {784269..785581}
do
    # Run qdel command with the current number
    qdel $i
done