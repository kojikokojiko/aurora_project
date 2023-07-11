#!/bin/bash

# Loop from 1 to 10
for i in {785472..785549}
do
    # Run qdel command with the current number
    qdel $i
done