#!/bin/bash

# Loop from 1 to 10
for i in {681421..681436}
do
    # Run qdel command with the current number
    qdel $i
done