#!/bin/bash

# Loop from 1 to 10
for i in {585469..585484}
do
    # Run qdel command with the current number
    qdel $i
done