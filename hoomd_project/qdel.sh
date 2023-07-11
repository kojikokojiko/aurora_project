#!/bin/bash

# Loop from 1 to 10
for i in {776427..776455}
do
    # Run qdel command with the current number
    qdel $i
done