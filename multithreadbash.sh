#!/bin/bash

for sd in {21..40}
do
	./Gillespie812 -o model812 --seed $sd --runs 50 &
done

wait
