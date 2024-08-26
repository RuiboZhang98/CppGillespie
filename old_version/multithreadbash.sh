#!/bin/bash

for sd in {21..40}
do
	./model5type -o toymodel5type --seed $sd --runs 100 &
done

wait
