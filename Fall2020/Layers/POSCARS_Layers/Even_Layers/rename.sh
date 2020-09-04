#!/bin/bash

for i in {2..20}
do
	mv POSCAR_${i}L POSCAR_L$i
done

