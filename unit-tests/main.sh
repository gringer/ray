#!/bin/bash

(

echo "Nothing reported = all unit tests passed"
echo ""

for i in $(ls test_*.sh)
do
	echo "Test Suite: $i"
	bash $i
done

) | tee UnitTests.txt

echo "Wrote results in UnitTests.txt"
