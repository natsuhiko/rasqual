#!/bin/bash

grep -v "^#" usage.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "#include<stdio.h>\nvoid usage(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' > usage.c
