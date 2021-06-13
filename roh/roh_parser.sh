#!/bin/bash

awk 'NR == 1 {print $0 "\t sample_id"; next;}{print $0 "\t" FILENAME;}' 