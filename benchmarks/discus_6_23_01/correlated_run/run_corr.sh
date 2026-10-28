#!/bin/bash
cd /mnt/c/Users/jmswe/Documents/Programming\ projects/discus-python/benchmarks/discus_6_23_01/correlated_run
/usr/local/bin/discus_suite < corr_model.mac > corr_run.log 2>&1
echo "Exit code: $?"
