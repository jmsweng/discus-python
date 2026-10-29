#!/bin/bash
cd /mnt/c/Users/jmswe/Documents/Programming\ projects/discus-python/benchmarks/discus_6_23_01/short_range_run
/usr/local/bin/discus_suite < short_range_mac.mac > sr_run.log 2>&1
echo "Exit code: $?"
