#!/bin/bash
cd "$(dirname "$0")"
/usr/local/bin/discus_suite < short_range_mac.mac > sr_run.log 2>&1
echo "Exit code: $?"
