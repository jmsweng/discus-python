#!/bin/bash
cd "$(dirname "$0")"
/usr/local/bin/discus_suite < corr_diffuse.mac > corr_diffuse_run.log 2>&1
echo "Exit code: $?"
