#!/bin/bash
cd "$(dirname "$0")"
/usr/local/bin/discus_suite < corr_model.mac > corr_run.log 2>&1
echo "Exit code: $?"
