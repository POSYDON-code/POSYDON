#!/bin/bash

INI_FILE="$1"
posydon-popsyn setup ${INI_FILE} --job_array=1 --walltime=00:40:00 --partition=partition --account=account --email=email@domain.com --mem_per_cpu=10G