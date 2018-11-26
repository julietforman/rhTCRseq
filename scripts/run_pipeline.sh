#!/bin/bash

cd `dirname $0`

RUN_NAME=$(python3 make_index_list.py)

python3 analyze_tcr.py --step master
wait
bash ../out/${RUN_NAME}/logs/master/master.sh
