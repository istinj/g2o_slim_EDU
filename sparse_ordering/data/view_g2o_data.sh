#!/bin/bash
project_dir="$(pwd)"
file="$1"

fullpath=${project_dir}/${file}

cd /usr/include/g2o/bin
./g2o_viewer ${fullpath}
