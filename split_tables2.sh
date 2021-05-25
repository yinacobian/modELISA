#!/bin/bash
mkdir -p table_parts
ls *txt | xargs -I{} -t sh -c ' FILE="{}"; iconv -f utf-16 -t utf-8 "{}" | grep -A 8 "Temperature" | cut -f3-15 | grep -v "^--$" | sed "s/\s*$//gi" | grep -v "^--$" | split -l 9 --numeric-suffixes=1 --additional-suffix="_${FILE%.*}.tab" - "table_parts/"'
