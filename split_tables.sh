#!/bin/bash
mkdir -p table_parts
ls *txt | xargs -I{} -t sh -c ' FILE="{}"; iconv -f utf-16 -t utf-8 "{}" | grep -B 8 "~End" | sed "s/^\s*//gi" | sed "s/\s*$//gi" | grep -v End | grep -v "^--$" | split -l 8 --numeric-suffixes=1 --additional-suffix="_${FILE%.*}.tab" - "table_parts/"'
