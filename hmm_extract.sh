#!/bin/bash

grep -v '^#' data/pfam.tblout | awk '{print $1"\t"$2"\t"$3"\t"$19" "$20" "$21" "$22" "$23" "$24" "$25" "$26" "$27" "$28" "$29" "$30" "$31" "$32" "$33" "$34" "$35" "$36}' > data/pfam.tsv