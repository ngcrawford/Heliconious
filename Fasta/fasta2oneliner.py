#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Nick Crawford on 2011-11-18.
Copyright (c) 2011

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""
import os
import sys
import glob
import argparse
import numpy as np
from pyfasta import Fasta

fastas = glob.glob("/Volumes/PROMISE_PEGASUS/Heliconius_Genome_Reseq_Data/MZH11020_Butterfly_re-sequencing/result_variation/snp/snp_result/*.fa")

window_size = 10000

for fa_count, fa in enumerate(fastas):
    fa = Fasta(fa)
    fa = np.array(fa["1"])
    slices = range(0,len(fa),window_size)
    print len(slices)
    
    for count, slice in enumerate(slices):
        if count == len(slices) -1: break
        start = slice
        stop = slices[count + 1]
        print start, stop, count
        print fa[start:stop]
        
    if fa_count == 0: break
