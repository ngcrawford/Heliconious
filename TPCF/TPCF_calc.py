import os
import sys
import glob
import time
import warnings
import argparse
import numpy as np
import random as rd
import multiprocessing 
from scipy import stats

# IGNORE RUNTIME WARNINGS DUE TO LOG ZERO
warnings.simplefilter("ignore", RuntimeWarning) 

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    
    parser.add_argument('input', 
                        help='Path to island file.')

    parser.add_argument('-c','--chromosome',
                        required=False,
                        default=None,
                        type=str,
                        help='Filter on selected chromosome.')

    parser.add_argument('-b', '--bootstraps',
                        required=False, 
                        default=0,
                        type=int,
                        help='Number of bootstraps to run.')
   
    parser.add_argument('-m','--multiplier',
                        required=True,
                        type=int, 
                        help='Multiplier value to determine \
                        the number of random points to generate.')

    parser.add_argument('-p','--processor-cores',
                        # auto detect the number of availble CPUs
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help="Defaults to the maximum number \
                        of available cores on the processor")
    
    args = parser.parse_args()
    return args


def round_down_bins(x):
    return x if x % 0.25 == 0 else x - (x % 0.25)

round_down_bins = np.vectorize(round_down_bins)


def calc_bins(diffs):
    """Uses numpy builtin fuctions to caculate bins."""

    logs = np.log10(np.abs(diffs))              # log transform diffs
    logs[np.isinf(logs)] = 0                    # -infs in diagonal to zero
    binned_logs = round_down_bins(logs).astype(float)    # bin values
    bin_ids, inv = np.unique(binned_logs.flatten(), return_inverse=True) # count bins
    counts = np.bincount(inv)
    # bin_ids = np.nonzero(counts)[0]             # make bin ids
    bin_counts = zip(bin_ids,counts)   # zip it up
    return bin_counts

def calc_TPCF(pts, rand_pts, N, M):

    # USE RESHAPE TO MAKE THE MATRICES OF DIFFERENCES
    DDr_diffs = pts - pts.reshape(len(pts),1)       
    DRr_diffs = rand_pts - pts.reshape(len(pts),1)
    RRr_diffs = rand_pts - rand_pts.reshape(len(rand_pts),1)

    # CACULATE DDr, DRr, AND RRr
    DDr = dict([(k, (1.0/N*(N-1))*v) for k,v in calc_bins(DDr_diffs)])
    DRr = dict([(k, (1.0/N*M)*v)     for k,v in calc_bins(DRr_diffs)])
    RRr = dict([(k, (1.0/M*(M-1))*v) for k,v in calc_bins(RRr_diffs)])

    # CACULATE 1 + SQUIGGLY
    for r in RRr:
        if r in DDr.keys() and r in DRr.keys(): # make sure the bin is present in all comparisions
            squig = ((DDr[r] - 2*DRr[r] + 2*RRr[r]) / RRr[r])
            yield r, squig

def bootstrap(d):
    """Bootstrap function."""

    pts, reps, multiplier, chrm_name, chrm_len = d
    N = float(len(pts))
    boot_values = {}
    
    start_time = time.time() # replace with timeit?
    while reps != 0:

        pts = np.array([rd.choice(pts) for p in xrange(len(pts))])              # sample with replacement the pts
        rand_pts = np.random.randint(1, chrm_len, len(pts)*multiplier)          # create new set of random pts

        # Do TPCF calculations
        M = float(len(rand_pts))
        for r, squig in calc_TPCF(pts, rand_pts, N, M):
            if boot_values.has_key(r) == False:
                boot_values[r] = [squig]
            else: 
                boot_values[r].append(squig)

        elapsed_time = time.time() - start_time
        #print "{} sec, bootstrap iteration {}".format(elapsed_time, reps)
        reps -= 1

    # summarize results (bin, mean, stdev, sterr)
    results = []
    for r,values in boot_values.iteritems():
        results.append([chrm_name, r, len(values), np.mean(values), np.std(values), stats.sem(values)])
    return results 


def parse_islands(fin):
    fin = open(fin, 'rU')
    chrms = {}
    for line in fin:
        line_parts = line.strip().split("\t")
        chrm, start, stop = line_parts
        start = float(start)
        if chrms.has_key(chrm) == False:
            chrms[chrm] = [start]
        else:
            chrms[chrm].append(start)

    fin.close()
    return chrms 

def main(args):

    chrm_sizes = {'chr1':       15890349.0,
        'chr10':      17598034.0,
        'chr11':      11318171.0,
        'chr11_unmapped':   322530.0,
        'chr12':      15942582.0,
        'chr12_unmapped':   21850.0,
        'chr13':      14025408.0,
        'chr13_unmapped':   328782.0,
        'chr14':      6768282.0,
        'chr15':      8075189.0,
        'chr15_unmapped':   119851.0,
        'chr16':      9422647.0,
        'chr16_unmapped':   138787.0,
        'chr17':      14112702.0,
        'chr17_unmapped':   47415.0,
        'chr18':      15591294.0,
        'chr18_unmapped':   197652.0,
        'chr19':      15084489.0,
        'chr19_unmapped':   364302.0,
        'chr1_unmapped':  396482.0,
        'chr2':       3608632.0,
        'chr20':      5856246.0,
        'chr20_unmapped':   887595.0,
        'chr2_unmapped':  93281.0,
        'chr3':       9185411.0,
        'chr3_unmapped':  180637.0,
        'chr4':       6784399.0,
        'chr5':       8121235.0,
        'chr5_unmapped':  127335.0,
        'chr6':       13324770.0,
        'chr6_unmapped':  183240.0,
        'chr7':       12031174.0,
        'chr7_unmapped':  334715.0,
        'chr8':       7067414.0,
        'chr8_unmapped':  304413.0,
        'chr9':       8360004.0,
        'chrZ':       4131199.0,
        'chrZ_unmapped':  479346.0}


    print 'Pair', 'Chrm','Bin','Count','Mean','STDev','SEM'
    for fin in glob.glob(args.input + "/*.txt"):    
        chrms = parse_islands(fin)

        data = []
        #print 'Summary:'
        for chrm in chrms.keys():

            if 'unmapped' in chrm:   # skip unmapped chromosomes
                continue
            
            if args.chromosome != None and args.chromosome != chrm:
                    continue

            pts = chrms[chrm]       
            
            if len(pts) == 1:    # skip chromosomes with singleton islands
                continue

            chrm_info = (chrms[chrm], args.bootstraps, args.multiplier, chrm, chrm_sizes[chrm])
            data.append(chrm_info)

        # print '\nResults:'
        
        pool = multiprocessing.Pool(args.processor_cores)
        results = map(bootstrap, data) 
        #results = pool.map(bootstrap, data)    
        for j in results:
            for r in j:
                print os.path.split(fin)[-1].split('.')[0], " ".join(map(str, r))


if __name__ == '__main__':
    args = get_args()
    main(args)


