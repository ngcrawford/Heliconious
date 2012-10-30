from pylab import *
import scipy.stats

def getlnLs(fin):
    lnLs = []
    bps = []
    for line in open(fin,'rU'):
        parts = line.split("=")[0].strip("tree ").split(",")
        parts_dict = dict([item.split(":") for item in parts])
        bps.append(int(parts_dict['start']))
        lnLs.append(float(parts_dict['lnL']))
    bp_lnLs = array(zip(bps,lnLs))
    sorted_bp_lnLs = bp_lnLs[bp_lnLs[:,0].argsort()]
    return sorted_bp_lnLs


unconstrained = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/trees/32.butterflies.1000bootstraps.trimmed.10000bp.lnL.trees"
constrained = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice/trees/32.butterflies.1000bootstraps.trimmed.10000bp.lnL.constraint.trees"
unconstrained_lnLs = getlnLs(unconstrained)
constrained_lnLs = getlnLs(constrained)


# bp_lnLs = array(zip(bps,lnLs))
# sorted_bp_lnLs = bp_lnLs[bp_lnLs[:,0].argsort()]
# plot(sorted_bp_lnLs[:,0],sorted_bp_lnLs[:,1])
# lnLs = abs(sorted_bp_lnLs[:,1])
# lnLs_mean = mean(lnLs)
# pvalues =  array([scipy.stats.chisqprob(lnL - lnLs_mean, 2) for lnL in lnLs])


