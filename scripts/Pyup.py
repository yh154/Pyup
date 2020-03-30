#!/usr/bin/python
from bisect import bisect_left, insort
import numpy as np
import pandas as pd
import sys
import os
import csv
import argparse
import pymc3 as pm
import arviz as az
from scipy.special import gammaln
from scipy.special import psi
from scipy.special import factorial
from scipy.optimize import fmin_l_bfgs_b as optim

seed = 123
np.random.seed(seed)

'''
def clean_zero(df):
    i_1 = list(range(0,53))
    i_2 = list(range(df.shape[0]-53, df.shape[0]))
    df.drop(i_1+i_2, inplace=True)
    df= df.loc[df['depth'] > 0]
    df.reset_index(inplace=True)
    return df
'''
def clean_zero(df):
    rr = df.shape[0]
    df = df.loc[df['coverage'] > 0]
    df.sort_values(['chr', 'pos'], ascending=[1, 1])
    df.reset_index(inplace=True, drop=True)
    print("{} zero coverage removed\n".format(rr-df.shape[0]))
    return df

def running_median(X, width, container=list):
    rm = RunningMedian(X, width, container=list)
    return rm.run()

class RunningMedian:

    def __init__(self, data, width, container=list):

        if (width % 2) != 1:
            print("Warning[sequana]:: window length should be odd. Added +1.")
            width += 1

        self.container = container
        self.W = width
        self.data = data

    def __call__(self):
        return self.run()

    def run(self):

        # initialise with first W values and sort the values
        lc = self.container(self.data[:self.W])
        lc.sort()  #

        mididx = (self.W - 1) // 2
        result = np.empty_like(self.data)

        # We initialise the first element
        idx = mididx
        result[idx] = lc[mididx]

        # We start at position W removing first element in lc
        # and adding a new one. We do not use enumerate since we do not
        # start at zero.
        for new_elem in self.data[self.W:]:
            old_elem = self.data[idx-mididx]
            del lc[bisect_left(lc, old_elem)]
            insort(lc, new_elem)
            idx += 1
            result[idx] = lc[mididx]

        # We decided to keep the first W/2 and last W/2 values as in the
        # original data. Ideally, they should not be used for post processing
        result[0:mididx] = self.data[0:mididx]
        result[-mididx:] = self.data[-mididx:]

        return result

def nb_estimate(data, mu, alpha, **kwargs):
    if len(data) > 100000:
        import random
        #random.seed(123)
        indices = random.sample(range(len(data)), 100000)
        cv = [data[i] for i in indices]
    else:
        cv=data

    md = ModelEstimate(cv, mu, alpha, **kwargs)
    if len(kwargs) > 0:
        return md.nb_mixture()
    else:
        return md.nb_model()

class ModelEstimate:

    def __init__(self, data, mu, alpha, **kwargs):
        self.mu = mu
        self.alpha = alpha
        self.data = data
        self.kwargs = kwargs

    def nb_mixture(self, N=1000, tune=1000):
        dat = np.asarray(self.data)
        #kwargs = self.kwargs
        #if len(kwargs) < 1:
        #    print("Missing args for nb mixture model estimation")
        #    sys.exit(2)
        kwargs = self.kwargs
        if len(kwargs)<1:
            print("missing args")
            sys.exit(0)
        print(np.max(dat))
        with pm.Model() as model:
            mu1 = pm.Uniform('mu1', lower=1, upper=self.mu)
            mu2 = pm.Uniform('mu2', lower=1, upper=self.kwargs['mu2'])

            alpha1 = pm.Uniform('alpha1', lower=0, upper=self.alpha)
            alpha2 = pm.Uniform('alpha2', lower=0, upper=self.kwargs['alpha2'])

            w = pm.Dirichlet('w', a=np.array([1, 1]))

            nb1 = pm.NegativeBinomial.dist(mu=mu1, alpha=alpha1)
            nb2 = pm.NegativeBinomial.dist(mu=mu2, alpha=alpha2)

            like = pm.Mixture('like', w=w, comp_dists=[nb1, nb2], observed=dat)
            trace_n = pm.sample(N, tune=tune,cores=2)
        return trace_n

    def nb_model(self, N=1000, tune=1000):
        dat = self.data
        mu = self.mu
        alpha = self.alpha

        dat = np.asarray(dat)
        dat[dat > 10] = 0
        dat = dat[dat > 0]
        print(np.max(dat))
        with pm.Model() as model_n:
            mu = pm.Uniform('mu', lower=0, upper=mu)
            alpha = pm.Uniform('alpha', lower=0, upper=alpha)
            # y_pred = pm.NegativeBinomial('y_pred', mu=mu, alpha=alpha)
            y_est = pm.NegativeBinomial('y_est', mu=mu, alpha=alpha, observed=dat)
            trace_n = pm.sample(N, tune=tune, cores=2)

        return trace_n

def fit_nbinom(X, initial_params=None):
    infinitesimal = np.finfo(np.float).eps
    X=X[X>0]

    def log_likelihood(params, *args):
        r, p = params
        X = args[0]
        N = X.size

        result = np.sum(gammaln(X + r)) - np.sum(np.log(factorial(X))) \
            - N*(gammaln(r)) + N*r*np.log(p) + np.sum(X*np.log(1-(p if p < 1 else 1-infinitesimal)))

        return -result

    if initial_params is None:
        #reasonable initial values (from fitdistr function in R)
        m = np.mean(X)
        v = np.var(X)
        size = (m**2)/(v-m) if v > m else 10

        #convert mu/size parameterization to prob/size
        p0 = size / ((size+m) if size+m != 0 else 1)
        r0 = size
        initial_params = np.array([r0, p0])

    bounds = [(infinitesimal, None), (infinitesimal, 1)]
    optimres = optim(log_likelihood,
                     x0=initial_params,
                     args=(X,),
                     approx_grad=True,
                     bounds=bounds)

    params = optimres[0]
    mu = params[0]*(1-params[1])/(params[1] if params[1]>0 else infinitesimal)
    sigmasqr = params[0]*(1-params[1])/((params[1]**2) if params[1]>0 else infinitesimal)

    return {'mu': mu, 'size': params[0], 'prob': params[1], 'var': sigmasqr}

def cal_residual(X, mu=None, var=None):
    if mu is None or var is None:
        sys.exit('calculate residual missing mean or variances!')

    residuals = (X - mu)/ (var if var > 0 else infinitesimal)

    return residuals.round(3)

def get_region(DF, site_width=26):
    df = DF.sort_values(['chr', 'pos'], ascending=[1, 1])
    roi = list()
    for index, row in df.iterrows():

        ichr, ipos, idepth, irm, irate, ires = row.iteritems()

        if index == 0:
            chr = int(ichr[1]) if isinstance(ichr[1],float) else ichr[1]
            start = end = ipos[1]
            imd_dp = [idepth[1]]
            imd_res = [ires[1]]
        else:
            if (ipos[1] - end == 1):
                end = ipos[1]
                imd_dp.append(idepth[1])
                imd_res.append(ires[1])
            else:
                roi.append([str(chr), int(start), int(end), np.median(imd_dp), int(end-start+1), np.median(imd_res)])
                chr = int(ichr[1]) if isinstance(ichr[1],float) else ichr[1]
                start = end = ipos[1]
                imd_dp = []
                imd_res = []

    roi.append([str(chr), int(start), int(end), '.', np.median(imd_res), '.', int(end - start + 1), np.median(imd_dp)])

    regions = pd.DataFrame(roi, columns =['chr', 'start', 'end', '.', 'median_zscore','.','width', 'median_coverage'])
    regions = regions.loc[regions['width']>=site_width]
    regions = regions.sort_values(['median_coverage', 'chr','start'], ascending=[0, 1, 1])
    regions.reset_index(inplace=True, drop=True)

    #return regions
    return regions

class AllNanAxisWarning(RuntimeWarning):
    pass

def identify_roi(DF, width, rate=2, rm=10):

    df = DF.loc[(DF['Rate'+str(width)] > rate) & (DF['RM'+str(width)] > rm)]
    df.reset_index(inplace=True, drop=True)
    print("Summary Regions ....")
    roi = get_region(df,site_width)
    roi.reset_index(inplace=True, drop=True)
    return {'coverage': df, 'region': roi}


def main(args):
    fout = args.output
    if os.path.isfile(fout+".csv") or os.path.isfile(fout+".roi.txt"):
        overwrite = input('\nFile already exists. Overwrite? Y/N ')
        if overwrite.lower() == 'n':
            exit(0)

    fin = args.input
    if not os.path.isfile(fin):
        sys.exit('Input file does not exit!')

    W = args.window

    df = pd.read_table(fin, header='infer', names=["chr","pos","coverage"])
    df=clean_zero(df)

    print("Calculate Running Median ...")
    RM = running_median(X=df['coverage'], width=W)
    print("Estimate background coverage ...")
    # For 10x Single Cell we consider C0 no more than 10, while may not be true for MT.
    # However, this would affect the Residual values, but not region of interest.
    cvb = fit_nbinom(X=df.loc[df['coverage'] <= 10]['coverage'])
    for key in cvb:
        cvb[key]=round(cvb.get(key), 3)
    print(cvb)

    Res = cal_residual(X = df['coverage'], mu=cvb['mu'], var=cvb['var'])
    df['RM'+str(W)] = RM
    df['Rate'+str(W)] = df['coverage']/df['RM'+str(W)]
    df['Zscore'] = Res
    df.to_csv(fout + ".details.csv", index=False)

    print("Identify Regions of Interest ...")
    roi = identify_roi(df, W)
    roi['region'].to_csv(fout+".roi.txt", sep="\t", index=False)
    roi['coverage'].to_csv(fout+".csv", index=False)

    print("Done!\n")


fout=snakemake.output[0].split('.')[1]

if os.path.isfile("chr."+fout+".csv") or os.path.isfile("chr."+fout+".roi.txt"):
        overwrite = input('\nFile already exists. Overwrite? Y/N ')
        if overwrite.lower() == 'n':
            exit(0)

fin = snakemake.input[0]
if not os.path.isfile(fin):
        sys.exit('Input file does not exit!')

W = snakemake.params[0] * 2 + 1
site_width = snakemake.params[0]
ratio = snakemake.params[1]
min_depth = snakemake.params[2]



df = pd.read_table(fin, header='infer', names=["chr","pos","coverage"])
df=clean_zero(df)

print("Calculate Running Median ...")
RM = running_median(X=df['coverage'], width=W)
print("Estimate background coverage ...")
# For 10x Single Cell we consider C0 no more than 10, while may not be true for MT.
# However, this would affect the Residual values, but not region of interest.
cvb = fit_nbinom(X=df.loc[df['coverage'] <= 10]['coverage'])
for key in cvb:
    cvb[key]=round(cvb.get(key), 3)
print(cvb)
Res = cal_residual(X = df['coverage'], mu=cvb['mu'], var=cvb['var'])

df['RM'+str(W)] = RM
df['Rate'+str(W)] = df['coverage']/df['RM'+str(W)]
df['Zscore'] = Res
df.to_csv("chr."+fout + ".details.csv", index=False)

print("Identify Regions of Interest ...")
roi = identify_roi(df, W, rate=ratio, rm=min_depth)
roi['region'].to_csv("chr."+fout+".roi.txt", sep="\t", index=False, header=False)
roi['coverage'].to_csv("chr."+fout+".csv", index=False)
print("Done!\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input', help='Input coverage, tab-delimited with three columns: chr, pos, coverage', default=None)
    parser.add_argument('--output', dest='output', help='output csv file', default="output.csv")
    parser.add_argument('--window', dest='window', help='Running median window size.', type=int, default=None)

    if len(sys.argv) == 1:
        #parser.print_help()
        sys.exit()
    args = parser.parse_args()
