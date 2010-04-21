import numpy as np
import swignifit as sf
import operator as op
import re


def extract_subclasses(base):
    to_visit = base.__subclasses__()
    subclasses = dict()
    for cl in to_visit:
        descriptor = cl.getDescriptor()
        if descriptor not in subclasses.keys():
            subclasses[descriptor] = cl
            to_visit.extend(cl.__subclasses__())
    return subclasses

sig_dict = extract_subclasses(sf.PsiSigmoid)
core_dict = extract_subclasses(sf.PsiCore)

class PsignifitException(Exception):
    pass

def get_sigmoid(descriptor):
    if not sig_dict.has_key(descriptor):
        raise PsignifitException("The sigmoid \'"+str(descriptor)+"\' you requested, is not available.")
    return sig_dict[descriptor]()

def get_core(descriptor, data, sigmoid_type):
    descriptor, parameter = re.match('([a-z]+)([\d\.]*)', descriptor).groups()
    if descriptor not in core_dict.keys():
        raise PsignifitException("The core \'"\
                +str(descriptor)\
                +"\' you requested, is not available.")
    if len(parameter) > 0:
        return core_dict[descriptor](data, sigmoid_type, float(parameter))
    else:
        return core_dict[descriptor](data, sigmoid_type)

def get_prior(prior):
    try:
        prior = "sf."+"Prior(".join(prior.split('('))
        return eval(prior)
    except Exception, e:
        return None

def get_cuts(cuts):
    if cuts is None:
        return sf.vector_double([0.5])
    elif op.isNumberType(cuts):
        return sf.vector_double([cuts])
    elif op.isSequenceType(cuts) and np.array([op.isNumberType(a) for a in cuts]).all():
        return sf.vector_double(cuts)
    else:
        raise PsignifitException("'cuts' must be either None, a number or a "+\
                "sequence of numbers.")

def available_sigmoids():
    print "The following sigmoids are available:"
    print sig_dict.keys()


def bootstrap(data, start=None, nsamples=2000, nafc=2, sigmoid="logistic",
        core="ab", priors=None, cuts=None, parametric=True ):

    data = np.array(data).T
    x = sf.vector_double(data[0])
    k = sf.vector_int(data[1].astype(int))
    N = sf.vector_int(data[2].astype(int))
    data = sf.PsiData(x,N,k,nafc)
    sigmoid = get_sigmoid(sigmoid)
    core = get_core(core, data, sigmoid.getcode())
    pmf = sf.PsiPsychometric(nafc, core, sigmoid)
    nparams = pmf.getNparams()
    if priors is not None:
        if len(priors) != nparams:
            raise PsignifitException("You specified \'"+str(len(priors))+\
                    "\' priors, but there are \'"+str(nparams)+ "\' parameters.")
        for (i,p) in enumerate((get_prior(p) for p in priors)):
            if p is not None:
                pmf.setPrior(i, p)

    cuts = get_cuts(cuts)
    ncuts = len(cuts)
    if start is not None:
        if len(start) != nparams:
            raise PsignifitException("You specified \'"+str(len(start))+\
                    "\' starting value(s), but there are \'"+str(nparams)+ "\' parameters.")
        start = sf.vector_double(start)

    bs_list = sf.bootstrap(nsamples, data, pmf, cuts, start, True, parametric)
    jk_list = sf.jackknifedata(data, pmf)

    nblocks = data.getNblocks()

    # construct the massive tuple of return values
    samples = np.zeros((nsamples, nblocks), dtype=np.int32)
    estimates = np.zeros((nsamples, nparams))
    deviance = np.zeros((nsamples))
    thres = np.zeros((nsamples, ncuts))
    Rpd = np.zeros((nsamples))
    Rkd = np.zeros((nsamples))
    for row_index in xrange(nsamples):
        samples[row_index] = bs_list.getData(row_index)
        estimates[row_index] = bs_list.getEst(row_index)
        deviance[row_index] = bs_list.getdeviance(row_index)
        thres[row_index] = [bs_list.getThres_byPos(row_index, j) for j in xrange(ncuts)]
        Rpd[row_index] = bs_list.getRpd(row_index)
        Rkd[row_index] = bs_list.getRkd(row_index)

    acc = np.zeros((ncuts))
    bias = np.zeros((ncuts))
    for cut in xrange(ncuts):
        acc[cut] = bs_list.getAcc(cut)
        bias[cut] = bs_list.getBias(cut)

    ci_lower = sf.vector_double(nparams)
    ci_upper = sf.vector_double(nparams)

    for param in xrange(nparams):
        ci_lower[param] = bs_list.getPercentile(0.025, param)
        ci_upper[param] = bs_list.getPercentile(0.975, param)

    outliers = np.zeros((nblocks), dtype=np.bool)
    influential = np.zeros((nblocks))

    for block in xrange(nblocks):
        outliers[block] = jk_list.outlier(block)
        influential[block] = jk_list.influential(block, ci_lower, ci_upper)

    return samples, estimates, deviance, thres, bias, acc, Rpd, Rkd, outliers, influential

