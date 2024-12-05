import numpy
# import pylab
import oskar_theory
import corrections
from scipy.optimize import newton, brentq
from math import exp, log
import coarse_grained_theory
import sfs_utils
import coalescent_sim

c = 4.0
p = 1.0


def condition(x, Ns, NUd, NR):
    variance = oskar_theory.calculate_var_k(Ns, NUd*x)*Ns*Ns
    return (1+(corrections.minimum_t2(variance)*NR/c)**p)*(x**p)-1


def calculate_asexual_fraction(Ns, NUd, NR):
    return brentq(lambda x: condition(x, Ns, NUd, NR), 0.00001, 0.99999)


def calculate_effective_params(Ns, NUd, NR, x):

    # Note: these effective parameters will only be accurate
    # in the interference selection regime when the linkage block
    # approximation is valid.

    Ns = Ns*1.0
    NUd = NUd*1.0
    NR = NR*1.0
    x = x*1.0

    Nseff, NUdeff = coarse_grained_theory.calculate_critical_params(
        Ns*Ns*oskar_theory.calculate_var_k(Ns, NUd*x))
    
    if ((Ns*exp(-2*NUd / (2*Ns + NR))) < 5): 
        NReff = 0
        return Nseff, NUdeff, NReff
    else:
        return(Ns, NUd, NR)


def inverse_NUeff(Ns, NUeff, NR):
    return newton(lambda x: calculate_asexual_fraction(Ns, x, NR)*x-NUeff, NUeff*NR)


def calculate_t2(Ns, NUd, NR, x, use_coalescent_sim=False):
    # Calculates mean pairwise coalescent time, T2/N
    # also equal to pi/pi0 from the text
    Nseff, NUdeff, NReff = calculate_effective_params(Ns, NUd, NR, x)
    return corrections.calculate_t2(Nseff, NUdeff, NReff)


def calculate_Q(Ns, NUd, NR, n, num_samples=100000):
    # calculates the average site frequency spectrum
    # normalized by the number of singletons
    # This is Q_n(i) from the text
    Nseff, NUdeff = calculate_effective_params(Ns, NUd, NR)
    Q = coalescent_sim.calculate_time_spectrum(
        Nseff, NUdeff, n, num_samples=num_samples)
    return Q/Q[0]


def calculate_maf(Ns, NUd, NR, n):
    # Calculates the average minor allele frequency
    return sfs_utils.calculate_maf(calculate_Q(Ns, NUd, NR, n, num_samples=10000))


def calculate_schaeffers_D(Ns, NUd, NR, n):
    # Calculates the Schaeffers D,
    # defined as Tajimas D/Dmin
    return sfs_utils.calculate_schaeffers_D(calculate_Q(Ns, NUd, NR, n, num_samples=10000))


if __name__ == '__main__':
    Ns = 10
    NU = 300
    NR = 100
    print("Calculating effective parameters for Ns=%g, NU=%g, NR=%g" % (Ns, NU, NR))
    maf = calculate_maf(Ns, NU, NR, 10)
    print("Success!")
