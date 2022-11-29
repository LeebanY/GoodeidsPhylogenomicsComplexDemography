seed =  -1

seqfile = bpp_aligned_loci_Goodeids.phy
Imapfile = Imap.txt
outfile = output_Run2.txt
mcmcfile = mcmc_Run2.txt

# fixed number of species/populations 
speciesdelimitation = 0

# fixed species tree
speciestree = 0

species&tree = 9  CB GM AT CL GA XC AS IF XR
                  1  1  1  1  1  1  1  1  1
                 (CB,((GM,AT),(CL,((GA,(XC,AS)),(IF,XR)))));
# unphased data for all 4 populations
phase =   1  1  1  1  1  1  1  1  1

# use sequence likelihood
usedata = 1

nloci = 2740

# do not remove sites with ambiguity data
cleandata = 0

# MCMC samples, locusrate, heredityscalars, Genetrees
thetaprior = 3 0.003 e  # Inv-gamma(a, b) for theta (integrated out by default; add E to also sample theta)
tauprior = 3 0.03     # Inv-gamma(a, b) for root tau
phiprior = 1 1  # Beta(a, b) for root tau & Dirichlet(a) for other tau'
heredity = 1 4 4
finetune =  1: 3 0.003 0.002 0.00002 0.005 0.9 0.001 0.001 # finetune for GBtj, GBspr, theta, tau, mix

print = 1 0 0 0   * 
burnin = 8000
sampfreq = 2
nsample = 100000
