import numpy as np
import os
from bisect import bisect_left # for BilinearInterpolation
import matplotlib.pyplot as plt

matrix_Logdelta_LogT_H2       = 'matrix_modif_Logdelta_LogT_H2.dat'
matrix_Logdelta_LogT_H2_tcool = 'matrix_modif_Logdelta_LogT_tcool.dat'
path_out                      = '/scratch2/dicerbo/plot_test/'
# global variables
redshift          = 19.0
Hubble            = 0.72
h2                = Hubble * Hubble
Omega0m           = 0.3
Omega0l           = 0.7
Omega0r           = 0.0
BOLTZMANN         = 1.3806e-16  # erg/K
PROTONMASS        = 1.6726e-24  # g
GAMMA_MINUS1      = 5./3. - 1.
HYDROGEN_MASSFRAC = 0.76
mu_h              = 4./ (5. * HYDROGEN_MASSFRAC + 3.)  #  molecular weight of the hot phase
mu_c              = 4./ (3. * HYDROGEN_MASSFRAC + 1.)  # molecular weight of the cold phase
FracC             = 0.9
NPCLASS           = 500
rho_cr            = 1.9e-29 * ((1-Omega0m-Omega0l)*pow((1+redshift),2) + Omega0m*pow((1+redshift),3) + Omega0r*pow((1+redshift),4) + Omega0l) * h2


UnitMass_in_g            = 1.989e43
UnitLength_in_cm         = 3.085678e21
UnitVelocity_in_cm_per_s = 1.0e5
UnitTime_in_s            = UnitLength_in_cm / UnitVelocity_in_cm_per_s
UnitEnergy_in_cgs        = UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2)
UnitDensity_in_cgs       = UnitMass_in_g / pow(UnitLength_in_cm,3)

SE_TO_T                  = GAMMA_MINUS1 * UnitEnergy_in_cgs * mu_h * PROTONMASS / BOLTZMANN / UnitMass_in_g
T_TO_SE                  = 1.0/SE_TO_T
Th2                      = 10.0
FTHR                     = 1.e-6

# global arrays: Temperature, H2OverDensity, H2Fraction, tcool to load UM's tables
#                T in K, tcool in Gyr
T          = None          # dimension 1x50
Dens       = None          # dimension 1x50
FH2        = None          # dimension 50x50
t_cool     = None          # dimension 50x50

def main():
    # functions to call
    options = {1:'OneT',2:'OneDelta',3:'Fit_T_Delta',4:'MolecularProfileI',5:'MolecularProfileTc',6:'plot'}
    num=1
    while num in range(1,7):
        # select what you want to do...
        print '\n\t    Select                    Action'
        print   '\t    1      H2 fraction as a function of the overdensity at a fix T'
        print   '\t    2      H2 fraction as a function of the temperature at a fixed overdensity'
        print   '\t    3      H2 fraction given one temperature and one overdensity (fitted)'
        print   '\t    4      H2 molecular profile as a function of the pressure, given an SPH density (istantaneous)'
        print   '\t    5      H2 molecular profile as a function of the pressure, given an SPH density (converging Tc)'
        print   '\t    6      make the plots'
        print   '\t    whatever else             Exit'
        num = int(raw_input("\n\t Select Action :> "))
        #num = 5
        if num in options.keys():
            eval(options[num])()
            num = 10

    print '\n\t END OF GAME!!!\n'


def LoadMatrix(filename=False):
    """
    This function loads one Umberto's file,
    returns the matrix and the corresponding edges

    """

    global matrix_Logdelta_LogT_H2
    global matrix_Logdelta_LogT_H2_tcool

    if filename==False:
        raise IOError('\n\t filename MUST be provided\n')
    
    # store the path of this module
    # locate = inspect.getfile(LoadMatrix)
    # dir_file = locate.replace('H2fraction.py','')
    # filename = dir_file+filename
    if not os.path.isfile(filename):
        raise IOError('\n\t filename ',filename,' NOT found\n')

    # load file
    matrix = np.loadtxt(filename,comments='#')

    # OverDensity edges
    global Dens ; global T ; global FH2 ; global t_cool
    Dens = matrix[0,:]
    # Temperature edges
    T = matrix[1,:]

    if filename == matrix_Logdelta_LogT_H2:
        FH2 = matrix[2:,:]
    elif filename == matrix_Logdelta_LogT_H2_tcool:
        t_cool = matrix[2:,:]
    else:
        raise IOError('\n\t It seems that ',filename,' does not exist\n')



def plot():
    """
    This function makes the following plots:
    - T vs Density and as color the H2Fraction
    - T vs Density and as color the cooling time
    """

    # first plot
    global matrix_Logdelta_LogT_H2
    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    global T ; global Dens ; global FH2
    H2 = FH2
    H2[H2 > 0.] = np.log10(H2[H2 > 0.])
    vmin = -6
    vmax = -2.
    H2[H2 == 0.] = vmin
    H2[H2 > vmax] = vmax
    H2[H2 < vmin] = vmin
    nlev = 15
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin
    plt.figure()
    figura = plt.contourf(Dens,T,H2,levels,extend='both')
    figura0 = plt.contour(Dens,T,H2,levels,colors = ('k',),linewidths = (0.3,))
    plt.title('H$_{2}$ fraction')
    plt.xlabel('log10 $\delta$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('H$_{2}$ fraction',fontsize=20)
    plt.savefig('fraction_H2.pdf')
    plt.close('all')
    print '\n\t fraction_H2.jpg done\n'

    # second plot
    global matrix_Logdelta_LogT_H2_tcool
    LoadMatrix(filename=matrix_Logdelta_LogT_H2_tcool)
    global t_cool
    cool = t_cool
    cool[cool > 0.] = np.log10(cool[cool > 0.])
    vmin = -5
    vmax = 7.
    cool[cool == 0.] = vmin
    cool[cool > vmax] = vmax
    cool[cool < vmin] = vmin
    nlev = 15
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin
    plt.figure()
    figura = plt.contourf(Dens,T,cool,levels,extend='both')
    figura0 = plt.contour(Dens,T,cool,levels,colors = ('k',),linewidths = (0.3,))
    plt.title('H$_{2}$ fraction')
    plt.xlabel('log10 $\delta$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('tcool   [Gyr]',fontsize=20)
    plt.savefig('fraction_H2_tcool.pdf')
    plt.close('all')
    print '\n\t fraction_H2_tcool.jpg done\n'




main()