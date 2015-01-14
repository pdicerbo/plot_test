import numpy as np
import os
import string
from bisect import bisect_left # for BilinearInterpolation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata

matrix_Logdelta_LogT_H2       = 'matrix_modif_Logdelta_LogT_H2.dat'
matrix_Logdelta_LogT_H2_tcool = 'matrix_modif_Logdelta_LogT_tcool.dat'
path_in                       = '/scratch2/dicerbo/plot_test/toprint/'
path_read                     = '/scratch2/dicerbo/plot_test/firstmat/'
path_out                      = '/scratch2/dicerbo/plot_test/boundary/'
path_matrix                   = '/scratch2/dicerbo/plot_test/matrix/'
path_file                     = '/scratch2/dicerbo/plot_test/toprint/t700/time_evolution_log10P4.7.dat'
path_file2                    = '/scratch2/dicerbo/plot_test/toprint/t700/time_evolution_log10P3.5.dat'

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
    global path_in
    '''
    directory = '/home/dicerbo/output/scratch2/plot_test/prova/prova1'
    if not os.path.exists(directory):
        os.makedirs(directory)
        print '%s created successfully'%(directory)
    '''
    dirs = os.listdir(path_in)
    for d in dirs:
        print '\n\tStart working on '+ d
        adjust(path_in+d+'/')
        plot_def(d)
        print '\n\tEnd working on ' + d

    print '\n\tFinally End\n'

def adjust(path):

    print '\n\tWithin adjust function\n'

    global matrix_Logdelta_LogT_H2
    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    global T ; global Dens ; global FH2
    global path_out

    tmin = T.min(); tmax = T.max()
    dmin = Dens.min(); dmax = Dens.max()
    #print str(tmin)+'\t'+str(tmax)
    #print str(dmin)+'\t'+str(dmax)
    files = os.listdir(path)
    directory = path[36:-1]
    for name in files:
        if string.count(name, 'time') != 0:
            print '\n\tWorking on '+name
            matrix = np.loadtxt(path+name, comments = '#')
            if matrix.size < 4:
                print '\n\t'+name+' empty\n'
            else:
                matrix = matrix.T
                time = matrix[0,:]
                rho_atom = matrix[1,:]
                t_atom = matrix[2,:]
                control = 0
                #newmat = np.array([[0.]*time.size for x in range(3)])
                newmat = np.zeros((3, time.size),dtype=float)
                lst = range(time.size)
                for j in lst:
                    if rho_atom[j] > dmin and rho_atom[j] < dmax and t_atom[j] >= tmin and t_atom[j] <= tmax:
                        newmat[0][j] += time[j]
                        newmat[1][j] += rho_atom[j]
                        newmat[2][j] += t_atom[j]
                        control += 1
                if control > 2.:
                    #toprint = np.zeros((3, control),dtype=float)
                    toprint = newmat.T[:control,:]
                    wr = open(path_out+directory+'_boundary_ctrl_T'+name[-7:], 'w')
                    np.savetxt(wr, toprint, fmt = '%e', delimiter = '\t', newline = '\n')
                    wr.close()

        else:
            print "\n\tFile " + name + " is for Blitz&Rosolowsky's plot -> Continue\n"
            continue


def plot_def(directory):
    print '\n\tWithin plot function\n'
    global path_in; global path_read; global path_out
    # first plot
    global matrix_Logdelta_LogT_H2
    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    global T ; global Dens ; global FH2

    H2 = FH2
    H2[H2 > 0.] = np.log10(H2[H2 > 0.])
    v_min = -6
    v_max = -2.
    H2[H2 == 0.] = v_min
    H2[H2 > v_max] = v_max
    H2[H2 < v_min] = v_min
    numlev = 15
    dmag0 = (v_max - v_min) / float(numlev)
    levels0 = np.arange(numlev) * dmag0 + v_min


    plt.figure()    
    plt.title('Paths')
    figura = plt.contourf(Dens,T,H2,levels0,extend='both')
    cbar = plt.colorbar(figura,format='%3.1f', shrink=0.7)
    cbar.set_ticks(np.linspace(v_min,v_max,num=levels0.size,endpoint=True))
    cbar.set_label('H$_{2}$ fraction',fontsize=20)
    
    #path's plot
    files = os.listdir(path_out)
    for name in files:
        if string.count(name, 'boundary') != 0 and string.count(name, directory) != 0:
            print '\n\tPlotting ' + name + ' file'
            plt.plotfile(path_out+name, delimiter = '\t', cols=(1, 2), comments='#',  marker='.', newfig=False)

    plt.xlabel('log10 $Rho$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    newname = path_out + 'path_' + directory + '.jpg'
    plt.savefig(newname)
    plt.close('all')
    print '\n\t'+newname[len(path_out):]+' done\n'

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


main()
