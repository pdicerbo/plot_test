import numpy as np
import os
from bisect import bisect_left # for BilinearInterpolation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata

matrix_Logdelta_LogT_H2       = 'matrix_modif_Logdelta_LogT_H2.dat'
matrix_Logdelta_LogT_H2_tcool = 'matrix_modif_Logdelta_LogT_tcool.dat'
path_out                      = '/scratch2/dicerbo/plot_test/'
path_read                     = '/home/dicerbo/output/scratch2/plot_test/toprint/'
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
PS_path    = [[0.]*50 for x in range(50)]

def main():
    global PS_path
    '''
    lst = np.arange(500, 1500, 200)
    lst2 = ['a', 'b', 'c', 'd', 'e']
    i = 0
    for n in lst:
        lst2[i] = path_read+'t'+str(n)
        i += 1
    '''
    dirs = os.listdir(path_read)
    i = 0
    for a in dirs:
        dirs[i] = path_read+'/'+a
        i += 1

    j = 0
    for a in dirs:
        files = os.listdir(a)
        for b in files:
            get_matrix(a+'/'+b)

    fp = open('total_phase_space_path.dat','w')
    np.savetxt(fp, PS_path, fmt='%e', delimiter='\t', newline='\n')
    fp.flush()
    fp.close()

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


def get_matrix(fname):
    global Dens ; global T ; global FH2 ; global t_cool; global PS_path
    if T==None or Dens==None or FH2==None:
    # load the matrix
	    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    if t_cool==None:
        LoadMatrix(filename=matrix_Logdelta_LogT_H2_tcool)

    #fname = 'time_evolution_log10P4.494.dat'
    #fname = 'time_evolution_log10P3.996.dat'
    matrix = np.loadtxt(fname,comments='#') #contiene: t - M_a - M_h2 - Log10Rho_a - Log10T - frac_h2 - tcool
    matrix = matrix.T

    time = matrix[0,:]
    rho_atom = matrix[3,:]
    t_atom = matrix[4,:]
    #tmat = FH2

    d = Dens.size
    b = time.size
    count = 0.
    time_def = 0.
    lst_dns = np.arange(0, d, 1)
    lst_tmp = np.arange(0, d, 1)
    big_list = np.arange(0, b, 1)
    gap = (Dens[1] - Dens[0])/2.
    for i in lst_dns:
        for j in lst_tmp:
            for k in big_list:
                if np.abs(Dens[i] - rho_atom[k]) <= gap and np.abs(T[j] - t_atom[k]) <= gap:
                    time_def += time[k]
                    count += 1.
            if count < 5e-1:
                    count = 1.
            PS_path[j][i] = time_def/count
            time_def = 0.
            count = 0.
    '''
    fp = open('phase_space_path.dat','w')
    np.savetxt(fp, tmat, fmt='%e', delimiter='\t', newline='\n')
    fp.flush()
    fp.close()
    #plot_newmatrix('phase_space_path.dat')
    plot_def('phase_space_path.dat')
    print '\n\tFinally done!\n'
    '''

def plot_interp():
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

    fname = 'time_evolution_log10P4.494.dat'
    matrix = np.loadtxt(fname,comments='#') #contiene: t - M_a - M_h2 - Log10Rho_a - Log10T - frac_h2 - tcool
    matrix = matrix.T

    #tmp1 = np.linspace(-27.2, -21, num=300, endpoint=True)
    #tmp2 = np.linspace(0, 5, num=300, endpoint=True)
    #tmp3 = np.linspace(0, 1000, num=300, endpoint=True)

    dmax = Dens.max()
    dmin = Dens.min()
    tmax = T.max()
    tmin = T.min()

    time = matrix[0,:]
    time[time > 0.] = np.log10(time[time > 0.])

    rho_atom = matrix[3,:]
    rho_atom[rho_atom > dmax] = dmax
    rho_atom[rho_atom < dmin] = dmin

    t_atom = matrix[4,:]
    t_atom[t_atom > tmax] = tmax
    t_atom[t_atom < tmin] = tmin

    #tmat = griddata(Dens, T, time, rho_atom, t_atom)#, interp='linear')
    tmat = griddata((rho_atom, t_atom), time, (Dens, T), method='cubic')
    #tmat = griddata(tmp1, tmp2, tmp3, Dens, T, interp='linear')

    wr = open('tmat.dat', 'w')
    np.savetxt(wr, tmat, fmt='%e',delimiter='\t',newline='\n')
    wr.flush(); wr.close()
    vmin = 0.
    vmax = 7
    tmat[tmat == 0.] = vmin
    tmat[tmat > vmax] = vmax; tmat[tmat < vmin] = vmin

    nlev = 6
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin

    plt.figure()

    plt.contour(Dens, T, tmat, levels, linewidths=0.5, colors='k')

    plt.title('Path')
    plt.xlabel('log10 $Rho$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('time',fontsize=20)

    figura = plt.contourf(Dens,T,H2,levels0,extend='both')
    cbar = plt.colorbar(figura,format='%3.1f', orientation='horizontal', shrink=0.7)
    #cbar = plt.colorbar(figura,format='%3.1f', shrink=0.7)
    cbar.set_ticks(np.linspace(v_min,v_max,num=levels0.size,endpoint=True))
    cbar.set_label('H$_{2}$ fraction',fontsize=20)

    plt.savefig('upgrade_path.pdf')
    plt.close('all')



def plot_def(fname):
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
    v_min = -6
    v_max = -2.
    H2[H2 == 0.] = v_min
    H2[H2 > v_max] = v_max
    H2[H2 < v_min] = v_min
    numlev = 15
    dmag0 = (v_max - v_min) / float(numlev)
    levels0 = np.arange(numlev) * dmag0 + v_min

    path = np.loadtxt(fname,comments='#')
    path[path > 0.] = np.log10(path[path > 0.])
    vmin = 0.
    vmax = 6.5
    path[path == 0.] = vmin
    path[path > vmax] = vmax
    path[path < vmin] = vmin
    nlev = 6
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin

    path2 = np.loadtxt('phase_space_path.dat.old',comments='#')
    path2[path2 > 0.] = np.log10(path2[path2 > 0.])
    path2[path2 == 0.] = vmin
    path2[path2 > vmax] = vmax
    path2[path2 < vmin] = vmin

    path3 = np.loadtxt('phase_space_pathP4.992.dat',comments='#')
    path3[path3 > 0.] = np.log10(path3[path3 > 0.])
    path3[path3 == 0.] = vmin
    path3[path3 > vmax] = vmax
    path3[path3 < vmin] = vmin

    plt.figure()

    figura = plt.contour(Dens,T,path,levels,extend='both', linewidths=0.3, cmap=cm.gray)
    #figura0 = plt.contour(Dens,T,path,levels,colors = ('k',),linewidths = (0.3,))
    plt.title('Path')

    plt.xlabel('log10 $Rho$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)

    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('time',fontsize=20)
    # figura = plt.contour(Dens,T,path2,levels,extend='both', alpha=0.8)
    figura = plt.contour(Dens,T,path2,levels,extend='both', linewidths=0.3, cmap=cm.gray)
    figura = plt.contour(Dens,T,path3,levels,extend='both', linewidths=0.3, cmap=cm.gray)

    # figura = plt.contourf(Dens,T,H2,levels0,extend='both', linewidths=0.8, cmap=cm.gray, alpha=0.5)
    figura = plt.contourf(Dens,T,H2,levels0,extend='both')
    cbar = plt.colorbar(figura,format='%3.1f', orientation='horizontal', shrink=0.7)
    #cbar = plt.colorbar(figura,format='%3.1f', shrink=0.7)
    cbar.set_ticks(np.linspace(v_min,v_max,num=levels0.size,endpoint=True))
    cbar.set_label('H$_{2}$ fraction',fontsize=20)

    plt.savefig('composite_path.pdf')
    plt.close('all')

    print '\n\t composite_path.pdf done\n'



def plot_newmatrix(fname):
    """
    This function makes the following plots:
    - T vs Density and as color the H2Fraction
    - T vs Density and as color the cooling time
    """

    # first plot
    global matrix_Logdelta_LogT_H2
    LoadMatrix(filename=matrix_Logdelta_LogT_H2)
    global T ; global Dens ; global FH2
    path = np.loadtxt(fname,comments='#')
    #H2 = FH2
    path[path > 0.] = np.log10(path[path > 0.])
    vmin = 0.
    vmax = 6.5
    path[path == 0.] = vmin
    path[path > vmax] = vmax
    path[path < vmin] = vmin
    nlev = 15
    dmag = (vmax - vmin) / float(nlev)
    levels = np.arange(nlev) * dmag + vmin
    plt.figure()
    figura = plt.contourf(Dens,T,path,levels,extend='both')
    figura0 = plt.contour(Dens,T,path,levels,colors = ('k',),linewidths = (0.3,))
    plt.title('Path')
    plt.xlabel('log10 $\delta$',fontsize=20) ; plt.ylabel('Log10 T[k]',fontsize=20)
    cbar = plt.colorbar(figura,format='%3.1f')
    cbar.set_ticks(np.linspace(vmin,vmax,num=levels.size,endpoint=True))
    cbar.set_label('time',fontsize=20)
    plt.savefig('path.pdf')
    plt.close('all')
    print '\n\t path.pdf done\n'


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
    print '\n\t fraction_H2.pdf done\n'

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
    print '\n\t fraction_H2_tcool.pdf done\n'




main()