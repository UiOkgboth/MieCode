import collections as col
from  scipy.interpolate import interp1d
import numpy as np
import scipy.special as sps
import warnings 
#Class for particle elements
class Elements:
    elementCount = 0
    
    def __init__(self, chemSymb, wp, gamma, vf, A, RefracFile):
        self.chemSymb = chemSymb
        self.wp = wp #E15 Hz
        self.gamma = gamma #s
        self.vf = vf #m/s
        self.A = A
        self.RefracFile = RefracFile
        Elements.elementCount += 1

#Class for matrix compounds		
class Compounds:
    compoundsCount = 0
    
    def __init__(self,chemFormula,RefracFile):
        self.chemFormula = chemFormula
        self.RefracFile = RefracFile
        Compounds.compoundsCount += 1
		
#named Tuple, contains information about the particle geometry, and physical constants
def new_particle(name, CoreParticle,Shape, ClusterS ,R ,wp ,gamma ,vf ,A):
    particle = col.namedtuple(name,['CoreParticle', 'Shape', 'ClusterS', 'R', 'wp', 'gamma', 'vf', 'A'])
    particle.CoreParticle = CoreParticle
    particle.Shape = Shape
    particle.ClusterS = ClusterS
    particle.R = R #E-9 m
    particle.wp = wp #E15 Hz
    particle.gamma = gamma #E15 s
    particle.vf = vf #E6 m/s
    particle.A = A
    return particle

#new calculation, passing on the important variables used dufing the calculation
def new_calculation(name, ParticleMaterial, MatrixMaterial, radii =10, N_of_calc =1, LowerL = 124, HigherL = 1240):
    Ncal = col.namedtuple(name,['ParticleMaterial', 'MatrixMaterial', 'radii', 'N_of_calc', 'LowerL', 'HigherL'])
    Ncal.ParticleMaterial = ParticleMaterial
    Ncal.MatrixMaterial = MatrixMaterial
    Ncal.radii = radii
    Ncal.N_of_calc = N_of_calc
    Ncal.LowerL = LowerL
    Ncal.HigherL = HigherL
    return Ncal
	
#Returns the entries of a file with one line header and all other lines entries, separated with (,) if multiple entries per line
def loadfile(FileName):
    data = np.loadtxt(FileName,skiprows=1,delimiter=",")
          
    data[:,0] = data[:,0]*1e3 #from um to nm (um default used on refractiveindex.info)
    return data

	#Extrapolates the initial data to fit within the desired range 
def extrapolate(Idata,LowerL = 124, UpperL = 1240):
    Isize = (3*Idata[:,0].size == Idata.size)
    
    data = np.zeros([Idata[:,0].size,3])
    if Isize == False:
        data[:,0] = Idata[:,0]
        data[:,1] = Idata[:,1]
    else:
        data = Idata
    
    
    step = 1
    points = 4    
    
    
    if data[0,0] > LowerL:
        minrange = np.zeros([points,3])
        for i in range (0,points , 1):
            minrange[i,:] = data[0,:]
            data = np.delete(data,0,0)
        x_val = np.arange(LowerL, minrange[points-1,0], 1)
        
        nmin_f = np.polyfit(minrange[:,0],minrange[:,1],2)
        nminp = np.poly1d(nmin_f)
        nmin = nminp(x_val)
        if Isize:
            kmin_f = np.polyfit(minrange[:,0],minrange[:,2],2)
            kminp = np.poly1d(kmin_f)
            kmin = kminp(x_val)
        else:
            kmin = np.zeros(nmin.size)
        min_arr = np.transpose(np.stack((x_val,nmin,kmin),axis=0))
        extrapolate_temp = np.concatenate((min_arr, data))
    if data[-1,0] < UpperL:
        maxtange = np.zeros([points,3])
        for i in range (points,0, 1):
            maxrange[i-1,:] = data[-1,:]
            data = np.delete(data,-1,0)
        max_val = np.arange(UpperL, maxrange[0,0],1)
        
        nmax_f = np.polyfit(maxrange[:,0],maxrange[:,1],2)
        nmaxp = np.poly1d(nmax_f)
        nmax = nmaxp(max_val)
        if isize:
            kmax_f = np.polyfit(maxrange[:,0], maxrange[:,2],2)
            kmaxp = np.poly1d(kmax_f)
            kmax = kmaxp(max_val)
        else:
            kmin = np.zeros(nmin.size)
        

        max_arr = np.transpose(np.stack((max_val, nmax, kmax),axis=0))
        extrapolate_data = np.concatenate((data, max_arr))
    else:
        extrapolate_data = extrapolate_temp
        
    return extrapolate_data
	
#interpolates a data array [:,3] to have a value for every nm, by default 
def interpolate(data, lowerL = 200, higherL = 1240, step = 1):
        
    x_val = np.arange(lowerL, higherL, step, dtype = int)
    
    nval = interp1d(data[:,0], data[:,1], kind = 'cubic')
    ninp = nval(x_val)
    
    if 3*data[:,0].size == data.size:
        kval = interp1d(data[:,0], data[:,2], kind = 'cubic')
        kinp = kval(x_val)
    else:
        kinp = np.zeros(ninp.size)
    
    interp_data = np.transpose(np.stack((x_val,ninp,kinp),axis = 0))

    return interp_data

#calculates the real and imaginary part of a dielectric function based on n and k, where n is data[:,0] and k data[:,1]
def calc_diel_f(data0, data1):
    eps = np.zeros([data0.size,2])
    eps[:,0] = data0**2 - data1**2
    eps[:,1] = 2*data0*data1
    return eps

#Applies changes to the dielectric function of the particle based on the 
def ClusterSizeEffects(Eps, wl, R, wp, gamma, vf, A):
    
    freq = 1883.651565/wl
    freq2 = freq**2
    
    dgamma = gamma + A*vf/R
    dgamma2 = dgamma
    
    Eps[:,0] = Eps[:,0] + wp**2/(freq2+gamma**2)-wp**2/(freq2+dgamma2)
    Eps[:,1] = Eps[:,1] + wp**2/freq*(dgamma/(freq2+dgamma2)-gamma/(freq2+gamma**2))
    
    return Eps
	
def bessel_jn(l,x,der = False):
    jn = sps.spherical_jn(l, x, derivative = der)
    return jn

def bessel_yn(l,x,der = False):
    yn = sps.spherical_yn(l, x, derivative = der)
    return yn

def CoreSphereParticle(wl, epsP, epsM, R):
    X0 = 2*np.pi*R
    epsM = np.sqrt(epsM)
    
    epsP0 = np.sqrt(epsP[:,0]**2+epsP[:,1]**2)/2
    Y0 = np.sqrt(epsP[:,0]/2+epsP0)+1j*np.sqrt(-epsP[:,0]/2+epsP0)
    
    X = X0/wl*epsM[:,0]
    Y = X0/wl*Y0
    M = Y/X
    
    A = np.zeros(X.size)
    B = np.zeros(X.size)
    Qext = np.zeros(X.size)
    Qsca = np.zeros(X.size)
    
    
    L = 1
    while True:
        bjY  = bessel_jn(L,Y)
        bjX  = bessel_jn(L,X)
        byX  = bessel_yn(L,X)
        bjXp = bessel_jn(L,X,True)
         
        pmx  = Y*bjY
        pmxp = bjY + Y*bessel_jn(L,Y,True)
            
        px   = X*bjX
        pxp  = bjX + X*bjXp
        
        ex   = X*bjX +1j*X*byX
        exp  = bjX + X*bjXp + 1j*byX + 1j*X*bessel_yn(L,X,True)
                
        A = (M*pmx*pxp-pmxp*px)/(M*pmx*exp-pmxp*ex)
        B = (pmx*pxp-M*pmxp*px)/(pmx*exp-M*pmxp*ex)
        warnings.filterwarnings("ignore")
        if np.isnan(np.sum(A)) or np.isnan(np.sum(B)):
            print("NaN value at L = ", L)

            break
            
        AB = (A+B).real*(2*L+1)
        ABabs = (A*np.conj(A))+(B*np.conj(B))*(2*L+1)

        Qext = Qext + AB
        Qsca = Qsca + ABabs
        L += 1
        
    Qext = 2*Qext/X**2
    Qsca = 2*Qsca/X**2 
   
    return Qext, Qsca
