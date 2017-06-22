import numpy as np
import pandas as pd
import statistics as stat

PROTON_MASS = 0.938272


def momentum_to_energy(p):
    """Return E [GeV/c^2] from P [GeV/c] (proton)."""
    return 1000*(np.sqrt(p**2+PROTON_MASS**2)-PROTON_MASS)


def momentum_to_brho(p):
    """Return BRHO [T.m] from P [GeV/c] (proton)."""
    return 3.33564 * p


def energy_to_brho(e):
    """Return BRHO [T.m] from E [GeV] (proton)."""
    return 3.33564 * energy_to_momentum(e)


def energy_to_momentum(ekin):
    """Return P [GeV] from E [GeV/c^2] (proton)."""
    E = PROTON_MASS + ekin/1000
    return np.sqrt(E**2-PROTON_MASS**2)


def energy_to_beta(ekin):
    """Return beta relativistic from E [GeV/c^2] (proton)."""
    gamma = (PROTON_MASS + ekin) / PROTON_MASS
    return np.sqrt((gamma ** 2 - 1) / gamma ** 2)


def compute_emittance(x, px):
    """Return the beam emittance [mm.mrad] from x [mm] and px [mrad]."""
    return x * px


def compute_beta(x, px):
    """Return the Twiss beta [m] from beam shape x [mm] and px [mrad]."""
    return x / px


def range_to_energy(r):
    """Return the range [g/cm^2] from the kinetic energy [MeV]."""
    a = 0.00169; b = -0.00490; c = 0.56137; d = 3.46405
    return np.exp(
        a * np.log(r)**3 + b * np.log(r)**2 + c * np.log(r) + d
    )


def energy_to_range(e):
    """Return the kinetic energy [MeV] from the range [g/cm^2]."""
    b = 0.008539; c = 0.5271; d = 3.4917
    return np.exp((-c + np.sqrt(c**2 - 4 * b * (d - np.log(e))))/(2*b))


def compute_ess_transmission(beam_sigma, slits, dispersion):
    """Compute the transmission as a function of the momentum offset (in %) from a simple analytical model."""
    n_steps = 10000
    dx = 3.0/n_steps
    sigma = beam_sigma/2.8
    slits_at = slits/dispersion
    error = np.arange(-1.5, 1.5, dx)
    slits = np.zeros(n_steps)
    slits[np.where((error < slits_at) & (error > -slits_at))] = 1.0
    beam = np.exp(-(error/sigma)**2/2)
    return np.roll(np.convolve(slits, beam, mode="same"), -1)/np.trapz(beam)
	
	
def compute_energy_divergence(ProfileTable):
    """ Compute the energy and divergence of dataframe obtained with G4BeamLine. """
    protonMass=PROTON_MASS*1000
    Px2 = ProfileTable.Px**2
    Py2 = ProfileTable.Py**2
    Pz2 = ProfileTable.Pz**2
    ProfileTable['Ptot'] = np.sqrt(Px2+Py2+Pz2)
    ProfileTable['Energy'] = np.sqrt((protonMass*protonMass)+Px2+Py2+Pz2)-protonMass
    ProfileTable['xp'] = 1000*ProfileTable['Px']/ProfileTable['Ptot']
    ProfileTable['yp'] = 1000*ProfileTable['Py']/ProfileTable['Ptot']
    
    # Change here and take the mean of Ptot in place of energy_to_momentum()
    Pmean=ProfileTable['Ptot'].mean()
    ProfileTable['dP_P'] = (ProfileTable['Ptot']-Pmean)/(Pmean)	
    return ProfileTable
	
def compute_meanAndsigma(Data):
    """ Compute useful parameters of the beam : mean, sigma, .... """
    # For the std, use N-1 because it's an independent sample from a distributed population 
	
    xmean=Data['x'].mean()
    sigmax=Data['x'].std(ddof=1)
    
    ymean=Data['x'].mean()
    sigmay=Data['y'].std(ddof=1)
    
    xpmean=Data['xp'].mean()
    sigmaxp=Data['xp'].std(ddof=1)
    
    ypmean=Data['xp'].mean()
    sigmayp=Data['yp'].std(ddof=1)
    
    Pmean=Data['Ptot'].mean()
    sigmaP=Data['Ptot'].std(ddof=1)
    
    Emean=Data['Energy'].mean()
    sigmaE=Data['Energy'].std(ddof=1)
    
    dp_pmean=Data['dP_P'].mean()
    dp_psigma=Data['dP_P'].std(ddof=1)
    
    columnsName=['xmean','ymean','xpmean','ypmean','sigmax','sigmay','sigmaxp','sigmayp','Pmean','sigmaP','Emean','sigmaE','dp_pmean','dp_psigma']
    DataBeam=[xmean,ymean,xpmean,ypmean,sigmax,sigmay,sigmaxp,sigmayp,Pmean,sigmaP,Emean,sigmaE,dp_pmean,dp_psigma]
    DataBeam=np.array(DataBeam).reshape(1,len(DataBeam))
    Beamparameter=pd.DataFrame(DataBeam,columns=columnsName)
    
    return Beamparameter


def compute_twiss_parameter(Data):
    """ Compute TWISS parameters of the beam : alpha, beta, emittance, enveloppe, .... """ 
    
    # Data for emittance calculation 
    xmean=Data['x'].mean() 
    xpmean=Data['xp'].mean()
    ymean=Data['x'].mean()
    ypmean=Data['xp'].mean()
    # Diagonal elements
    sigma_xx = ((Data['x']-xmean)*(Data['x']-xmean)).mean()
    sigma_yy = ((Data['y']-ymean)*(Data['y']-ymean)).mean()
    sigma_xpxp = ((Data['xp']-xpmean)*(Data['xp']-xpmean)).mean()
    sigma_ypyp = ((Data['yp']-ypmean)*(Data['yp']-ypmean)).mean()
    # Off-diagonal elements
    sigma_xxp=((Data['x']-xmean)*(Data['xp']-xpmean)).mean()
    sigma_yyp=((Data['y']-ymean)*(Data['yp']-ypmean)).mean()
    # Compute the emittance
    EmitHOR=np.sqrt(sigma_xx*sigma_xpxp-sigma_xxp**2)
    EmitVER=np.sqrt(sigma_yy*sigma_ypyp-sigma_yyp**2)
    # Twiss parameter
    BetaHOR=sigma_xx/EmitHOR
    BetaVER=sigma_yy/EmitVER
    AlphaHOR=-1*sigma_xxp/EmitHOR
    AlphaVER=-1*sigma_yyp/EmitVER
    #GammaHOR=sigma_xpxp/EmitHOR 
    #GammaVER=sigma_ypyp/EmitVER 
    
    columnsName=['EmitHOR','EmitVER','BetaHOR','BetaVER','AlphaHOR','AlphaVER'] 
    DataBeam=[EmitHOR,EmitVER,BetaHOR,BetaVER,AlphaHOR,AlphaVER] 
    DataBeam=np.array(DataBeam).reshape(1,len(DataBeam)) 
    Beamparameter=pd.DataFrame(DataBeam,columns=columnsName) 
    
    return Beamparameter 
