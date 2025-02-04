## Function name: popt_planck
## Purpose: This function loads the PLANCK dust emission polarization map and calculates the expected optical polarization (V-Band) based on the ratio from:
##       Planck 2018 (XII) - A&A 641, A12 (2020)  // Planck Collaboration Int. XXI. 2015, A&A, 576, A106
## Input:
#     ra and dec of target in degrees 
## Optional input:
#     ebvside        size (arcmin) of field to get error on EBV from Schlafly [def: 7' from FORS]



from astropy.coordinates import SkyCoord
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os

def get_planck(ra,dec,planckfile=None):
       
    ## --- Read map
    if planckfile is None:
        home = os.path.expanduser('~')
        planckdir='/media/joaomfras/GalaxyPol/Pol-Gal/PLANCK/'#home+
        planckfile = planckdir+'COM_CompMap_IQU-thermaldust-gnilc-varres_2048_R3.00.fits'#unires/varres
    plmap = hp.read_map(planckfile,field=0)#(0,1,2,3,6,8))#I,Q,U,erI*2,erQ*2,erU*2

    ## --- Mollview
    # hp.mollview(plmap)

    ## --- Coordinates (galactic): pix to ang
    nside = hp.get_nside(plmap)
    npix = hp.nside2npix(nside)
    #glons, glats = hp.pix2ang(nside, np.arange(npix), lonlat=True)

    ## --- Resolution
    #resol = hp.nside2resol(nside, arcmin=True)

    ## --- Coordinates: ang to pix
    pos = SkyCoord(ra,dec,frame='fk5',unit='deg')
    pix = hp.ang2pix(nside,pos.galactic.l.deg,pos.galactic.b.deg,lonlat=True)

    ## -- individual values to save memory
    intens = plmap[pix]
    del plmap
    qmap = hp.read_map(planckfile,field=1)
    q = qmap[pix]
    del qmap
    umap = hp.read_map(planckfile,field=2)
    u = umap[pix]
    del umap
    erimap = hp.read_map(planckfile,field=3)
    eri = np.sqrt(erimap[pix])
    del erimap
    erqmap = hp.read_map(planckfile,field=6)
    erq = np.sqrt(erqmap[pix])
    del erqmap
    erumap = hp.read_map(planckfile,field=8)
    eru = np.sqrt(erumap[pix])
    del erumap
    
    return q,erq,u,eru,intens,eri

### This gets from q/u directly from Fig. 7 right, Eq. 15 (planck15) 
def quopt_planck(ra,dec,planckfile=None,corr_ang=-3.1):

    ## get q,u,i from Planck
    q,erq,u,eru,intens,eri = get_planck(ra,dec,planckfile=planckfile)

    ## conversion from K_CMB to MJy/sr (see https://irsa.ipac.caltech.edu/docs/knowledgebase/ercsc-validation-conversions.pdf)
    fac = 296.877

    q,erq,u,eru,intens,eri = fac*q,fac*erq,fac*u,fac*eru,fac*intens,fac*eri
    
    ## pol and error
    #p = np.sqrt(q**2+u**2)
    #sig_p = np.sqrt(((q*erq)**2 + (u*eru)**2)/(q**2+u**2))

    ## fit slopes /intercepts (MJy/sr)
    sl,ersl = -5.4,np.sqrt(0.3**2+0.2**2)#-5.32,0.06
    interc,erinterc = 0.002,0.0009
    qopt = (q - interc)/sl
    uopt = (u - interc)/sl
    erqopt = np.sqrt((erq/sl)**2 + (erinterc/sl)**2 + ((q-interc)*ersl/sl**2)**2)
    eruopt = np.sqrt((eru/sl)**2 + (erinterc/sl)**2 + ((q-interc)*ersl/sl**2)**2)

    ## get p and angle
    popt = np.sqrt(qopt**2+uopt**2)
    erpopt = np.sqrt((qopt**2*erqopt**2+uopt**2*eruopt**2)/\
                     (qopt**2+uopt**2))
    print("p = %f +/- %f" %(popt,erpopt))
    angopt = 0.5*np.arctan2(-uopt,qopt)/np.pi*180
    
    if corr_ang is not None:
        angopt += corr_ang
    erangopt = 0.5*np.sqrt((qopt*eruopt)**2+(uopt*erqopt)**2)/\
            ((1+(uopt/qopt)**2)*qopt**2)/np.pi*180
    print("ang(deg) = %f +/- %f" %(angopt,erangopt))

     
    return popt,erpopt,angopt,erangopt
    #return qopt,erqopt,uopt,eruopt
    
### This gets from R relation of P_submm/P_opt and EBV instead of q/u
def popt_planck(ra,dec,ebv=None,ebverr=None,planckfile=None,ebvside=7,corr_ang=-3.1):

    pos = SkyCoord(ra,dec,frame='fk5',unit='deg')
    
    ## get q,u,i from Planck
    q,erq,u,eru,intens,eri = get_planck(ra,dec,planckfile=planckfile)
    
    ## pol and error
    p = np.sqrt(q**2+u**2)
    sig_p = np.sqrt(((q*erq)**2 + (u*eru)**2)/(q**2+u**2))

    ## tau_v = Av/1.086 = EBV*Rv/1.086
    if ebv is None:
        from dustmaps.sfd import SFDQuery
        sfd = SFDQuery()
        ebv = sfd(pos)
        nside = ebvside/1 #one per arcmin
        ar_ra = np.arange(nside)*0.0166667-0.5*nside*0.0166667+ra
        ar_dec = np.arange(nside)*0.0166667-0.5*nside*0.0166667+dec
        m_ra,m_dec = np.meshgrid(ar_ra,ar_dec)
        m_pos = SkyCoord(m_ra,m_dec,frame='fk5',unit='deg')
        m_ebv = sfd(m_pos)
        ebverr = np.std(m_ebv)
        print("EBV = %f +/- %f" %(ebv,ebverr))
    Rv = 3.1
    tau_v = ebv*Rv/1.086
    sig_tau_v = ebverr*Rv/1.086
    
    ## Ratio of Planck/optical = 4.2 +/- 0.5 /// 4.31 +/- 0.04
    ratio,sig_ratio = 4.2,np.sqrt(0.3**2+0.2**2)#4.31,0.04
    
    p_opt = (p/intens)/(ratio/tau_v)
    sig_p_opt = np.sqrt(((tau_v*sig_p)/(intens*ratio))**2 + ((p*sig_tau_v)/(intens*ratio))**2 + \
                    ((p*tau_v*eri)/(intens**2*ratio))**2 + ((p*tau_v*sig_ratio)/(intens*ratio**2))**2)
    print("p = %f +/- %f" %(p_opt,sig_p_opt))

    ## Rotate pol angle 90deg
    angle_opt = 0.5*np.arctan2(-u,q)/np.pi*180 + 90
    sig_angle_opt = 0.5*np.sqrt((q*eru)**2+(u*erq)**2)/\
                    ((1+(u/q)**2)*q**2)/np.pi*180
    if corr_ang is not None:
        angle_opt += corr_ang
    print("ang(deg) = %f +/- %f" %(angle_opt,sig_angle_opt))
        
    return (p_opt,sig_p_opt,angle_opt,sig_angle_opt)

