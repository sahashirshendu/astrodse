import numpy as np
#  This program will calculate a static stellar model using the
#  stellar structure equations.  The user is expected to supply the
#  star's mass, luminosity, effective temperature, and composition
#  (X and Z).  If the choices for these quantities are not consistent
#  with the central boundary conditions, an error message will be
#  generated and the user will then need to supply a different set of
#  initial values.
def STARTMDL(deltar, X, Z, mu, Rs, r_i, M_ri, L_ri, tog_bf,irc,cst):
      r = r_i + deltar
      M_rip1 = M_ri
      L_rip1 = L_ri
      if (irc == 0):
          T_ip1 = cst.G*M_rip1*mu*cst.m_H/(4.25e0*cst.k_B)*(1.0e0/r - 1.0e0/Rs)
          A_bf = 4.34e25*Z*(1.0e0 + X)/tog_bf
          A_ff = 3.68e22*cst.g_ff*(1.0e0 - Z)*(1.0e0 + X)
          Afac = A_bf + A_ff
          P_ip1 = np.sqrt((1.0e0/4.25e0)*(16.0e0/3.0e0*np.pi*cst.a*cst.c)*(cst.G*M_rip1/L_rip1)*(cst.k_B/(Afac*mu*cst.m_H)))*T_ip1**4.25e0
      else:
          T_ip1 = cst.G*M_rip1*mu*cst.m_H/cst.k_B*(1.0e0/r - 1.0e0/Rs)/cst.gamrat
          P_ip1 = cst.kPad*T_ip1**cst.gamrat
      return r,P_ip1, M_rip1, L_rip1, T_ip1
def   EOS(X, Z, XCNO, mu, P, T,izone,cst):
      if ((T < 0.0e0) or (P < 0.0e0)):
          print(' Something is a little wrong here.')
          print(' You are asking me to deal with either a negative temperature')
          print(' or a negative pressure.  I am sorry but that is not in my')
          print(' contract! You will have to try again with different') 
          print(' initial conditions.')
          print(' In case it helps, I detected the problem in zone ')
          print(' with the following conditions:')
          print( 'T = ',T,' K') 
          print('P = ',P,' dynes/cm**2')
          return (0.0, 0.0, 0.0,0.0,1)
      Prad = cst.a*T**4/3.0e0
      Pgas = P - Prad
      rho=(mu*cst.m_H/cst.k_B)*(Pgas/T)
      if (rho < 0.0e0):
          print(' I am sorry, but a negative density was detected.')
          print(' my equation-of-state routine is a bit baffled by this new')
          print(' physical system you have created.  The radiation pressure')
          print(' is probably too great, implying that the star is unstable.')
          print(' Please try something a little less radical next time.')
          print('In case it helps, I detected the problem in zone ')
          print(' with the following conditions:')
          print('T = {0:12.5E} K'.format(T))
          print('P_total = {0:12.5E} dynes/cm**2'.format(P))
          print('P_rad   = {0:12.5E} dynes/cm**2'.format(Prad))
          print('P_gas= {0:12.5E} dynes/cm**2'.format(Pgas))
          print('rho = {0:12.5E} g/cm**3'.format(rho))
          return (rho, 0.0 , 0.0 ,0.0,1)
      tog_bf = 2.82e0*(rho*(1.0e0 + X))**0.2e0
      k_bf = 4.34e25/tog_bf*Z*(1.0e0 + X)*rho/T**3.5e0
      k_ff = 3.68e22*cst.g_ff*(1.0e0 - Z)*(1.0e0 + X)*rho/T**3.5e0
      k_e = 0.2e0*(1.0e0 + X)
      kappa = k_bf + k_ff + k_e
      oneo3=0.333333333e0
      twoo3=0.666666667e0
      T6 = T*1.0e-06  
      fx = 0.133e0*X*np.sqrt((3.0e0 + X)*rho)/T6**1.5e0
      fpp = 1.0e0 + fx*X
      psipp = 1.0e0 + 1.412e8*(1.0e0/X - 1.0e0)*np.exp(-49.98*T6**((-1.0)*oneo3))
      Cpp = 1.0e0 + 0.0123e0*T6**oneo3 + 0.0109e0*T6**twoo3 + 0.000938e0*T6
      epspp = 2.38e6*rho*X*X*fpp*psipp*Cpp*T6**(-twoo3)*np.exp(-33.80e0*T6**(-oneo3))
      CCNO = 1.0e0 + 0.0027e0*T6**oneo3 - 0.00778e0*T6**twoo3- 0.000149e0*T6
      epsCNO = 8.67e27*rho*X*XCNO*CCNO*T6**(-twoo3)*np.exp(-152.28e0*T6**(-oneo3))
      epslon = epspp + epsCNO
      return (rho, kappa, epslon,tog_bf,0)
def dPdr(r, M_r, rho,cst):
      return -cst.G*rho*M_r/r**2
def dMdr(r, rho, cst):
      return (4.0e0*np.pi*rho*r**2)
def dLdr(r, rho, epslon,cst):
      return (4.0e0*np.pi*rho*epslon*r**2)
def dTdr(r, M_r, L_r, T, rho, kappa, mu, irc,cst):
      if (irc == 0):
          return (-(3.0e0/(16.0e0*np.pi*cst.a*cst.c))*kappa*rho/T**3*L_r/r**2)
      else:
          return (-1.0e0/cst.gamrat*cst.G*M_r/r**2*mu*cst.m_H/cst.k_B)
def RUNGE(f_im1, dfdr, r_im1, deltar, irc, X, Z, XCNO,mu, izone,cst):
      f_temp=np.zeros(4)
      f_i=np.zeros(4)
      dr12 = deltar/2.0e0
      dr16 = deltar/6.0e0
      r12  = r_im1 + dr12
      r_i  = r_im1 + deltar
      for i in range(0,4):
          f_temp[i] = f_im1[i] + dr12*dfdr[i]
      df1, ierr = FUNDEQ(r12, f_temp, irc, X, Z, XCNO, mu, izone,cst)
      if (ierr != 0):
          return f_i,ierr
      for i in range(0,4):
          f_temp[i] = f_im1[i] + dr12*df1[i]
      df2,ierr = FUNDEQ(r12, f_temp, irc, X, Z, XCNO, mu, izone,cst)
      if (ierr != 0):
          return f_i,ierr
      for i in range(0,4):
          f_temp[i] = f_im1[i] + deltar*df2[i]
      df3,ierr=FUNDEQ(r_i, f_temp, irc, X, Z, XCNO, mu, izone,cst)
      if (ierr != 0):
          return f_i,ierr
      for i in range(0,4):
          f_i[i] = f_im1[i] + dr16*(dfdr[i] + 2.0e0*df1[i] + 2.0e0*df2[i] + df3[i])
      return f_i,0
def FUNDEQ(r,f,irc,X,Z,XCNO,mu,izone,cst):
      dfdr=np.zeros(4)
      P   = f[0]
      M_r = f[1]
      L_r = f[2]
      T   = f[3]
      rho,kappa,epslon,tog_bf,ierr = EOS(X, Z, XCNO, mu, P, T, izone,cst)
      dfdr[0] = dPdr(r, M_r, rho,cst)
      dfdr[1] = dMdr(r, rho,cst)
      dfdr[2] = dLdr(r, rho, epslon,cst)
      dfdr[3] = dTdr(r, M_r, L_r, T, rho, kappa, mu, irc,cst)
      return (dfdr,ierr)
def StatStar(Msolar,Lsolar,Te,X,Z):
      nsh=999
      r=np.zeros(nsh,float)
      P=np.zeros(nsh,float)
      M_r=np.zeros(nsh,float)
      L_r=np.zeros(nsh,float)
      T=np.zeros(nsh,float)
      rho=np.zeros(nsh,float)
      kappa=np.zeros(nsh,float)
      epslon=np.zeros(nsh,float)
      tog_bf=np.zeros(nsh,float)
      dlPdlT=np.zeros(nsh,float)
      deltar=0.0
      XCNO=0.0
      mu=0.0
      Ms=0.0
      Ls=0.0
      Rs=0.0
      T0=0.0
      P0=0.0
      Pcore=0.0   
      Tcore=0.0   
      rhocor=0.0  
      epscor=0.0  
      rhomax=0.0  
      Rsolar=0.0
      Y=1.0-(X+Z)
      tog_bf0=0.01
      f_im1=np.zeros(4,float)
      dfdr=np.zeros(4,float)
      f_i=np.zeros(4,float)
      Nstart=10
      Nstop=999
      Igoof=-1
      ierr=0
      P0=0.0
      T0=0.0
      dlPlim=99.9
      debug=0
      Rsun=6.9599e10
      Msun=1.989e33
      Lsun=3.826e33
      class Constants:
           pass
      cst = Constants()
      cst.sigma=5.67051e-5
      cst.c=2.99792458e10
      cst.a=7.56591e-15
      cst.G=6.67259e-8
      cst.k_B=1.380658e-16
      cst.m_H=1.673534e-24 
      cst.gamma= 5.0e0/3
      cst.g_ff= 1.0e0   
      XCNO = Z/2.0e0
      Ms = Msolar*Msun
      Ls = Lsolar*Lsun
      Rs = np.sqrt(Ls/(4.e0*np.pi*cst.sigma))/Te**2
      Rsolar = Rs/Rsun
      deltar = -Rs/1000.0e0
      idrflg = 0
      mu = 1.0e0/(2.0*X + 0.75*Y + 0.5*Z)
      cst.gamrat = cst.gamma/(cst.gamma - 1.0e0)
      initsh=0
      r[initsh]   = Rs
      M_r[initsh] = Ms
      L_r[initsh] = Ls
      T[initsh]   = T0
      P[initsh]   = P0
      tog_bf[initsh]   = tog_bf0
      if (P0 <= 0.0) or (T0 <= 0.0):
          rho[initsh]    = 0.0
          kappa[initsh]  = 0.0
          epslon[initsh] = 0.0
          tog_bf[initsh] = 0.01
      else:
          rho[initsh],kappa[initsh],epslon[initsh],tog_bf[initsh],ierr=EOS(X, Z, XCNO, mu, P[initsh], T[initsh], 0 ,cst)
          if ierr != 0:
              print ("we're stopping now")
              istop=0
      cst.kPad = 0.3e0
      irc = 0
      dlPdlT[initsh] = 4.25e0
      for i in range(0,Nstart):
          ip1 = i + 1
          r[ip1],P[ip1],M_r[ip1],L_r[ip1],T[ip1]=STARTMDL(deltar, X, Z, mu, Rs, r[i], M_r[i], L_r[i], tog_bf[i], irc,cst)
          rho[ip1],kappa[ip1],epslon[ip1],tog_bf[ip1],ierr=EOS(X, Z, XCNO, mu, P[ip1], T[ip1], ip1 ,cst)
          if ierr != 0:
              print('Values from the previous zone are:')
              print('r/Rs      = {0:12.5E}'.format(r[i]/Rs))
              print('rho       = {0:12.5E}  g/cm**3'.format(rho[i]))
              print('M_r/Ms    = {0:12.5E}'.format(M_r[i]/Ms))
              print('kappa     = {0:12.5E}  cm**2/g'.format(kappa[i]))
              print('T         = {0:12.5E}  K'.format(T[i]))
              print('epsilon   = {0:12.5E}  ergs/g/s'.format(epslon[i]))
              print('P         = {0:12.5E}  dynes/cm**2'.format(P[i]))
              print('L_r/Ls    = {0:12.5E}'.format(L_r[i]/Ls))
              break
          if (i > initsh):
              dlPdlT[ip1] = np.log(P[ip1]/P[i])/np.log(T[ip1]/T[i])
          else:
              dlPdlT[ip1] = dlPdlT[i]
          if (dlPdlT[ip1] < cst.gamrat):
              irc = 1
          else:
              irc = 0
              cst.kPad = P[ip1]/T[ip1]**cst.gamrat
          deltaM = deltar*dMdr(r[ip1], rho[ip1],cst)
          M_r[ip1] = M_r[i] + deltaM
          if (np.abs(deltaM) > (0.001e0*Ms)):
              if (ip1 > 1):
                  ip1 = ip1 - 1
                  print(' The variation in mass has become larger than 0.001*Ms')
                  print(' leaving the approximation loop before Nstart was reached')
                  break
      Nsrtp1 = ip1 + 1
      if (ierr != 0):
          Nstop=Nsrtp1-1
          istop=Nstop
      for i in range(Nsrtp1,Nstop):
          im1 = i - 1
          f_im1[0] = P[im1]
          f_im1[1] = M_r[im1]
          f_im1[2] = L_r[im1]
          f_im1[3] = T[im1]
          dfdr[0]  = dPdr(r[im1], M_r[im1], rho[im1],cst)
          dfdr[1]  = dMdr(r[im1], rho[im1],cst)
          dfdr[2]  = dLdr(r[im1], rho[im1], epslon[im1],cst)
          dfdr[3]  = dTdr(r[im1], M_r[im1], L_r[im1], T[im1], rho[im1],kappa[im1], mu, irc,cst)
          f_i,ierr=RUNGE(f_im1, dfdr, r[im1], deltar, irc, X, Z, XCNO, mu, i,cst)
          if (ierr != 0):
              print(' The problem occurred in the Runge-Kutta routine')
              print(' Values from the previous zone are:')
              print('r/Rs    = {0:12.5e}'.format(r[im1]/Rs))
              print('rho     = {0:12.5e} g/cm**3'.format(rho[im1]))
              print('M_r/Ms  = {0:12.5e}'.format(M_r[im1]/Ms))
              print('kappa   = {0:12.5e} cm**2/g'.format(kappa[im1]))
              print('T       = {0:12.5e} K'.format(T[im1]))
              print('epsilon = {0:12.5e} ergs/g/s'.format(epslon[im1]))
              print('P       = {0:12.5e} dynes/cm**2'.format(P[im1]))
              print('L_r/Ls  = {0:12.5e}'.format(L_r[im1]/Ls))
              break
          r[i]   = r[im1] + deltar
          P[i]   = f_i[0]
          M_r[i] = f_i[1]
          L_r[i] = f_i[2]
          T[i]   = f_i[3]
          rho[i],kappa[i],epslon[i],tog_bf[i],ierr=EOS(X, Z, XCNO, mu, P[i], T[i], i, cst)
          if (ierr != 0):           
              print( ' Values from the previous zone are:')
              print( 'r/Rs    = {0:12.5e}'.format(r[im1]/Rs))
              print( 'rho     = {0:12.5e} g/cm**3'.format(rho[im1]))
              print( 'M_r/Ms  = {0:12.5e}'.format(M_r[im1]/Ms))
              print('kappa   = {0:12.5e} cm**2/g'.format(kappa[im1]))
              print('T       = {0:12.5e} K'.format(T[im1]))
              print('epsilon = {0:12.5e} ergs/g/s'.format(epslon[im1]))
              print('P       = {0:12.5e} dynes/cm**2'.format(P[im1]))
              print('L_r/Ls  = {0:12.5e}'.format(L_r[im1]/Ls))
              istop = i
              break
          if (debug == 1):
            print (i,r[i],M_r[i],L_r[i],T[i],P[i],rho[i],kappa[i],epslon[i],tog_bf[i])
          dlPdlT[i] = np.log(P[i]/P[im1])/np.log(T[i]/T[im1])
          if (dlPdlT[i] < cst.gamrat):
              irc = 1
          else:
              irc = 0
          if ((r[i] <= np.abs(deltar)) and ((L_r[i] >= (0.1e0*Ls)) or (M_r[i] >= (0.01e0*Ms)))):
              Igoof = 6   
          elif (L_r[i] <= 0.0e0):
              Igoof = 5   
              rhocor = M_r[i]/(4.0e0/3.0e0*np.pi*r[i]**3)               
              if (M_r[i] != 0.0e0):
                  epscor = L_r[i]/M_r[i]                               
              else:
                  epscor = 0.0e0
              Pcore = P[i] + 2.0e0/3.0e0*np.pi*cst.G*rhocor**2*r[i]**2  
              Tcore = Pcore*mu*cst.m_H/(rhocor*cst.k_B)                 
          elif (M_r[i] <= 0.0e0):
              Igoof  = 4
              Rhocor = 0.0e0
              epscor = 0.0e0
              Pcore  = 0.0e0
              Tcore  = 0.0e0
          elif ((r[i] < (0.02e0*Rs)) and ((M_r[i] < (0.01e0*Ms)) and ((L_r[i] < 0.1e0*Ls)))):
              rhocor = M_r[i]/(4./3.*np.pi*r[i]**3)                 
              rhomax = 10.0e0*(rho[i]/rho[im1])*rho[i]              
              epscor = L_r[i]/M_r[i]          
              Pcore  = P[i] + 2.0e0/3.0e0*np.pi*cst.G*rhocor**2*r[i]**2   
              Tcore  = Pcore*mu*cst.m_H/(rhocor*cst.k_B)                  
              if ((rhocor < rho[i]) or (rhocor > rhomax)):          
                  Igoof = 1                             
              elif (epscor < epslon[i]):
                  Igoof = 2                           
              elif (Tcore < T[i]):
                  Igoof = 3                            
              else:
                  Igoof = 0
          if (Igoof != -1):
              istop = i
              break
          if ((idrflg == 0) and (M_r[i] < (0.99e0*Ms))):
               deltar = (-1.0)*Rs/100.0e0
               idrflg = 1
          if ((idrflg == 1) and (deltar >= (0.5*r[i]))):
               deltar = (-1.0)*Rs/5000.0e0
               idrflg = 2
          istop = i
      rhocor = M_r[istop]/(4.0e0/3.0e0*np.pi*r[istop]**3)
      epscor = L_r[istop]/M_r[istop]
      Pcore  = P[istop] + 2.0e0/3.0e0*np.pi*cst.G*rhocor**2*r[istop]**2
      Tcore  = Pcore*mu*cst.m_H/(rhocor*cst.k_B)
      if  (Igoof != 0):
          if (Igoof == -1):
              print('Sorry to be the bearer of bad news, but...')
              print('Your model has some problems')
              print('The number of allowed shells has been exceeded')
          if (Igoof == 1):
              print('It looks like you are getting close,')
              print('however, there are still a few minor errors')
              print('The core density seems a bit off')
              print(' density should increase smoothly toward the center.')
              print(' The density of the last zone calculated was rho = ',rho[istop],' gm/cm**3')
              print(rhocor,rhomax)
          if (rhocor > 1e10):
              print('It looks like you will need a degenerate')
              print(' neutron gas and general relativity')
              print(' to solve this core.  Who do you think I am, Einstein?')
          if (Igoof == 2):
              print('It looks like you are getting close')
              print('however, there are still a few minor errors')
              print('The core epsilon seems a bit off')
              print(' epsilon should vary smoothly near the center.')
              print(' The value calculated for the last zone was eps =',epslon[istop],' ergs/g/s')
          if (Igoof == 3):
              print('It looks like you are getting close,')
              print('however, there are still a few minor errors')
              print(' Your extrapolated central temperature is too low')
              print(' a little more fine tuning ought to do it.')
              print(' The value calculated for the last zone was T = ',T[istop],' K')
          if (Igoof == 4):
              print('Sorry to be the bearer of bad news, but...')
              print('Your model has some problems')
              print('You created a star with a hole in the center!')
          if (Igoof == 5):
              print('Sorry to be the bearer of bad news, but...')
              print('Your model has some problems')
              print('This star has a negative central luminosity!')
          if (Igoof == 6):
              print('Sorry to be the bearer of bad news, but...')
              print('Your model has some problems')
              print('You hit the center before the mass and/or ')
              print('luminosity were depleted!')
      else:
          print('CONGRATULATIONS, I THINK YOU FOUND IT!')
          print('However, be sure to look at your model carefully.')
      Rcrat = r[istop]/Rs
      if (Rcrat < -9.999e0): Rcrat = -9.999e0
      Mcrat = M_r[istop]/Ms
      if (Mcrat < -9.999e0): Mcrat = -9.999e0
      Lcrat = L_r[istop]/Ls
      if (Lcrat < -9.999e0): Lcrat = -9.999e0
      f=open('starmodl.dat','w')
      f.write('A Homogeneous Main-Sequence Model\n')
      f.write(' The surface conditions are:        The central conditions are:\n')
      f.write(' Mtot = {0:13.6E} Msun          Mc/Mtot     = {1:12.5E}\n'.format(Msolar,Mcrat))
      f.write(' Rtot = {0:13.6E} Rsun          Rc/Rtot     = {1:12.5E}\n'.format(Rsolar,Rcrat))
      f.write(' Ltot = {0:13.6E} Lsun          Lc/Ltot     = {1:12.5E}\n'.format(Lsolar,Lcrat))
      f.write(' Teff = {0:13.6E} K             Density     = {1:12.5E}\n'.format(Te,rhocor))
      f.write(' X    = {0:13.6E}               Temperature = {1:12.5E}\n'.format(X,Tcore))
      f.write(' Y    = {0:13.6E}               Pressure    = {1:12.5E} dynes/cm**2\n'.format(Y,Pcore))
      f.write(' Z    = {0:13.6E}               epsilon     = {1:12.5E} ergs/s/g\n'.format(Z,epscor))
      f.write('                                    dlnP/dlnT   = {0:12.5E}\n'.format(dlPdlT[istop]))
      f.write('Notes:\n')
      f.write(' (1) Mass is listed as Qm = 1.0 - M_r/Mtot, where Mtot = {0:13.6}\n'.format(Msun))
      f.write(' (2) Convective zones are indicated by c, radiative zones by r\n')
      f.write(' (3) dlnP/dlnT may be limited to +99.9 or -99.9# if so it is\n')
      f.write(' labeled by *\n')
      f.write('   r        Qm       L_r       T        P        rho      kap      eps     dlPdlT\n')
      for ic in range(0,istop+1):
          i = istop - ic 
          Qm = 1.0e0 - M_r[i]/Ms    # Total mass fraction down to radius
          if (dlPdlT[i] < cst.gamrat):
              rcf = 'c'
          else:
              rcf = 'r'
          if (np.abs(dlPdlT[i]) > dlPlim):
              dlPdlT[i] = np.copysign(dlPlim,dlPdlT[i])   
              clim = '*'
          else:
              clim = ' '
          s='{0:7.2E} {1:7.2E} {2:7.2E} {3:7.2E} {4:7.2E} {5:7.2E} {6:7.2E} {7:6.2E}{8:1s}{9:1s} {10:5.1f}\n'.format(r[i], Qm, L_r[i], T[i], P[i], rho[i], kappa[i],epslon[i], clim, rcf, dlPdlT[i])  
          f.write(s)
      print ('***** The integration has been completed *****')
      print ('      The model has been stored in starmodl.dat')
      return Igoof,ierr,istop
def main():
      getinp=1  # read in input
      if (getinp == 1):
           Msolar=float(input(' Enter the mass of the star (in solar units):'))
           Lsolar=float(input(' Enter the luminosity of the star (in solar units):'))
           Te=float(input(' Enter the effective temperature of the star (in K):'))
           Y=-1.0
           while (Y < 0.0):
               X=float(input(' Enter the mass fraction of hydrogen (X):'))
               Z=float(input(' Enter the mass fraction of metals (Z):'))
               Y = 1.e0 - X - Z
               if Y < 0:
                     print('You must have X + Z <= 1. Please reenter composition.')
      Igoof,ierr,istop=StatStar(Msolar,Lsolar,Te,X,Z)
main()
