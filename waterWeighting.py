
import sys,os,re
import subprocess as sp
import numpy as np
import fileinput
import math
from decimal import Decimal
import time

initPARAMS = np.loadtxt(sys.argv[1])
exptdata = np.loadtxt(sys.argv[2])
PARAMS = np.zeros((2,2), np.float64)
GAFFPARAMS = np.zeros((2,2), np.float64)
GAFFDIFFO = np.zeros((2,2), np.float64)
GAFFDIFFN = np.zeros((2,2), np.float64)
PARAMSNEW = np.zeros((2,2), np.float64)
ACOEFF = np.zeros((2,2), np.float64)
BCOEFF = np.zeros((2,2), np.float64)
NEWACOEFF = np.zeros((2,2), np.float64)
NEWBCOEFF = np.zeros((2,2), np.float64)
molrats = [0.00]
FRAMES = float(sys.argv[3])
U_z = np.zeros((len(molrats),FRAMES), np.float64)
U_z_dens = np.zeros(FRAMES, np.float64)
U_z_prime = np.zeros((len(molrats),FRAMES), np.float64)
U_z_prime_dens = np.zeros((len(molrats),FRAMES), np.float64)
U_z_diff = np.zeros((len(molrats),FRAMES), np.float64)
U_z_diffdens = np.zeros((len(molrats),FRAMES), np.float64)
U_z_rewe = np.zeros((len(molrats),FRAMES), np.float64)
U_z_tmp = np.zeros(FRAMES, np.float64)
U_z_tmp2 = np.zeros(FRAMES, np.float64)
STREDRE = np.zeros(len(molrats),np.float64)
STREDRD = np.zeros(len(molrats),np.float64)
STREDREe = np.zeros(len(molrats),np.float64)
STREDRDe = np.zeros(len(molrats),np.float64)
UZTMP = np.zeros(FRAMES, np.float64)
OLDSTORE = np.zeros(len(molrats),np.float64)
denstore = np.zeros(len(molrats),np.float64)
mixheatU_z = np.zeros(len(molrats),np.float64)
mixheatU_z_prime = np.zeros(len(molrats),np.float64)
mixheatU_z_er = np.zeros(len(molrats),np.float64)
varU_z = np.zeros(len(molrats),np.float64)
varU_z_prime = np.zeros(len(molrats),np.float64)
CVD = sys.argv[4]
scales = np.loadtxt(sys.argv[5])
GAFFDIFFAO = np.zeros(2, np.float64)
GAFFDIFFBO = np.zeros(2, np.float64)
GAFFDIFFAN = np.zeros(2, np.float64)
GAFFDIFFBN = np.zeros(2, np.float64)



def replc(NEWACOEFF, ACOEFF, NEWBCOEFF, BCOEFF, topopath):
   for i in range(2):
    for n in range(2):
      textToReplacea = "%.8E" % (NEWACOEFF[i,n])
      #print textToReplacea
      textToSearcha = "%.8E" % (ACOEFF[i,n])
      textToReplaceb = "%.8E" % (NEWBCOEFF[i,n])
      textToSearchb = "%.8E" % (BCOEFF[i,n])
      fileToSearch = topopath
      #print "Test"
      tempFile = open( fileToSearch, 'r+' )

      for line in fileinput.input( fileToSearch):
        tempFile.write (line.replace( textToSearcha, textToReplacea ) )
      tempFile.close()

      tempFile = open( fileToSearch, 'r+' )

      for line in fileinput.input( fileToSearch):
        tempFile.write (line.replace( textToSearchb, textToReplaceb ) )
      tempFile.close()

for q in range(len(molrats)):
  xpath = 'x'+str(molrats[q])
  #lmpstore = np.loadtxt(xpath+'/dens.dat')
  tmpstore = np.loadtxt(xpath+'/eptot.dat')
  for j in range(int(FRAMES)):
    U_z[q,j] = tmpstore[j]
  #for j in range(int(FRAMES)):
    #UZTMP[j] = U_z[q,j]
  #OLDSTORE[q] = np.mean(UZTMP)
  #lmpstore = np.loadtxt(xpath+'/dens.dat')
  #denstore[q] = np.mean(lmpstore)
  #print xpath+'run', OLDSTORE[q]


for i in range(2):
  PARAMS[i,0] = initPARAMS[i,0]
  #print PARAMS[i,0]
  PARAMS[i,1] = initPARAMS[i,1]
  #print PARAMS[i,1]

for i in range(2):
  GAFFPARAMS[i,0] = initPARAMS[i,0]
  #print PARAMS[i,0]
  GAFFPARAMS[i,1] = initPARAMS[i,1]

for i in range(len(molrats)):
  xpath = 'x'+str(molrats[i])
  sp.call(['cp','cpp.in',xpath])
  sp.call('mv full2.topo full.topo', cwd = xpath, shell = True)
  sp.call('cp full.topo full2.topo', cwd = xpath, shell = True)
  sp.call('cpptraj -i cpp.in', cwd = xpath, shell = True)

for l in range(len(scales)):
 for i in range(len(PARAMS[0])):
    for n in range(len(PARAMS[0])):
      sig = (PARAMS[i,0]+PARAMS[n,0])
      eps = np.sqrt(PARAMS[i,1]*PARAMS[n,1])
      ACOEFF[n,i] = (eps*(math.pow(sig,12)))
      #print ACOEFF[n,i]
      BCOEFF[n,i] = (eps*(math.pow(sig,6)))
      #print BCOEFF[n,i]
  #print "Assigned Old LJ Parameters"

 for i in range(2):
   if initPARAMS[i,2] == 1:
     PARAMSNEW[i,0] = abs(initPARAMS[i,0]+(initPARAMS[i,0]*scales[l,0]))
     print PARAMSNEW[i,0], scales[l,0], l
   if initPARAMS[i,2] == 0:
     PARAMSNEW[i,0] = initPARAMS[i,0]
     print PARAMSNEW[i,0]
   if initPARAMS[i,3] == 1:
     PARAMSNEW[i,1] = abs(initPARAMS[i,1]+(initPARAMS[i,1]*scales[l,1]))
     print PARAMSNEW[i,1], scales[l,1], l
   if initPARAMS[i,3] == 0:
     PARAMSNEW[i,1] = initPARAMS[i,1]
     print PARAMSNEW[i,1]

 for i in range(len(PARAMS[0])):
    for n in range(len(PARAMS[0])):
      sig = (PARAMSNEW[i,0]+PARAMSNEW[n,0])
      eps = np.sqrt(PARAMSNEW[i,1]*PARAMSNEW[n,1])
      NEWACOEFF[n,i] = (eps*(math.pow(sig,12)))
      #print NEWACOEFF[n,i]
      NEWBCOEFF[n,i] = (eps*(math.pow(sig,6)))
      #print NEWBCOEFF[n,i]
  #print "Assigned New LJ Parmaters"

 for q in range(len(molrats)):
      xpath = 'x'+str(molrats[q])
      topopath = 'x'+str(molrats[q])+'/full.topo'
      replc(NEWACOEFF, ACOEFF, NEWBCOEFF, BCOEFF, topopath)
      f=open('filename.sh', 'a')
      hoo = str("sander -O -i minsan.in -o "+xpath+"/minsan.out -p "+xpath+"/full.topo -c "+xpath+"/rst.001 -r
 new.rst -y "+xpath+"/traj.nc")
      f.write(hoo+'\n')
      f.close()
      print "Beginning ",molrats[q]
      sp.call('chmod u+x filename.sh', shell = True)
      sp.call('./filename.sh', shell = True)
      sp.call('rm -r filename.sh', shell = True)
      print "Done!"


 for q in range(len(molrats)):
   xpath = 'x'+str(molrats[q])
   sp.call('more minsan.out | grep ENE= | awk \'{print $4}\' > tmp.dat', cwd = xpath, shell = True)
   #sp.call('more minsan.out | grep ENE= | awk \'{print $4}\' > tmp.dat', cwd = xpath, shell = True)
   tmpstore = np.loadtxt(xpath+'/tmp.dat')
   lmpstore = np.loadtxt(xpath+'/dens.dat')

   for j in range(int(FRAMES)):
      U_z_prime[q,j] = tmpstore[j]
      U_z_prime_dens[q,j] = lmpstore[j]

   for b in range(int(FRAMES)):
      U_z_diff[q,b] = U_z_prime[q,b]*np.exp((U_z[q,b]-U_z_prime[q,b])*(1/(-8.314*300)))
      U_z_diffdens[q,b] = U_z_prime_dens[q,b]*np.exp((U_z[q,b]-U_z_prime[q,b])*(1/(-8.314*300)))
      U_z_rewe[q,b] = np.exp((U_z[q,b]-U_z_prime[q,b])*(1/(-8.314*300)))
      U_z_tmp[b] = U_z_diff[q,b]
      U_z_dens[b] = U_z_diffdens[q,b]
      U_z_tmp2[b] = U_z_rewe[q,b]

   STREDRE[q] = (np.mean(U_z_tmp))/(np.mean(U_z_tmp2))
   STREDRD[q] = (np.mean(U_z_dens))/(np.mean(U_z_tmp2))
   #print STREDRE[q], STREDRD[q]

 for w in range(len(molrats)):
  xpath = 'x'+str(molrats[w])
  topopath = 'x'+str(molrats[w])+'/full.topo'
  #replc(ACOEFF, NEWACOEFF, BCOEFF, NEWBCOEFF, topopath)
  #sp.call('pmemd.cuda -O -p full.topo -c full.crds -i mini.in -o mini2.out -r mini2.rst -inf /dev/null', cwd=
xpath, env=dict(os.environ, CUDA_VISIBLE_DEVICES=str(CVD)), shell=True)
  sp.call('pmemd.cuda -O -p full.topo -c mini.rst -i eq1.in -o eq12.out -r eq12.rst -e eq12.mden -inf eq12.inf
',cwd=xpath, env=dict(os.environ, CUDA_VISIBLE_DEVICES=str(CVD)), shell=True)
  sp.call('pmemd.cuda -O -p full.topo -c eq12.rst -i eq2.in -o eq22.out -r eq22.rst -e eq22.mden -inf eq22.inf
',cwd=xpath, env=dict(os.environ, CUDA_VISIBLE_DEVICES=str(CVD)), shell=True)
  sp.call('pmemd.cuda -O -p full.topo -c eq22.rst -i mdin -o mdout2.001 -r rst2.001 -x traj2.001 -e mden2.001
-inf mdinfo2.001',cwd=xpath, env=dict(os.environ, CUDA_VISIBLE_DEVICES=str(CVD)), shell=True)
  sp.call('grep L3 mden2.001 | tail -n +2 | awk \'{printf "%.5f\\n",$2}\' > vol2.dat', cwd=xpath, shell=True)
  sp.call('cat vol2.dat | python /home/henrikse/pe/reblock.py >& vol2.err.dat', cwd=xpath, shell=True)
  sp.call('grep L6 mden2.001 | tail -n +2 | awk \'{printf "%.5f\\n",$3}\' > eptot2.dat', cwd=xpath, shell=True
)
  sp.call('cat eptot2.dat | python /home/henrikse/pe/reblock.py >& eptot2.err.dat', cwd=xpath, shell=True)
  sp.call('grep L9 mden2.001 | tail -n +2 | awk \'{printf "%.5f\\n",$5}\' > dens2.dat', cwd=xpath, shell=True)
  sp.call('cat dens2.dat | python /home/henrikse/pe/reblock.py >& dens2.err.dat', cwd=xpath, shell=True)
  sp.call('more dens2.err.dat | grep BlkAvg | awk \'{print $7}\' > tmp.dat', cwd = xpath, shell = True)
  denstore[w] = np.loadtxt(xpath+'/tmp.dat')
  sp.call('more dens2.err.dat | grep BlkAvg | awk \'{print $9}\' > tmp.dat', cwd = xpath, shell = True)
  STREDRDe[q] = np.loadtxt(xpath+'/tmp.dat')
  sp.call('more eptot2.err.dat | grep BlkAvg | awk \'{print $7}\' > tmp.dat', cwd = xpath, shell = True)
  OLDSTORE[w] = np.loadtxt(xpath+'/tmp.dat')
  sp.call('more eptot2.err.dat | grep BlkAvg | awk \'{print $9}\' > tmp.dat', cwd = xpath, shell = True)
  STREDREe[q] = np.loadtxt(xpath+'/tmp.dat')

 for i in range(len(molrats)):
   print "MD DENSITY = ", denstore[i], "MD ENERGY per WATER = ", OLDSTORE[i]/1000, "MD ENERGY ERR per WATER =
", STREDREe[i]/1000, "REWEIGHTING DENSITY = ",  STREDRD[i], "REWEIGHTING ENERGY per WATER = ", STREDRE[i]/1000
 for w in range(len(molrats)):
  xpath = 'x'+str(molrats[w])
  sp.call('cp full2.topo full.topo', cwd = xpath, shell = True)
