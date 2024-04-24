import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

PTMC=np.genfromtxt("MCpT_pi.txt",delimiter=',',dtype=float)
PTTrk=np.genfromtxt("TrkpT_pi_old.txt",delimiter=',',dtype=float)
PTTrkNew=np.genfromtxt("TrkpT_pi_new.txt", delimiter=',', dtype=float)

for pT in PTMC:
    pT=float(pT)
for pT in PTTrk:
    pT=float(pT)
for pT in PTTrkNew:
    pT=float(pT)

bins=np.linspace(0, 150,50)
#MChist=plt.hist(PTMC,bins=bins,color="lightskyblue",label="MC Charged Pi from Tau")
#TRKHist=plt.hist(PTTrk,bins=bins,color="darkviolet",label="Trks, Charged Pi from Tau")
#TRKHistNew=plt.hist(PTTrkNew, bins=bins,color="fuchsia", label="NewTrks, Charged Pi from Tau")
MChist=np.histogram(PTMC,bins=bins)
TRKHist=np.histogram(PTTrk,bins=bins)
TRKHistNew=np.histogram(PTTrkNew,bins=bins)

MCstds=[]
TRKstds=[]
NewTRKstds=[]
for i in range(len(bins)-1):
    binvec=[]
    for j in range(len(PTMC)):
        if PTMC[j]>= bins[i] and PTMC[j]<bins[i+1]:
            binvec.append(PTMC[j])
    MCstds.append(np.std(binvec))
    binvec.clear()
    for j in range(len(PTTrk)):
        if PTTrk[j]>= bins[i] and PTTrk[j] < bins[i+1]:
            binvec.append(PTTrk[j])
    TRKstds.append(np.std(binvec))
    binvec.clear()
    for j in range(len(PTTrkNew)):
        if PTTrkNew[j]>= bins[i] and PTTrkNew[j] < bins[i+1]:
            binvec.append(PTTrkNew[j])
    NewTRKstds.append(np.std(binvec))
        

plt.xlabel("pT [GeV]")
#plt.legend()
plt.ylabel('$\epsilon$')
plt.title("Trk Efficiency for Charged Pions from Taus")
#plt.title("Charged Pion pT Spectrum, Trk and MC")
#plt.savefig("hist_compare.png")

npT_xaxis=[]
opT_xaxis=[]
nEff=[]
oEff=[]
nxerrs=[]
oxerrs=[]
nyerrs=[]
oyerrs=[]
#get the pT markers
for i in range(len(bins)-1):
    opT_xaxis.append(bins[i]+(bins[i+1]-bins[i])/2)
    npT_xaxis.append(opT_xaxis[i])
    oxerrs.append((bins[i+1]-bins[i])/2)
    nxerrs.append(oxerrs[i])


for i in range(len(MChist[0])):
    if MChist[0][i]!=0 and TRKHist[0][i]!=0:
        efficiency=TRKHist[0][i]/MChist[0][i]
        oEff.append(efficiency)
        oyerrs.append(efficiency*np.sqrt((TRKstds[i]/TRKHist[0][i])**2+(MCstds[i]/MChist[0][i])**2)) #add in quadrature
    elif MChist[0][i]==0 and TRKHist[0][i]==0:
        oEff.append(1)
        oyerrs.append(0)
    elif MChist[0][i]==0 and TRKHist[0][i]>0:
        opT_xaxis.pop(i)
        oxerrs.pop(i)
    elif TRKHist[0][i]==0 and MChist[0][i]>0:
        oEff.append(0)
        oyerrs.append(0)
    
    if MChist[0][i]!=0 and TRKHistNew[0][i]!=0:
        efficiency=TRKHistNew[0][i]/MChist[0][i]
        nEff.append(efficiency)
        nyerrs.append(efficiency*np.sqrt((NewTRKstds[i]/TRKHistNew[0][i])**2+(MCstds[i]/MChist[0][i])**2)) #add in quadrature
    elif MChist[0][i]==0 and TRKHistNew[0][i]==0:
        nEff.append(1)
        nyerrs.append(0)
    elif MChist[0][i]==0 and TRKHistNew[0][i]>0:
        npT_xaxis.pop(i)
        nxerrs.pop(i)
    elif TRKHistNew[0][i]==0 and MChist[0][i]>0:
        nEff.append(0)
        nyerrs.append(0)
        
#print(Eff)


plt.plot(opT_xaxis, oEff, 'o', color='b',markersize=2,label='$\pi\pm$ from taus (old)')
plt.errorbar(opT_xaxis,oEff, xerr=oxerrs, yerr=oyerrs, fmt = ' ', color='r')
plt.plot(npT_xaxis, nEff, 'o', color='orangered',markersize=2,label='$\pi\pm$ from taus (new)')
plt.errorbar(npT_xaxis,nEff, xerr=nxerrs, yerr=nyerrs, fmt = ' ', color='plum')
plt.axhline(y=1,xmin=0,xmax=100,color='darkgreen',ls='--',lw=2)
plt.legend()
plt.savefig("compare_piEff.png")


