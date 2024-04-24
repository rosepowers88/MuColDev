#read in from the efficiency file, split the lines, plot efficiency vs pT
import matplotlib
import matplotlib.pyplot as plt
pT=[]
eff=[]
reader=open('EffFile.txt', 'r')
for line in reader:
    p_eff = line.split(',')
    p=float(p_eff[0])
    pT.append(p)
    eff.append(p_eff[1])


plt.plot(pT, eff, 'o',color='k',markersize='2')
plt.xlabel('pT [GeV]')
plt.ylabel('Tau Efficiency')
plt.savefig('efficiency.png')
