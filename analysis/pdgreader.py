#read from list of all pdg's
reader=open('pdgFile.txt','r')
pdglist=[]
pdglist.clear()
for pdg in reader:
    if pdg not in pdglist:
        pdglist.append(pdg)

for i in range(len(pdglist)):
    pdglist[i]=pdglist[i].strip('\n')
   
#create a dictionary with all the pdgs and their rates

pdgRateDict = {}
with open('pdgFile.txt', 'r') as reader:
    pdglist_all=[]
    for line in reader:
        pdglist_all.append(line.strip('\n')) #remove newline character

#convert to ints and floats for math purposes
for i in range(len(pdglist_all)):
    pdglist_all[i]=float(pdglist_all[i])
for i in range(len(pdglist)):
    pdglist[i]=int(pdglist[i])

total=len(pdglist_all)

for pdg in pdglist:
    numParticles = 0
    rate=0
    for i in range(total):
        if pdglist_all[i]==pdg:
            numParticles +=1
    rate = numParticles/total
    pdgRateDict[pdg]=rate

print(pdgRateDict)
#now write the dict to a file
writer=open('pdgRates.txt', 'w')
for rate in pdgRateDict:
    writer.write(str(rate) +': '+str(pdgRateDict[rate])+'\n')

