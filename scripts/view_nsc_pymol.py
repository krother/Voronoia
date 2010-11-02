
#
# view_nsc_pymol.py
# View NSC structures with PyMOL
#
# usage: run view_nsc_pymol.py (from PyMOL)
# new commands:
#    nsc <filename> - displays a sco file
#    nsc <filename> <number> - displays a patch 
#    nsc <filename> <number1> <number2> - displays an interface between two patches
#
# (c) Kristian Rother, 2002
# Python License
#

import os,re,string
from pymol import cmd

TEMP_FILE = os.sep+'tmp'+os.sep+'nsc.ent'

class nscCache:
    def __init__(self):
        self.cache=[]
        self.cachedSSNames=[]
        self.cachedName=""

    def loadNSC(self,name):
        if self.cachedName==name:
            return 1
        else:
            name=name
            if os.access(name,os.F_OK):
                # load protein file
                f=open(name,'r')
                fl=f.readlines()
                f.close()
                self.cache=[]

                # parse atom lines
                actSecNum=0
                for l in fl:
                    if re.search("^ATOM",l,re.IGNORECASE):
                        secNum=int(l[82:87])
                        if actSecNum!=secNum:
                            while len(self.cache)<actSecNum:
                                self.cache.append([])
                        self.cache[actSecNum].append(l)
                                                
                self.cachedName=Name
                return 1

            print name
            return 0


    def lineContactToSSE(self,line,sse,cutoff=2.8):
        if sse=="all":
            return 1
        # parse contact map from the line
        s=line[90:-1]
        s=re.sub("\A +","",s)
        s=re.sub(" +"," ",s)
        t=string.split(s," ")
        iz=0
        # loop through all contacts
        while (iz<len(t)-1):
            contactSSE=t[iz]
            contactDist=float(t[iz+1])
            if contactSSE==sse and contactDist<=cutoff:
                return 1
            iz += 2
        return 0
       

nsc=nscCache()


def nsc(name,sekstr1="all",sekstr2="all",cutoff="2.8"):
    if nsc.loadNSC(name):        
        out=[]
        co=float(cutoff)
        if sekstr1=='all' and sekstr2=='all':
            # convert entire sco file
            for sse in nsc.cache:
                for l in sse:
                    out.append(l[:90]+'\n')
                    # out.append(sc.scoLine2Pdb(l))    
        elif len(nsc.cache)>=int(sekstr1):
            sse=nsc.cache[int(sekstr1)]
            for l in sse:
                if nsc.lineContactToSSE(l,sekstr2,co):
                    out.append(l[:90]+'\n')
                    # out.append(sc.scoLine2Pdb(l))                            

        out.append("END\n")

        o=open(TEMP_FILE,'w')
        o.writelines(out)
        o.close()
        objname=name+"_"+sekstr1+"_"+sekstr2
        objname=re.sub("\_all","",objname)
        cmd.load(TEMP_FILE,objname)
        os.remove(TEMP_FILE)
        
    else:
        print "ERROR: no NSC file found for "+name
    
cmd.extend('nsc',nsc)

