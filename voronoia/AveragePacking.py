
import re,string,math,os
from Report import Report

def __init__(self):
    pass
    
class AveragePacking:

    def __init__(self,options):
        """
        There are two atom typing schemes,
           'protor' using 19 concise types (default)
           'native' using all 173 types occuring in proteins.
        There are also two modes.
           'packing_density' (default) and 'total_volume'.
        """
        # calculating average packing
        self.atomtyping=options['atomtyping']
        self.atomtype_file=options['atomtype_file']
        self.mode = options['mode']

        self.packing_data = {}        
        self.residue_data = {}

        # handling of reference files
        self.ref_data = {}

        # lookup table for ProtOr atom types
        self.assign = {}
        for l in open(self.atomtype_file).readlines():
            t=string.split(l[:-1],'\t')
            if len(t)>=2: self.assign[t[0]]=t[1]

        # read average volumes/packing densities
        self.ref_file = options['reference_file']
        if self.ref_file and os.access(self.ref_file,os.F_OK):
            self.parse_avg_file(self.ref_file)

        self.options = options


    def parse_avg_file(self,filename):
        """
        Parses an .avg file containing average packing values.
        """
        self.ref_data = {}
        for l in open(filename,'r').readlines()[1:]:
            if not re.search("\AName",l) and not re.search("\Aatomtype",l):
                tok=string.split(l,"\t")
                if len(tok)>=4:
                    # tok[0]:atomtype
                    # tok[1]:n, tok[2]:average, tok[3]:stddev
                    self.ref_data[tok[0]] = (int(tok[1]),float(tok[2]),float(tok[3]))


    def report(self):
        """
        Returns a detailed report on the packing of this structure.
        """
        report = Report()
        report.add_html_body(self.html_report())
        if self.options['reportmode'] == 'html':
            return [report.get_html(),'']
        elif self.options['reportmode'] == 'text':
            return [report.get_text(),'']
        elif self.options['reportmode'] == 'both':
            return [report.get_text(),report.get_html()]
        else:
            raise "unknown report mode option '%s'"%options['reportmode']

    def html_report(self):
        """
        Returns a string with a atom-type-wise tabular report of
        average packing data.
        """
        # print atom info for all proteins
        out="""<table>
<tr><th>atomtype</th><th># atoms</th><th>average %s</th><th>std.dev</th><th>std.dev [%%]</th></tr>\n"""%(self.options['mode'])
        keys = self.ref_data.keys()
        keys.sort()
        for k in keys:
            rec=self.ref_data[k]
            if rec[1]>0:
                out += "<tr><td>%s</td><td>%5i</td><td>%6.3f</td><td>%6.3f</td><td>%6.3f</td></tr>\n"\
                       %(k,rec[0],rec[1],rec[2],rec[2]*100.0/rec[1])
            else:
                out += "<tr><td>%s</td><td>-</td><td>-</td><td>-</td><td>-</td></tr>\n"%(k)
        out += "</table>\n"
        return out
    

    def clear_counter(self):
        """
        Clears the memory of stored packing values.
        """
        self.packing_data = {}
        self.residue_data = {}


    def count_molecule(self,molecule_pack):
        """
        Counts packing data in the MoleculePacking object given.
        """
        for i in range(len(molecule_pack.atoms)):
            if not molecule_pack.valid[i]: continue
            atom = molecule_pack.get_atom(i)
            #print
            #for k in atom.keys(): print k,atom[k]
            
            # assign an atomtype to each atom
            atomtype = self.get_atomtype(atom)

            # get a packing density or volume value
            if self.mode in ['packing_density','total_volume']:
                value = atom[self.mode]
            else:
                value = 0.0

            # store the value
            if not self.packing_data.has_key(atomtype): self.packing_data[atomtype] = []
            self.packing_data[atomtype].append(value)



    def calc_average_packing(self):
        """
        Calculates average packing values from all
        molecules counted so far. 
        """
        self.ref_data = {}

        keys = self.packing_data.keys()
        for atomtype in keys:
            data = self.packing_data[atomtype]
            
            # calculate average
            n = len(data)
            total = 0.0
            for d in data: total += d
            average = total/n

            # calculate standard deviation
            total = 0.0
            for d in data: total += ((d-average)**2)
            if n>1:
                stddev = math.sqrt(total/(n-1))
            elif n==1:
                stddev = 0
            else:
                stddev = -1

            self.ref_data[atomtype] = (n,average,stddev)
            
        return self


    def get_atomtype(self,atom):
        """
        In the default atomtyping mode, residue+atom names are combined, e.g. TYR_CA.
        In the 'protor' atomtyping, 19 short types, like C4H3u are used.
        
        In a lookup table, the ProtOr atomtypes are stored.
        For each regular amino acid atom, the according of
        19 concise types is looked up.
        All other atom types are returned unchanged.
        """
        atomtype = atom['residue_type']+'_'+atom['atom_type']
        if self.atomtyping == 'protor':
            if self.assign.has_key(atomtype):
                atomtype = self.assign[atomtype]
        return atomtype


    def get_avg_value(self,atom):
        """Returns a tuple of (average,stddev) for the given atom)."""
        atomtype = self.get_atomtype(atom)
        if self.ref_data.has_key(atomtype):
            return (self.ref_data[atomtype][1],self.ref_data[atomtype][2])
        return None


    def count_residue(self,atom):
        code = atom['chain_id']+'_%i'%(atom['residue_number'])+'_'+atom['residue_type']
        if not self.residue_data.has_key(code): self.residue_data[code] = []
        self.residue_data[code].append(atom)
        
    def get_residue_results(self):
        results = []
        keys = self.residue_data.keys()
        keys.sort()        
        for code in keys:
            record = self.sum_atom_set(self.residue_data[code])
            t = string.split(code,'_')
            results.append((t[0],t[1],t[2],record[0],record[1],record[2],record[3],record[4],record[5]))

        return results
            
                
    def count_subset(self,molecule_pack,subset):
        """
        Returns tuple of n_atoms_total\tn_atoms\t<VdW>\t<SolvExc>\t<totVol>\tSD(<totVol>)\t<PD>\tSD(<PD>) for a defined subset.
        """
        
        set = []
        for i in range(len(molecule_pack.atoms)):
            atom = molecule_pack.get_atom(i)
            valid = 0
            if subset == "total atoms in structure": valid = 1
            elif molecule_pack.valid[i]:
                if atom['hetero']:
                    if subset == "hetero atoms": valid = 1
                    elif subset == "water molecules" and atom['residue_type']=="HOH": valid = 1
                    elif subset == "other hetero atoms" and not atom['residue_type']=="HOH": valid = 1
                else:
                    if subset == "peptidic atoms": valid = 1
                    if atom['surface']:
                        if subset == "surface": valid = 1
                        elif subset == "surface cavity neighbors" and atom['cavity_nb']: valid = 1
                        elif subset == "other surface atoms" and not atom['cavity_nb']: valid = 1
                        
                    else:
                        if subset == "buried": valid = 1     
                        elif subset == "buried cavity neighbors" and atom['cavity_nb']: valid = 1
                        elif subset == "other buried atoms" and not atom['cavity_nb']: valid = 1
            
            if valid:
                set.append(atom)

        return self.count_atom_set(set)


    def count_atom_set(self,set):
        n_atoms       = len(set)
        sum_vdw       = 0.0
        sum_solvexc   = 0.0
        sum_totvol    = 0.0
        sum_pd        = 0.0
        rtotvol = []
        rpd     = []

        for atom in set:
            sum_vdw += atom['vdw_volume']
            sum_solvexc += atom['solv_ex_volume']
            sum_totvol += atom['total_volume']
            sum_pd += atom['packing_density']
            rtotvol.append( atom['total_volume'])
            rpd.append(atom['packing_density'])

        if n_atoms>0:
            vdw = sum_vdw/n_atoms
            solvexc = sum_solvexc/n_atoms
            totvol = sum_totvol/n_atoms
            pd = sum_pd/n_atoms
        else:
            vdw = 0
            solvexc = 0
            totvol = 0
            pd = 0

        # calculate standard deviations
        sum_tot = 0.0
        sum_pd = 0.0
        for i in range(len(rtotvol)):
            sum_tot = (rtotvol[i]-totvol)**2
            sum_pd  = (rpd[i]-pd)**2
            
        if n_atoms>1:
            sd_vol = math.sqrt(sum_tot/(n_atoms-1))
            sd_pd = math.sqrt(sum_pd/(n_atoms-1))
        else:
            sd_vol = 0.0
            sd_pd  = 0.0
            
        record = [n_atoms,vdw,solvexc,totvol,sd_vol,pd,sd_pd]
        return record

    def sum_atom_set(self,set):
        n_atoms       = len(set)
        sum_vdw       = 0.0
        sum_solvexc   = 0.0
        sum_totvol    = 0.0
        sum_pd        = 0.0
        rtotvol = []
        rpd     = []

        for atom in set:
            sum_vdw += atom['vdw_volume']
            sum_solvexc += atom['solv_ex_volume']
            sum_totvol += atom['total_volume']
            sum_pd += atom['packing_density']
            rtotvol.append( atom['total_volume'])
            rpd.append(atom['packing_density'])

        if n_atoms>0: pd = sum_pd/n_atoms
        else: pd = 0

        # calculate standard deviations
        sum_pd = 0.0
        for i in range(len(rtotvol)):
            sum_pd  = (rpd[i]-pd)**2
            
        if n_atoms>1:
            sd_pd = math.sqrt(sum_pd/(n_atoms-1))
        else:
            sd_pd  = 0.0
            
        record = [n_atoms,sum_vdw,sum_solvexc,sum_totvol,pd,sd_pd]
        return record




