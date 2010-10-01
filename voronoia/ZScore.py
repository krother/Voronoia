
# ZScore.py
# 
# calculates deviations from average volumes for one protein
#
# calculates how well the packing of buried atoms in a 
# MolecularPacking object fits to reference values
# from an AveragePacking object.
#
# 0         : perfect
# 0.9 - 1.2 : nominal
# 1.5+      : weird data
# 3.0+      : error
#

import re,string,math
from Report import Report

def __init__(self):
    pass

class ZScore:

    def __init__(self,average_pack):
        self.reference = average_pack
        self.molecule  = None

        self.atom_z_scores    = []
        self.residue_z_scores = []
        self.z_score_rms = 0.0

    def get_z_score(self,atom,options):
        zscore = -1
        packing = atom[options['mode']]
        ref = self.reference.get_avg_value(atom)
        if not ref:
            # if not atom['hetero']:
            #    print "atomtype %s not in reference set"%(self.reference.get_atomtype(atom))
            return None
        else:
            ref_packing = ref[0]
            ref_deviation = ref[1]
            
            # normalize
            if ref_deviation > 0: return ((packing-ref_packing)/ref_deviation)**2
            else: return None


    def compare_molecule(self,molecule_pack,options):
        """
        Calculates the z-score (normalized difference between
        packing in a sample protein and a reference) for each
        residue in the VolParser object. Also calculates
        a z-score-rms for the entire molecule. Returns a
        ZScore object.
        """
        self.molecule = molecule_pack
        self.atom_z_scores    = []
        self.residue_z_scores = []
        self.z_score_rms = 0.0

        total_sum = 0.0
        n_atoms   = 0

        # loop over all atoms of the molecule
        i = 0        
        while i<len(self.molecule.atoms):            

            atom = self.molecule.get_atom(i)
            i += 1
            if not self.molecule.valid[i-1]:continue                        
            atomtype = atom['residue_type']+'_'+atom['atom_type']
            atomtype = re.sub(" ","",atomtype)

            zscore = self.get_z_score(atom,options)
            if zscore:
                total_sum += zscore
                n_atoms += 1


        # calculate z-score-rms
        if n_atoms > 0:
            self.z_score_rms = math.sqrt(total_sum / n_atoms)
        else:
            self.z_score_rms = -1

        return self.z_score_rms


    def report(self,options):
        report = Report()
        report.add_html_body(self.html_report())
        if options['reportmode'] == 'html':
            return [report.get_html()]
        elif options['reportmode'] == 'text':
            return [report.get_text()]
        elif options['reportmode'] == 'both':
            return [report.get_text(),report.get_html()]
        else:
            raise "unknown report mode option '%s'"%options['reportmode']

    def html_report(self):
        """
        Returns a string containing a report on the
        per-residue z-score plus the z-score-rms.
        """
        out = """
<h2>Comparison with reference set</h2>
<p>Molecule file: %s<br>
Reference file  : %s</p>\n"""%(self.molecule.vol_file,self.reference.ref_file)
        out += "<p>z-score-rms: %s</p>\n"%(self.z_score_rms)
        return out
    


