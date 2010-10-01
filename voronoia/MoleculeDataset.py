
import os,re,string
from Voronoia import *
from VolParser import *

def __init__(self):
    pass

class MoleculeDataset:

    def __init__(self,input_path,out_path,options):
        """
        Creates a MoleculeDataset object.
        Looks for files in the input dir automatically.
        """
        self.out_path = out_path        
        self.n_molecules = 0
        self.pdb_files = []
        self.vol_files = []
        
        for fn in os.listdir(input_path):
            filename = input_path + fn
            if re.search("(\.ent\Z)|(\.pdb\Z)",fn):
                self.pdb_files.append(filename)
            if re.search("(\.vol\Z)",fn):
                self.vol_files.append(filename)
                
    def get_molecule(self,filename,options):
        # create MoleculePacking object
        molecule_pack = VolParser(filename)
        molecule_pack.parse_vol_file(options)
        return molecule_pack

    def calculate_packing(self, options):
        from Voronoia import calculate_packing_from_pdb
        message = '%i packing files were calculated:\n'
        i_vol = 0
        for pdb_file in self.pdb_files:
            vol_file = calculate_packing_from_pdb(pdb_file,self.out_path,options)
            if vol_file:
                message += '%s\tfrom\t%s\n'%(vol_file,pdb_file) 
                i_vol += 1
        message = message%i_vol
        return message, None

    def report_packing(self,options):
        """
        Writes packing report files for each molecule in the
        dataset. Calculates packing .vol files on demand.
        """
        message = "Report files written for dataset:\n"
        for filename in self.vol_files:
            molecule_pack = self.get_molecule(filename,options)
            if not molecule_pack: continue

            reports = molecule_pack.report(options)
            report_file = filename[:-4]+'.report'
            if options['reportmode'] == 'text':
                open(report_file+'.txt','w').writelines(reports[0])
                message += "\t%s\n"%(report_file+'.txt')
            elif options['reportmode'] == 'html':
                open(report_file+'.html','w').writelines(reports[0])
                message += "\t%s\n"%(report_file+'.html')
            elif options['reportmode'] == 'both':
                open(report_file+'.txt','w').writelines(reports[0])
                message += "\t%s\n"%(report_file+'.txt')
                open(report_file+'.html','w').writelines(reports[1])
                message += "\t%s\n"%(report_file+'.html')
                    
        return message,message

    def compare_to_reference(self,average_pack,options):
        """
        Returns a string report of z-score rms values
        for a dataset. Calculates packing .vol files on demand.
        """
        report = Report()
        report.add_html_body("""<h2>Packing of dataset compared to reference values</h2>\n
<p>Reference file: %s</p>        
<table><tr><th>packing file</th><th>z-score-rms</th></tr>\n"""%options['reference_file'])

        zscore = ZScore(average_pack)
        for filename in self.vol_files:
            molecule_pack = self.get_molecule(filename,options)
            if not molecule_pack: continue
            zscore.compare_molecule(molecule_pack,options)
            report.add_html_body("<tr><td>%s</td><td>%4.2f</td></tr>\n"%(filename,zscore.z_score_rms))
        report.add_html_body("</table>")
        if options['reportmode'] == 'html':
            return [report.get_html(),'']
        elif options['reportmode'] == 'text':
            return [report.get_text(),'']
        elif options['reportmode'] == 'both':
            return [report.get_text(),report.get_html()]

    def calculate_average_packing(self,options):
        """
        Returns a AveragePacking object for the data in the
        entire dataset. Calculates packing .vol files on demand.
        """
        avg_pack = AveragePacking(options)
        for filename in self.vol_files:
            molecule_pack = self.get_molecule(filename,options)
            if not molecule_pack: continue
            avg_pack.count_molecule(molecule_pack)
        avg_pack.calc_average_packing()
        return avg_pack.report()


    def write_cavities(self,options):
        """
        Writes cavity center point PDB files for each
        molecule in the dataset. Calculates packing .vol
        files on demand.
        """
        message = "Cavity PDB files written for dataset:\n"
        for filename in self.vol_files:
            molecule_pack = self.get_molecule(filename,options)
            if not molecule_pack: continue
            
            cavpdb = molecule_pack.get_cavities()            
            cf = filename[:-4]+'.cav.ent'
            open(cf,'w').write(cavpdb)     
            message += "\t%s\n"%(cf)
        return [message,None]

    def write_cavity_neighbors(self,options):
        """
        Writes cavity neighbor PDB files for each
        molecule in the dataset. Calculates packing .vol
        files on demand.
        """
        message = "Cavity neighbor PDB files written for dataset:\n"
        for filename in self.vol_files:
            molecule_pack = self.get_molecule(filename,options)
            if not molecule_pack: continue
            
            cavpdb = molecule_pack.get_cavity_neighbors(options)
            cf = filename[:-4]+'.cavnb.ent'
            open(cf,'w').write(cavpdb)     
            message += "\t%s\n"%(cf)
        return [message,None]
            
