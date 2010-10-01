# -*- coding: cp1252 -*-

#
#
# a class for handling vol files.
#
#
import os,sys,re,string,math
from AveragePacking import AveragePacking
from ZScore import ZScore
from Report import Report

def __init__(self):
    pass

class VolParser:
    """
    The .vol format is very similar to the PDB format.
    It contains three additional types of data.
    1) In the coordinate lines, atomic volumes (VdW and Solvent
       excluded volume) and surface (binary) are listet.
    2) Lines starting with HOLE NUMBER list atom numbers participating
       in a particular cavity.
    3) Lines starting with SELECT will determine only parts of the
       file for calculation. The selections are cumulative. Examples;
       SELECT chain A AND resi 1-100
       SELECT name CA AND resn TYR
       SELECT resn HOH
       SELECT resi 65
       when no SELECT commands are specified, all atoms are selected
       automatically.
    """
    def __init__(self,vol_file):
        self.vol_file = vol_file
        self.atoms = []
        self.atom_dict = []
        self.cavities = []
        self.cav_dist = []
        self.cav_dict = {}        

    def parse_atom(self,l,options):
        hetero = 0
        if re.search("\AHETATM",l):
            hetero = 1

        # parse atom number
        atom_number    = int(l[6:11])
        # parse atom type, residue,chain
        atom_type      = re.sub('\s','',l[13:16])
        residue_type   = re.sub('\s','',l[17:20])
        chain_id       = l[21]
        residue_number = int(l[23:26])
        pdb_line       = l[:66]
                        
        # parse atomic volumes
        vdw_volume        = float(l[66:73])
        solv_ex_volume    = float(l[73:80])
        total_volume      = vdw_volume + solv_ex_volume
        if total_volume > 0:
            packing_density   = vdw_volume / total_volume
        else:
            packing_density = 0.0
        surface           = 1-int(l[80:-1])

        # check for cavity neighbors
        if self.cav_dict.has_key(str(atom_number)):
            cavity_nb = 1
            cavities = self.cav_dict[str(atom_number)]
        else:                
            cavities  = []
            cavity_nb = 0
            
        hetero_nb = 0

        # coordinates
        coord = [float(l[30:38]),float(l[38:46]),float(l[46:54])]
                

        atom = {
            'chain_id':chain_id,
            'residue_number':residue_number,
            'residue_type':residue_type,
            'atom_type':atom_type,
            'packing_density':packing_density,
            'vdw_volume':vdw_volume,
            'solv_ex_volume':solv_ex_volume,
            'total_volume':total_volume,
            'surface':surface,
            'cavity_nb':cavity_nb,
            'cavities':cavities,
            'hetero_nb':hetero_nb,
            'pdb_line':pdb_line,
            'index':atom_number,
            'hetero':hetero,
            'coord':coord,
        }

        # conditional parsing
        valid = 1
        if atom['hetero']:
            if atom['residue_type'] == 'HOH':
                if not options['include_water']: valid = 0
            elif not options['include_hetero']: valid = 0
        # if hetero_nb and options['discard_hetero_neighbors']: valid = 0
        if atom['cavity_nb'] and options['discard_cavity_neighbors']: valid = 0
        if not atom['cavity_nb'] and options['discard_non_cavity_neighbors']: valid = 0
        if atom['surface'] and options['discard_surface']: valid = 0
        if not atom['surface'] and options['discard_buried']: valid = 0
        atom['valid'] = valid
        
        return atom
    
    def parse_vol_file(self,options):

        self.atoms = []        
        self.cavities = []
        self.valid = []
        self.selects = []
        # each entry in this list is a tuple containing 
        # lots of information per atom

        i_atom    = 0
        i_cavity  = 0
        self.cav_dict  = {}
        if not os.access(self.vol_file,os.F_OK): return None
        
        # loop through all input lines
        for l in open(self.vol_file).readlines():
            # pares lines indicating atom selections
            if re.search("\ASELECT",l):
                ll = re.sub('\s+',' ',l[7:-1])                
                terms = string.split(ll,'AND')
                select = []
                for term in terms:
                    tokens = string.split(term,' ')
                    if len(tokens) == 2:
                        keyword = tokens[0]
                        value = tokens[1]
                        if keyword != 'resi':
                            select.append((keyword,value))
                        else:
                            vv = string.split(value,'-')
                            if len(vv) == 2:
                                select.append('min_resi',int(vv[0]))
                                select.append('max_resi',int(vv[1]))
                            else:
                                select.append('resi',int(vv[0]))
                                              
                self.selects.append(select)
                
            # parse lines indicating cavities
            if re.search("\AHOLE NUMBER",l):
                self.cavities.append([])
                ll = re.sub('\r','',l[12:-1])
                for a in string.split(ll,' ')[3:]:
                    if not self.cav_dict.has_key(a): self.cav_dict[a] = []
                    self.cav_dict[a].append(i_cavity)
                i_cavity += 1
            
            # parse coordinate lines
            if re.search("(\AATOM)|(\AHETATM)",l):
                atom = self.parse_atom(l,options)

                if atom['cavity_nb']:
                    cavities = self.cav_dict[str(atom['index'])]
                    for c in cavities:
                        self.cavities[c].append(i_atom)

                record = [
                    atom['chain_id'],
                    atom['residue_number'],
                    atom['residue_type'],
                    atom['atom_type'],
                    atom['packing_density'],
                    atom['vdw_volume'],
                    atom['solv_ex_volume'],
                    atom['total_volume'],
                    atom['surface'],
                    atom['cavity_nb'],
                    atom['cavities'],
                    atom['hetero_nb'],
                    atom['pdb_line'],
                    atom['index'],
                    atom['hetero'],
                    atom['coord']
                    ]
                

                self.atom_dict.append(atom)
                self.atoms.append(record)
                self.valid.append(atom['valid'])
                i_atom += 1


        # calculate medium distance of cavity neighbors to cavity center
        i_cav = 0
        i_cav_atom = 1
        for cav in self.cavities:
            i_cav += 1
            # calculate center of mass
            sumx = 0.0
            sumy = 0.0
            sumz = 0.0
            n = 0.0
            for iatom in cav:
                if iatom == None: continue
                n += 1.0
                coord = self.get_atom(iatom)['pdb_line']
                sumx += float(coord[30:38])
                sumy += float(coord[38:46])
                sumz += float(coord[46:54])
            if n > 0:
                mx = sumx/n
                my = sumy/n
                mz = sumz/n
                cavcenter = [mx,my,mz]
                # calculate medium euclidean distance to cavity centre
                eu_sum = 0.0
                for iatom in cav:
                    if iatom == None: continue
                    coord = self.get_atom(iatom)['coord']
                    eu_sum += self.distance(coord,cavcenter)
                avg_eu = eu_sum/n
                self.cav_dist.append(avg_eu)
            else:
                self.cav_dist.append(-1.0)

        #
        # parse the select statements
        # and switch some atoms off
        #
        # please note that the options like 'discard_buried' override
        # the SELECT statements
        #
        if self.selects != []:
            for i in range(len(self.atoms)):
                if not self.valid[i]:continue
                valid = 0
                for s in self.selects:
                    all_matched = 1
                    for c in s:
                        keyword = c[0]
                        value = c[1]
                        if keyword == 'resn':
                            if self.atoms[i][2] != value: all_matched = 0
                        elif keyword == 'resi':
                            if self.atoms[i][1] != value: all_matched = 0
                        elif keyword == 'min_resi':
                            if self.atoms[i][1] < value: all_matched = 0
                        elif keyword == 'max_resi':
                            if self.atoms[i][1] > value: all_matched = 0
                        elif keyword == 'name':
                            if self.atoms[i][3] != value: all_matched = 0
                        elif keyword == 'chain':
                            if self.atoms[i][0] != value: all_matched = 0
                    if all_matched: valid = 1
                self.valid[i] = valid

                        
    def report(self,options):
        """
        Returns a detailed report on the packing of this structure.
        """
        report = Report()
        report.add_html_body(self.html_report(options))
        if options['reportmode'] == 'html':
            return [report.get_html(),None]
        elif options['reportmode'] == 'text':
            return [report.get_text(),None]
        elif options['reportmode'] == 'both':
            return [report.get_text(),report.get_html()]
        else:
            raise "unknown report mode option '%s'"%options['reportmode']


    def html_report(self,options):
        report = """

<h1>Packing report</h1>
"""
        total_atoms = len(self.atoms)        
        calc_atoms  = 0
        sum_pd      = 0.0
        for i in range(len(self.atoms)):
                       if self.valid[i]:
                           calc_atoms += 1
                           sum_pd += self.get_atom(i)['packing_density']

        if calc_atoms > 0:
            avg_pd = sum_pd / calc_atoms
        else:
            avg_pd = 0
            
        n_cavi      = len(self.cavities)

        # calculate z-score-rms
        avg_pack = AveragePacking(options)
        z = ZScore(avg_pack)
        z_score = z.compare_molecule(self,options)


        report += """
<h2>1.SUMMARY</h2>

<table cellspacing=0 cellpadding=4 border=1>
<tr><td width='25%%'><b>structure file               </b></td><td width='75%%'>%s</td></tr>
<tr><td><b>reference file               </b></td><td>%s</td></tr>
<tr><td></td></tr>
<tr><td><b>total atoms in structure     </b></td><td>%5i</td></tr>
<tr><td><b>atoms used for calculation   </b></td><td>%5i</td></tr>
<tr><td><b>average packing density      </b></td><td>%5.3f</td></tr>
<tr><td><b>number of cavities           </b></td><td>%5i</td></tr>
<tr><td><b>z-score-rms                  </b></td><td>%5.2f</td></tr>
</table>
"""%(self.vol_file,options['reference_file'],total_atoms,calc_atoms,avg_pd,n_cavi,z_score)

        report += """

&nbsp;
<h2>2.ABBREVIATIONS USED</h2>

<table cellspacing=0 cellpadding=4 border=1>
<tr><td width='25%'><b>   chain         </b></td><td width='75%'> chain ID (one letter, e.g. H)</td></tr>
<tr><td><b>   resi_no       </b></td><td> residue number (integer, e.g. 131)</td></tr>
<tr><td><b>   resi_type     </b></td><td> residue type (e.g. TYR, HOH)</td></tr>
<tr><td><b>   n_atoms_total </b></td><td> number of atoms in the structure</td></tr>
<tr><td><b>   n_atoms       </b></td><td> number of atoms used in the calculation</td></tr>
<tr><td><b>   atom_no       </b></td><td> atom number in PDB and .vol file.</td></tr>
<tr><td><b>   atom_name     </b></td><td> atom name (three-letter, from PDB file, e.g. CA, CD2)</td></tr>
<tr><td><b>   atom_type     </b></td><td> classes by which average packing is calculated.<br>
This is either 'residuename_atomname' e.g. 'TYR_CA', or one of 18 concise ProtOr atom types, e.g. 'C4H3u' for methyl groups.</td></tr>
<tr><td><b>   &lt;VdW&gt;         </b></td><td> average Van-der-Waals volume in cubic Angstroms</td></tr>
<tr><td><b>   &lt;SolvExc&gt;     </b></td><td> average Solvent Excluded volume in cubic Angstroms</td></tr>
<tr><td><b>   &lt;totVol&gt;      </b></td><td> average total volume (VdW + SolvExc) in cubic Angstroms</td></tr>
<tr><td><b>   &lt;PD&gt;          </b></td><td> average packing density (VdW/totVol)</td></tr>
<tr><td><b>   SD            </b></td><td> standard deviation</td></tr>
<tr><td><b>   [&Aring;<sup>3</sup>]         </b></td><td> cubic Angstroms</td></tr>
</table>

&nbsp;
<h2>3.ATOM NUMBERS AND PACKING</h2>

<table cellspacing=0 cellpadding=4 border=1>
<tr><th width='20%'>atomic subset</th><th width='5%'>n_atoms</th><th width='10%'>&lt;VdW&gt; [&Aring;<sup>3</sup>]</th><th width='10%%'>&lt;SolvExc&gt; [&Aring;<sup>3</sup>]</th><th width='10%'>&lt;totVol&gt; [&Aring;<sup>3</sup>]</th><th width='10%'>SD(&lt;totVol&gt;) [&Aring;<sup>3</sup>]</th><th width='10%'>&lt;PD&gt;</th><th width='10%'>SD(&lt;PD&gt;)</th></tr>
"""
        subsets = ["total atoms in structure           ",
                   "   peptidic atoms                  ",
                   "      buried                       ",
                   "         buried cavity neighbors   ",
                   "         other buried atoms        ",
                   "      surface                      ",
                   "         surface cavity neighbors  ",
                   "         other surface atoms       ",
                   "   hetero atoms                    ",
                   "      water molecules              ",
                   "      other hetero atoms           "]
        avg_pack = AveragePacking(options)
        for s in subsets:            
            rec = avg_pack.count_subset(self,re.sub('(^\s+)|(\s+\Z)','',s))            
            line = '<tr align="right"><td align="left">%s</td><td>%i</td><td>%5.2f</td><td>%5.2f</td><td>%5.2f</td><td>%5.2f</td><td>%5.3f</td><td>%5.3f</td></tr>\n'%(s,rec[0],rec[1],rec[2],rec[3],rec[4],rec[5],rec[6])
            report += line
            
        report += """</table>

&nbsp; 
<h2>4.CAVITIES</h2>

<table cellspacing=0 cellpadding=4 border=1>
<tr><th width='15%'>cavity_no</th><th width='10%'>n_neighbors</th><th width='10%'>&lt;dist&gt; [&Aring;]</th><th width='10%'>%surface</th><th width='30%'>atom_types</th><th width='20%'>atom_indices</th></tr>
"""
        for i in range(len(self.cavities)):
            cav = self.cavities[i]
            nsurf = 0
            for c in cav:
                if self.get_atom(c)['surface']: nsurf +=1
            atypes = string.join([avg_pack.get_atomtype(self.get_atom(ai)) for ai in cav],', ') 
            anums = string.join([str(num) for num in cav],', ')
            eu_med = self.cav_dist[i]
            report += '<tr align="right"><td>%3i</td><td>%3i</td><td>%5.2f</td><td>%4.1f</td><td align="left">%s</td><td align="left">%s</td></tr>\n'%(i,len(cav),eu_med,nsurf*100.0/len(cav),atypes,anums)
            
        
        report += """</table>

&nbsp;
<h2>5.PACKING PER RESIDUE</h2>

<table cellspacing=0 cellpadding=4 border=1>
<tr><th width='5%'>chain</th><th width='10%'>resi_no</th><th width='10%'>resi_type</th><th width='5%'>n_atoms</th><th width='10%'>sum(VdW) [&Aring;<sup>3</sup>]</th><th width='10%'>sum(SolvExc) [&Aring;<sup>3</sup>]</th><th width='10%'>sum(totVol) [&Aring;<sup>3</sup>]</th><th width='10%'>&lt;PD&gt;</th><th width='10%'>SD(&lt;PD&gt;)</th></tr>
"""
        avg_pack = AveragePacking(options)
        for i in range(len(self.atoms)):
            if self.valid[i]:
                atom = self.get_atom(i)
                avg_pack.count_residue(atom)        
        for r in avg_pack.get_residue_results():
                report += '<tr align="right"><td>%s</td><td>%s</td><td>%s</td><td>%i</td><td>%5.2f</td><td>%5.2f</td><td>%5.2f</td><td>%5.2f</td><td>%5.3f</td></tr>\n'%r

        report += """</table>

&nbsp;  
<h2>6.PACKING PER ATOM</h2>

<table cellspacing=0 cellpadding=4 border=1>
<tr><th width='10%'>atom_no</th><th width='5%'>chain</th><th width='10%'>resi_no</th><th width='10%'>residue_name</th><th width='10%'>atom_name</th><th width='10%'>atom_type</th><th width='10%'>VdWVol [&Aring;<sup>3</sup>]</th><th width='10%'>SolvExcVol [&Aring;<sup>3</sup>]</th><th width='10%'>totVol [&Aring;<sup>3</sup>]</th><th width='10%'>PD</th></tr>
"""
        for i in range(len(self.atoms)):
            if self.valid[i]:
                atom = self.get_atom(i)
                atomtype = avg_pack.get_atomtype(atom)
                report += '<tr align="right"><td>%i</td><td>%s</td><td>%i</td><td>%s</td><td>%s</td><td>%s</td><td>%5.2f</td><td>%5.2f</td><td>%5.2f</td><td>%5.3f</td></tr>\n'%(atom['index'],atom['chain_id'],atom['residue_number'],atom['residue_type'],atom['atom_type'],atomtype,atom['vdw_volume'],atom['solv_ex_volume'],atom['total_volume'],atom['packing_density'])

        report += """</table>"""
        return report



    def summary(self,options):
        """
        Returns a short one-line report on the packing of this structure.
        """
        n_counted = 0
        sum_value = 0
        for i in range(len(self.atoms)):
            atom = self.get_atom(i)
            if atom['surface'] and options['discard_surface']: continue
            if not atom['surface'] and options['discard_buried']: continue
            if atom['hetero_nb'] and options['discard_hetero_neighbors']: continue
            if atom['cavity_nb'] and options['discard_cavity_neighbors']: continue
            if options['mode'] == 'packing_density': sum_value += atom['packing_density']
            else: sum_value += atom['total_volume']
            n_counted += 1
            
        avg_value = sum_value*1.0/n_counted

        n_cavities = len(self.cavities)
        avg_cavnb = 0.0
        for cav in self.cavities:
            avg_cavnb += len(cav)
        if n_cavities > 0:
            avg_cavnb = avg_cavnb / (n_cavities*1.0)
            
        return "%s\t%i tot.atoms\t%i count.atoms\t%5.2f avg.%s\t%i cavities\t%4.1f avg.cav.neighbors"%(self.vol_file,len(self.atoms),n_counted,avg_value,options['mode'],n_cavities,avg_cavnb)

    def write_pdb_file(self,filename,options):
        """
        Writes a pdb file containing all valid atoms. The bfactor option
        may be set to 'packdens' - substitutes bfactor with packing density.
        'volume' - same with total atomic volume
        'minuspackdens' - same with 1-packing density
        'zscore' - z-score (normalized difference to reference value)
        """
        bfactor = options['bfactor']
        out = []
        zsc = ZScore(AveragePacking(options))
        for i in range(len(self.atoms)):
            if not self.valid[i]:continue
            atom = self.get_atom(i)
            pl = atom['pdb_line']

            # now get the value to be inserted as bfactor
            #
            # For some programs they need to be normalized
            #
            # Hallo Fiete!
            #
            # unter 'packdens' kannst Du dich austoben.
            # Habe ein Beispiel reingeschrieben,
            # in Python sollte etwas Mathe nicht so schwer sein.
            # Aber frag mich trotzdem falls Du im Code versumpfst.
            # 11.7.06 KR
            #
            if bfactor == 'packdens':
                value = atom['packing_density']
                # normalize values from 0.4 to 0.9 to bfactors between 0.0 and 30.
                # 'packdens'
                vmin = 0.4
                vmax = 0.9
                bmin = 0.0
                bmax = 30.0
                value = ((value-vmin)/vmax) * (bmax-bmin) + bmin
                # 
            if bfactor == 'minuspackdens':
                value = 1-atom['packing_density']
                # normalize values from 0.4 to 0.9 to bfactors between 30.0 and 0.0
                # 'minuspackdens'
                vmin = 0.4
                vmax = 0.9
                bmin = 0.0
                bmax = 30.0
                value = ((value-vmin)/vmax) * (bmax-bmin) + bmin
                # 
            elif bfactor == 'volume':
                value = atom['total_volume']
                # nomalize
                # (insert code on demand)
            elif bfactor == 'zscore':
                value = zsc.get_z_score(atom,options)
                if not value: value=100.0
                # normalize
                # (insert code on demand)

            pl = pl[:60]+"%6.2f"%(value)+pl[66:]
            out.append(pl+'\n')
        open(filename,'w').writelines(out)

    def distance(self,coord1,coord2):
        x = (coord1[0]-coord2[0])**2
        y = (coord1[1]-coord2[1])**2
        z = (coord1[2]-coord2[2])**2
        dist = math.sqrt(x+y+z)
        return dist
    
    def get_surfdist(self):
        """
        Calculates a list of distances to the next surface atom.
        """
        surfdist = []
        for i in range(len(self.atoms)):
            atom = self.get_atom(i)
            if not self.valid[i]:
                surfdist.append(-1)
            elif atom['surface']:
                surfdist.append(0)
            else:
                min_dist = 9999.0
                for j in range(len(self.atoms)):
                    atom2 = self.get_atom(j)
                    if atom2['surface']:
                        # calculate distance
                        dist = self.distance(atom['coord'],atom2['coord'])
                        min_dist = (dist<min_dist) and dist or min_dist
                surfdist.append(min_dist)
        return surfdist

    def get_surfdist_report(self,surfdist,options):
        bins = options['surfdist_bins']
        maxdepth = options['surfdist_maxdepth']
        step = maxdepth/bins
        count = [0] * bins
        bigger = 0
        total_count = 0
        for s in surfdist:
            if s == -1 or s == 9999: continue
            for i in range(bins):
                if s<=(i+1)*step:
                    total_count += 1
                    count[i] += 1
                    break
            if s>maxdepth: bigger += 1

        # write table
        report = """<h2>Distances to next surface atom</h2>

<table>
<tr><th>bin</th><th>dist</th><th>#atoms</th></tr>
"""
        for i in range(bins):
            report += "<tr><td>%i</td><td>%5.2f</td><td>%i</td></tr>\n"%(i+1,step*(i+1),count[i])

        report += "</table><p>total atoms counted: %i<br>"%total_count
        report += "\natoms out of range: %i</p>\n"%bigger

        hreport = Report()
        hreport.add_html_body(report)
        if options['reportmode'] == 'html':
            return [hreport.get_html(),None]
        elif options['reportmode'] == 'text':
            return [hreport.get_text(),None]
        elif options['reportmode'] == 'both':
            return [hreport.get_text(),hreport.get_html()]
        else:
            raise "unknown report mode option '%s'"%options['reportmode']
        return hreport
                
        
    def get_cavities(self):
        """
        Returns a string containing a PDB file of all cavity
        centers.
        The PDB file contains one artificial atom
        for the center point of each cavity in a molecule.
        """
        pdb = ""
        i_cav = 0
        i_cav_atom = 1
        for cav in self.cavities:
            i_cav += 1
            # calculate center of mass
            sumx = 0.0
            sumy = 0.0
            sumz = 0.0
            n = 0.0
            for iatom in cav:
                if iatom == None: continue
                n += 1.0
                coord = self.get_atom(iatom)['pdb_line']
                sumx += float(coord[30:38])
                sumy += float(coord[38:46])
                sumz += float(coord[46:54])
            if n > 0:
                mx = sumx/n
                my = sumy/n
                mz = sumz/n
                # calculate medium euclidean distance to cavity centre
                eu_sum = 0.0
                for iatom in cav:
                    if iatom == None: continue
                    coord = self.get_atom(iatom)['pdb_line']
                    eu_sum += math.sqrt((mx-float(coord[30:38]))**2 + (my-float(coord[38:46]))**2 + (mz-float(coord[46:54]))**2) 
                avg_eu = eu_sum/n
                
                # create artificial PDB atom entry
                atomline = "ATOM  %5i  O   %3i  %4i     %7.3f %7.3f %7.3f  1.00 %5.2f   CAV\n"%(i_cav_atom,i_cav,i_cav,mx,my,mz,avg_eu)
                pdb += atomline
                i_cav_atom += 1

        pdb += "END\n"
        return pdb


    def get_cavity_neighbors(self,cut_off=0.0):
        """
        Returns a string containing a PDB file with neighbor
        atoms to all cavities.
        """
        pdb = ""
        atom_no = 1
        i_cav = 0
        for cav in self.cavities:
            i_cav += 1
            # get all neighbors
            for i_atom in cav:
                if i_atom:
                    coord = self.get_atom(i_atom)['pdb_line']
                    # relabel the residue number as cavity no.
                    atomline = "ATOM  %5i"%(atom_no)+coord[11:23]+"%3i"%(i_cav)+coord[26:]+'\n'
                    pdb += atomline
                    atom_no += 1
                    
        pdb += "END\n"
        return pdb
        

    def get_atom(self,index):
        """
        Returns a dictionary with all properties of the indexed atom
        """
        if index<0 or index >= len(self.atoms): return None
        return self.atom_dict[index]

    def replace_bfactors(self,infile,outfile,options,check_valid=False):
        """
        Writes the .vol file 'infile' to 'outfile', but replaces the bfactor
        column with either packing density or volume.
        Also adds a short remark before the first atom.
        """
        extralines_pd="""REMARK 
REMARK   COLUMNS MODIFIED FROM ORIGINAL PDB FILE:
REMARK      PD     - LOCAL PACKING DENSITY
REMARK      VDWVOL - VOLUME INSIDE VAN-DER-WAALS SPHERE
REMARK      SEVOL  - VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
REMARK      BUR    - INDICATES BURIED ATOMS (1=BURIED, 0=SURFACE)
REMARK  
REMARK  AT# ELEM RES C RES#       X       Y       Z     PD   PD    VDWVOL  SEVOL BUR
"""
        extralines_vol="""REMARK 
REMARK   COLUMNS MODIFIED FROM ORIGINAL PDB FILE:
REMARK      VOL    - TOTAL ATOMIC VOLUME
REMARK      VDWVOL - VOLUME INSIDE VAN-DER-WAALS SPHERE
REMARK      SEVOL  - VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
REMARK      BUR    - INDICATES BURIED ATOMS (1=BURIED, 0=SURFACE)
REMARK  
REMARK  AT# ELEM RES C RES#       X       Y       Z     VOL  VOL   VDWVOL  SEVOL BUR
"""

        extralines_z="""REMARK 
REMARK   COLUMNS MODIFIED FROM ORIGINAL PDB FILE:
REMARK      Z      - Z-SCORE
REMARK      VDWVOL - VOLUME INSIDE VAN-DER-WAALS SPHERE
REMARK      SEVOL  - VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
REMARK      BUR    - INDICATES BURIED ATOMS (1=BURIED, 0=SURFACE)
REMARK  
REMARK  AT# ELEM RES C RES#       X       Y       Z     Z    Z     VDWVOL  SEVOL BUR
"""

        f_atom = 0
        out = ""
        mode = options['bfactor']
        
        if mode == 'zscore':
            z = ZScore(AveragePacking(options))
            
        for l in open(infile).readlines():
            if re.search('(^ATOM)|(^HETATM)',l):
                if not f_atom:
                    f_atom = 1
                    if mode == 'packing_density': out += extralines_pd
                    elif mode == 'volume': out += extralines_vol
                    elif mode == 'zscore': out += extralines_z
                # extract volume data
                atom = self.parse_atom(l,options)
                if check_valid and not atom['valid']: continue
                if mode == 'packing_density': value =  atom['packing_density']
                elif mode == 'zscore':
                    value = z.get_z_score(atom,options)
                    if not value: value = -1
                else: value = atom['total_volume']
                newline = l[:55] + "%5.2f %5.2f"%(value,value) + l[67:]
                out += newline
            else:
                out += l

        open(outfile,'w').writelines(out+'\n')
        
    def my_file(self):
        """
        Code von Fiete fuer die Website. 3/2007
        """
        pdb = ""
        das = "ging schief"
        atomno = 0  
        pack_dens = 0 
        if not os.access(self.vol_file,os.F_OK): return das
        for l in open(self.vol_file).readlines():
            if re.search("(\AATOM)|(\AHETATM)",l):
                atomno += 1
                # new
                pack_dens = self.get_atom(atomno-1)['packing_density']
                vdw_volume        = float(l[67:73])
                solv_ex_volume    = float(l[74:80])
                total_volume      = vdw_volume + solv_ex_volume
                #pack_dens  = vdw_volume / total_volume 
                madeline = l[0:60]
                madeline += "%5.2f "%(pack_dens)
                madeline += l[67:]
                #end new
                pdb += madeline  
            if re.search("\AHOLE",l): 
                pdb += l
            if re.search("\AEND",l):
                i_cav = 0
                for cav in self.cavities:
                    i_cav += 1
                    # calculate center of mass
                    sumx = 0.0
                    sumy = 0.0
                    sumz = 0.0
                    n = 0.0
                    for iatom in cav:
                        if iatom == None: continue
                        n += 1.0
                        coord = self.get_atom(iatom)['pdb_line']
                        sumx += float(coord[30:38])
                        sumy += float(coord[38:46])
                        sumz += float(coord[46:54])
                    if n > 0:
                        mx = sumx/n
                        my = sumy/n
                        mz = sumz/n
                        # calculate medium euclidean distance to cavity centre
                        eu_sum = 0.0
                        atomno += 1
                        for iatom in cav:
                            if iatom == None: continue
                            coord = self.get_atom(iatom)['pdb_line']
                            eu_sum += math.sqrt((mx-float(coord[30:38]))**2 + (my-float(coord[38:46]))**2 + (mz-float(coord[46:54]))**2) 
                        avg_eu = eu_sum/n 
                        # create artificial PDB atom entry
                        atomline = "ATOM  %5i  O   HOH  %4i     %7.3f %7.3f %7.3f  1.00 0 %5.2f   CAV\n"%(atomno,i_cav,mx,my,mz,eu_sum)
                        pdb += atomline 
             
#             pdb += l             
        return pdb


        
                
