#! /usr/bin/python

# Voronoia Copyright Notice
# ============================
#
# The Voronoia source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# Voronoia is Copyright (C) 2007 by Andrean Goede and Kristian Rother
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of the authors not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# A. GOEDE AND K. ROTHER DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL THEY BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

import Voronoia_GUI
from VolParser import *
from AveragePacking import *
from ZScore import *
from MoleculeDataset import *

import os,sys,re
import tempfile

INSTALL_DIR = Voronoia_GUI.__file__[:-16]
# patch paths on Py2Exe
INSTALL_DIR = re.sub('library\.zip','',INSTALL_DIR)
# there is another strange backslash bug here:
#   each time i launch the module after modifications, the
#   INSTALL_DIR variable comes up one char too short;
#   in the second run, it works.
if INSTALL_DIR[-1] != os.sep: INSTALL_DIR += os.sep

REFERENCE_FILE = INSTALL_DIR + 'data'+os.sep+"avg_scop_protor.avg"
EXECUTABLE = 'data'+os.sep+"get_volume.exe"
ATOMTYPE_FILE = INSTALL_DIR +'data'+os.sep+"protor_atomtypes.txt"
ABOUT_FILE = INSTALL_DIR + 'data'+os.sep+'about.txt'
HELP_FILE = INSTALL_DIR + 'data'+os.sep+'README.TXT'
APPLICATION_NAME = 'Voronoia Packing Toolbox'


def __init__(self):
    pass

def get_options():
    options = {
        'grid_dist':0.1,
        'include_hetero':False,
        'include_water':False,
        'executable':EXECUTABLE,
        'executable_path':INSTALL_DIR,
        'discard_surface':True,
        'discard_buried':False,
        'discard_hetero_neighbors':False,
        'discard_cavity_neighbors':False,
        'discard_non_cavity_neighbors':False,
        'reference_file':REFERENCE_FILE,
        'atomtyping':'protor',
        'atomtype_file':ATOMTYPE_FILE,
        'radii':'protor',
        'mode':'packing_density',
        'reportmode':'text',
        'dataset':None,
        'bfactor':'packing_density',
        'surfdist_bins':20,
        'surfdist_maxdepth':15.0
        }
    return options


def calculate_packing_from_pdb(pdb_file,out_path,options):
    """
    Invokes a binary application that calculates a .vol file
    containing packing data from a PDB file.
    """

    vol_file = pdb_file[:-4]+'.vol'
    vol_file = out_path + string.split(vol_file,os.sep)[-1]
    # make temporary pdb file without water and ligands.
    temp_filename = tempfile.mktemp()+'.pdb'
    p = []
    for l in open(pdb_file).readlines():        
        if (not options['include_water']) and l[17:20] == 'HOH': continue
        elif (not options['include_hetero']) and l[:6] == 'HETATM': continue
        else: p.append(l)
    open(temp_filename,'w').writelines(p)
    log_file = tempfile.mktemp()+'.log'

    old_cwd = os.getcwd()
    if options['executable_path']: os.chdir(options['executable_path'])
    
    if os.sep == '\\':
        # on Windows
        command = """%s rad:%s ex:%s x:yes i:"%s" l:"%s" o:"%s" """%(options['executable'],options['radii'],options['grid_dist'],temp_filename,log_file,vol_file)
    else:
        # other
        command = '%s rad:%s ex:%s x:yes i:%s o:%s'%(options['executable'],options['radii'],options['grid_dist'],temp_filename,vol_file)

    print "\nExecuting packing calculation\n",command
    os.system(command)
    
    if os.access(log_file,os.F_OK): os.remove(temp_filename)
    if os.access(log_file,os.F_OK): os.remove(log_file)
    
    if os.access(vol_file,os.F_OK):    
        # read packing file
        vol_obj = VolParser(vol_file)
        vol_obj.parse_vol_file(options)

        # replace b factors in .vol file by packing data
        print options
        print vol_file
        vol_obj.replace_bfactors(vol_file,vol_file,options)

    print "Calculation finished"

    os.chdir(old_cwd)    
    return vol_file


    
def calculate_packing_from_pymol(selection,out_path,options):
    """
    Invokes a binary application that calculates a .vol file
    containing packing data from a PyMOL selection or object.
    Uses calculate_packing_from_pdb via a temporary PDB file.
    """
    # write the pymol selection to a temporary file
    sel = re.sub('[\(\)\[\]\,\#\*\?]','',selection)
    temp_file = out_path + 'pymol_sel_%s.ent'%(sel)
    i = 1
    while os.access(temp_file,os.F_OK):
        temp_file = out_path + 'pymol_sel_%s%i.ent'%(sel,i)
        i += 1
    cmd.do('save %s,%s'%(temp_file,selection))

    # calculate on the pdb file
    vol_file = calculate_packing_from_pdb(temp_file,out_path,options)
    os.remove(temp_file)

    return vol_file


def process_command(command,input_file,out_path,options):
    #
    # calculate packing files
    #
    message = None
    avg_pack = AveragePacking(options)

    if options['dataset']: 
        dataset = MoleculeDataset(options['dataset'],out_path,options)
        if command == "calc": message = dataset.calculate_packing(options)
        elif command == "report": message = dataset.report_packing(options)
        elif command == "average": message = dataset.calculate_average_packing(options)
        elif command == "compare": message = dataset.compare_to_reference(avg_pack,options)
        elif command == "cavities": message = dataset.write_cavities(options)
        elif command == "cavneighbors": message = dataset.write_cavity_neighbors(options)
        elif command == "surfdist": message = dataset.write_surfdist(options)

    elif command == "calc":
        # calculate .vol from .pdb file
        vol_file = calculate_packing_from_pdb(input_file,out_path,options)
        message = ("Packing file %s calculated from PDB file\n%s"%(vol_file,input_file),vol_file)
        
    elif input_file:
            # read packing file
            vol_obj = VolParser(input_file)
            vol_obj.parse_vol_file(options)
        
            if command == "compare":
                # compares content of a .avg file to packing data.            
                zscore = ZScore(avg_pack)
                zscore.compare_molecule(vol_obj,options)
                message = zscore.report(options)

            elif command == "report":
                # create packing reports            
                message = vol_obj.report(options)

            elif command == "average":
                # calculate average packing           
                avg_pack.count_molecule(vol_obj)
                avg_pack.calc_average_packing()
                message = avg_pack.report()

            elif command == "surfdist":
                # write table of surface distances
                surf_dist = vol_obj.get_surfdist()
                message = vol_obj.get_surfdist_report(surf_dist,options)
                

            elif command == "cavities":
                # write cavities as PDB files
                pdb = vol_obj.get_cavities()
                out_file = vol_obj.vol_file[:-4]+'.cav.ent'
                open(out_file,'w').write(pdb)
                message = ["%i cavities written\n\nfrom %s\n\nto %s"%(len(vol_obj.cavities),input_file,out_file),out_file]
        
            elif command == "cavities_php":
                # write cavities as PDB files (for web site only)
                #
                # BAH dirty code!!!
                #
                pdb = vol_obj.my_file()
                fn = string.split(vol_obj.vol_file,os.sep)[-1]
                out_file = out_path + fn[:-4]+'_cav.pdb'
                open(out_file,'w').write(pdb)
                message = ("%i cavities written\n\nfrom %s\n\nto %s"%(len(vol_obj.cavities),input_file,out_file),out_file)
                
            elif command == "cavneighbors":
                # write cavity neighbors
                pdb = vol_obj.get_cavity_neighbors()
                out_file = vol_obj.vol_file[:-4]+'.cavnb.ent'
                open(out_file,'w').write(pdb)
                message = ["Cavity neighbor atoms for %i cavities written\n\nfrom %s\n\nto %s"%(len(vol_obj.cavities),input_file,out_file),out_file]

    return message

#
#
# Command-line interface
#
#
if __name__ == '__main__':
    doc = """

    Packing Tools (c) Andrean Goede and Kristian Rother 2006

    usage

    Packing.py <command> [options] <input file>

    <command> is one of:
    calc         - calculate packing .vol file(s)
    report       - write concise packing report from .vol file(s)
    average      - calculate average packing densities from .vol file(s)
    compare      - compare .vol file(s) to average packing densities
    surfdist     - calculate distances to next surface atom 
    cavities     - write cavities from .vol file(s) as PDB file(s)
    cavneighbors - write cavity neighbor atoms from .vol file(s) as PDB file(s)

    [options] for 'calc':
    -g <float> - grid distance for 'calc' (default 0.1)
    -stouten   - use Stouten radii instead of ProtOr radii set
    -e <file>  - location of executable (get_volume.exe)
    
    [options] for other commands:
    -d             - process all files in directory <input file>
    -o             - output directory
    -surface       - use surface atoms only
    -buried        - use non-surface atoms only
    -cavnb         - use cavity neighbor atoms only
    -no_cavnb      - discard atoms neighboring cavities
    -bins <num>    - number of bins for surfdist command [20]
    -maxd <value>  - maximum distance for surfdist command [15.0]
    -h             - include hetero groups 
    -w             - include water     
    -r <file>      - specify reference file for 'compare'
    -volume        - calculate atomic volumes instead of packing densities
    -html          - HTML output of reports.
    -native        - use 173 native residuename_atomname atom types
    """
    
    undocumented_options="""
    -packing   - calculate packing densities (default)
    -minuspackdens - substitute bfactors with 1-packing density.
    -protor        - use concise ProtOr atom typing scheme with 19 types
    -??????        - use ProtOr radii set (Tsai 1999) (default)       -no_hetnb      - discard atoms neighboring hetero groups (0) 
    """

    command = None
    input_file = None
    out_path = os.getcwd()+os.sep
    options = get_options()
    
    a = sys.argv
    if len(a) > 1:
        command = a[1]
        if len(a)>2: input_file = a[-1]
        i = 2
        while i < len(a)-1:
            if i+1 < len(a): next = a[i+1]
            else: next = None
            if a[i] == '-o' and next:
                    out_path = next
                    if out_path[-1] != os.sep: out_path += os.sep
                    i += 1
            elif a[i] == '-d' and next:
                    input_file = next
                    if input_file[-1] != os.sep:
                        input_file += os.sep 
                    options['dataset'] = input_file
                    i += 1
            elif a[i] == '-g' and next:
                    options['grid_dist'] = float(next)
                    i += 1
            elif a[i] == '-h': options['include_hetero'] = 1
            elif a[i] == '-nohet': options['include_hetero'] = 0
            elif a[i] == '-w': options['include_water'] = 1
            elif a[i] == '-buried':
                options['discard_surface'] = 1
                options['discard_buried'] = 0
            elif a[i] == '-surface':
                options['discard_surface'] = 0
                options['discard_buried'] = 1
            elif a[i] == '-no_hetnb': options['discard_hetero_neighbors'] = 1
            elif a[i] == '-cavnb': options['discard_non_cavity_neighbors'] = 1
            elif a[i] == '-no_cavnb': options['discard_cavity_neighbors'] = 1            
            elif a[i] == '-e' and next:
                    options['executable'] = next
                    i += 1
            elif a[i] == '-bins' and next:
                    options['surfdist_bins'] = int(next)
            elif a[i] == '-maxd' and next:
                    options['surfdist_maxdepth'] = float(next)
            elif a[i] == '-r' and next:
                    options['reference_file'] = next
                    i += 1
            elif a[i] == '-html': options['reportmode'] = 'html'
            elif a[i] == '-minuspackdens': options['bfactor'] = 'minuspackdens' 
            elif a[i] == '-packing': options['mode'] = 'packing_density'
            elif a[i] == '-volume': options['mode'] = 'total_volume'
            elif a[i] == '-stouten': options['radii'] = 'stouten'
            elif a[i] == '-protor': options['atomtyping'] = 'protor' 
            elif a[i] == '-native': options['atomtyping'] = 'native' 
            i += 1

    message = process_command(command,input_file,out_path,options)
    if message:
	print message[0]
    else:
        print doc


    
