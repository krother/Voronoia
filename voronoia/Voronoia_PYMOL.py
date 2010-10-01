#!/usr/bin/python

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

N_BINS = 11

from voro_modules.Voronoia_GUI import PackTools
from voro_modules import VolParser, AveragePacking, ZScore, MoleculeDataset

from pymol import cmd

import tempfile, os, string, re
import Tkinter
from Tkinter import *
import tkFileDialog, tkFont, Pmw

# set this to False if you don't want to have Voronoia pop up 
# each time you start PyMOL
AUTOMATICALLY_OPEN_PYMOL_PLUGIN = True
# AUTOMATICALLY_OPEN_PYMOL_PLUGIN = False


ptools = None
def start_voronoia(self):
    ptools = PackTools(self,pymol_gui)
    
def __init__(self):
    """
    This is the part that is automatically invoked on PyMOL startup.
    """    
    self.menuBar.addmenuitem('Plugin', 'command',
                      'Voronoia',
                      label='Voronoia',
                      command = lambda: start_voronoia(self))

    if AUTOMATICALLY_OPEN_PYMOL_PLUGIN:
        start_voronoia(self)

    cmd.extend('packing',show_packing)
    cmd.extend('packing_zscore',show_zscore)
    cmd.extend('show_cavities',show_cavities)
    cmd.extend('show_cav_neighbors',show_cav_neighbors)

class PymolGUI:

    def __init__(self):
        # pymol tab elements
        self.pymol_colorbar = None
        self.pymol_sel = 'all'
        self.pymol_sel_var = StringVar()
        self.pymol_sel_var.set(self.pymol_sel)
        self.obj_name = 'all'
        self.derived_vol = 0
        self.col = []
        self.coldesc = []

        self.subset = 'both'
        self.packing = None        

    def make_pymol_gui(self,packing,notebook):
            self.packing = packing
            self.get_colors(N_BINS)
            page = notebook.add('PyMOL')

            #if PIL_FOUND:
            #    # bg='white',
            #    can = Canvas(page,height='55',width='180')
            #    can.pack()
            #    self.drawn = can.create_image(10,10,image=self.pymollogo,anchor=NW)           

            group = Pmw.Group(page,tag_text='Structure to display:')
            group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)

            # Create button to launch the dialog.
            frm = Frame(group.interior())
            frm.pack(fill='both',expand=1, padx = 10, pady = 2)
            w = Tkinter.Button(frm,text = 'Choose selection',command = self.select_object)
            w.pack(side=LEFT)
            rl = Label(frm,textvariable=self.pymol_sel_var,font=tkFont.Font(family="Helvetica",size="9",weight="normal"),relief='sunken', padx = 5, pady = 2)
            rl.pack(side=LEFT,fill='x', expand=1)

            ll = Label(group.interior(),text = 'OR')
            ll.pack()
            
            # input .vol file(s)
            frm = Frame(group.interior())
            frm.pack(fill='both',expand=1, padx = 10, pady = 2)
            rb = Button(frm,text='Choose .vol file',command=self.pymol_volfiledialog)
            rb.pack(side=LEFT)
            rl = Label(frm,textvariable=packing.vol_file_var,font=tkFont.Font(family="Helvetica",size="9",weight="normal"),relief='sunken', padx = 5, pady = 2)
            rl.pack(side=LEFT,fill='x', expand=1)

            # colorbar
            frame=Frame(page,background='white')
            frame.pack(side=RIGHT)
            self.pymol_colorbar = frame
            self.colorbar(10,0.0,1.0)

            group = Pmw.Group(page,tag_text='Actions')
            group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)
            rb = Button(group.interior(),text='Display packing as color scale',command=self.pymol_pd_scale)
            rb.pack(fill='x',expand=1, padx = 10, pady = 2)
            rb = Button(group.interior(),text='Display z-score as color scale',command=self.pymol_z_scale)
            rb.pack(fill='x',expand=1, padx = 10, pady = 2)
            rb = Button(group.interior(),text='Highlight cavities',command=self.pymol_cav)
            rb.pack(fill='x',expand=1, padx = 10, pady = 2)
            rb = Button(group.interior(),text='Highlight cavity neighbors',command=self.pymol_cavnb)
            rb.pack(fill='x',expand=1, padx = 10, pady = 2)

            Label(page,text='PyMOL (TM) is Copyright of DeLano Scientific 2005. \n http://www.pymol.org').pack()

    def prepare_input(self):
        """
        Checks 
        1) whether a PyMOL selection has been entered. Writes it as input file.
        2) whether the input is in PDB format, and invokes calc_volume.
        """
        if self.approot: sel = self.pymol_sel.getvalue()
        else: sel = self.pymol_sel
        if sel != "":
            temp_file = self.packing.out_path + sel+'_packing.ent'
            pymol.cmd.save(temp_file,sel)
            self.packing.input_file = temp_file

    # -------------------------------------------------------------
    def select_object(self):
        # ComboBox dialog.
        sel_dialog = Pmw.ComboBoxDialog(self.packing.parent,
            title = 'Select PyMOL object',
            buttons = ('OK', 'Cancel'),
            defaultbutton = 'OK',
            combobox_labelpos = 'n',
            label_text = 'Select PyMOL object',
            scrolledlist_items = ['all']+cmd.get_names('all') # 'models','selections'
            )
        sel_dialog.withdraw()
        sel_dialog.activate()
        self.pymol_sel = sel_dialog.get()
        self.pymol_sel_var.set(self.pymol_sel)
        self.packing.vol_file = None
        self.packing.vol_file_var.set("")        

    def pymol_volfiledialog(self):
        self.derived_vol = 0
        return self.packing.vol_filedialog()
        
    def get_pymol_input(self):
        """Returns a .vol file for use by PyMOl functions"""
        if not self.packing.vol_file:
            # save selection as PDB
            temp_file = tempfile.mktemp()+'.pdb'
            cmd.save(temp_file,self.pymol_sel)
            # calculate .vol file of it.
            prev_inp = self.packing.input_file
            self.packing.input_file = temp_file
            self.packing.calc_packing(show_message=0)
            self.packing.input_file = prev_inp
            self.obj_name = self.pymol_sel + '_pack'
            self.packing.vol_file_var.set(self.obj_name)
            self.derived_vol = 1
            mode = 1
        else:
            if not self.derived_vol:
                self.obj_name = string.split(self.packing.vol_file,os.sep)[-1]
                self.obj_name = re.sub('\.vol\Z','',self.obj_name)
            mode = 0

        return self.packing.vol_file,mode
        
    def pymol_pd_scale(self):        
        """Colors a molecule in pymol according to packing density."""
        vol_file,from_selection = self.get_pymol_input()
        if not self.packing.vol_file:
            self.packing.select_first_message('PyMOL selection or .vol file')
            return
        
        self.packing.get_options()
        self.packing.options['bfactor'] = 'packing_density' 
        cmd.delete(self.obj_name)

        # copy options, to create some special selections
        diff_options = self.packing.options.copy()
        diff_options['discard_buried'] = False
        diff_options['discard_surface'] = False
        diff_options['discard_hetero_neighbors'] = False
        diff_options['discard_cavity_neighbors'] = False
        diff_options['discard_non_cavity_neighbors'] = False

        # 1) load complete object, if its not there yet 
        #if not from_selection:
        complete_name = self.obj_name
        self.load_vol_to_pymol(self.packing.vol_file,complete_name,diff_options)
        cmd.color('gray',complete_name)
        #else:
        #    complete_name = self.pymol_sel

        # 2) define surface selection
        diff_options['discard_buried'] = True
        self.load_vol_to_pymol(self.packing.vol_file,complete_name+'_surface',diff_options,complete_name)

        # 3) define buried selection
        diff_options['discard_buried'] = False
        diff_options['discard_surface'] = True
        self.load_vol_to_pymol(self.packing.vol_file,complete_name+'_buried',diff_options,complete_name)

        # 4) load set of atoms where packing should be displayed
        self.load_vol_to_pymol(self.packing.vol_file,self.obj_name,diff_options)
        cmd.show('spheres',self.obj_name)
        # cmd.set('sphere_scale',0.5,obj_name)
        cmd.color('gray',self.obj_name)

        # self.color_bfactor(self.obj_name,'packdens')
        self.color_bfactor(complete_name+'_buried','packdens')
        cmd.do('set sphere_scale,0.5,%s'%self.obj_name)
        cmd.do('indicate none')
        cmd.orient()        

    def load_vol_to_pymol(self,filename,obj_name,options,map_obj=None):
        temp_file = tempfile.mktemp()+'.pdb'
        # create PDB file with packing density as B-factor
        mol = VolParser.VolParser(self.packing.vol_file)        
        # mol.parse_vol_file(options)
        # mol.write_pdb_file(temp_file,options)
        mol.replace_bfactors(self.packing.vol_file,temp_file,options,True)
        if map_obj:
            cmd.load(temp_file,'dummy_vol')
            cmd.select(obj_name,'dummy_vol expand 0.01 and %s'%(map_obj))
            cmd.delete('dummy_vol')
        else:
            cmd.load(temp_file,obj_name)        
        os.remove(temp_file)
        

    def pymol_z_scale(self):
        vol_file,from_selection = self.get_pymol_input()
        if not self.packing.vol_file:
            self.packing.select_first_message('PyMOL selection or .vol file')
            return
        
        self.packing.get_options()
        self.packing.options['bfactor'] = 'zscore' 
        cmd.delete(self.obj_name)

        # copy options, to create some special selections
        diff_options = self.packing.options.copy()
        diff_options['discard_buried'] = False
        diff_options['discard_surface'] = False
        diff_options['discard_hetero_neighbors'] = False
        diff_options['discard_cavity_neighbors'] = False
        diff_options['discard_non_cavity_neighbors'] = False

        # 1) load complete object, if its not there yet 
        # if not from_selection:
        complete_name = self.obj_name
        self.load_vol_to_pymol(self.packing.vol_file,complete_name,diff_options)
        cmd.color('gray',complete_name)
        #else:
        #    complete_name = self.pymol_sel

        # 2) define surface selection
        diff_options['discard_buried'] = True
        self.load_vol_to_pymol(self.packing.vol_file,complete_name+'_surface',diff_options,complete_name)

        # 3) define buried selection
        diff_options['discard_buried'] = False
        diff_options['discard_surface'] = True
        self.load_vol_to_pymol(self.packing.vol_file,complete_name+'_buried',diff_options,complete_name)

        # 4) load set of atoms where packing should be displayed
        self.load_vol_to_pymol(self.packing.vol_file,self.obj_name,diff_options)
        cmd.show('spheres',self.obj_name)
        cmd.set('sphere_scale',0.5,self.obj_name)
        cmd.color('gray',self.obj_name)

        cmd.show('spheres',self.obj_name)
        cmd.set('sphere_scale',0.5,self.obj_name)
        cmd.color('gray',self.obj_name)
        
        self.color_bfactor(self.obj_name,'zscore')
        cmd.do('indicate none')

        cmd.orient()

    def pymol_cav(self):
        vol_file = self.get_pymol_input()
        if not self.packing.vol_file:
            self.packing.select_first_message('PyMOL selection or .vol file')
            return
        
        cavobject = "cavities"
        cmd.delete(cavobject)          
        cavfile = self.packing.write_cav(0)
        temp_file = tempfile.mktemp()+'.pdb'

        cmd.load(cavfile,cavobject)
        cmd.show('spheres',cavobject)
        cmd.set('sphere_scale',1.0,cavobject)        
        cmd.color('tv_orange',cavobject)            
        cmd.orient()    

    def pymol_cavnb(self):
        vol_file = self.get_pymol_input()
        if not self.packing.vol_file:
            self.packing.select_first_message('PyMOL selection or .vol file')
            return
        temp_file = tempfile.mktemp()+'.pdb'
        
        cavnbobject = "cav_neighbors"
        cmd.delete(cavnbobject)
        cavnbfile = self.packing.write_cavnb(0)
        cmd.load(cavnbfile,cavnbobject)
        cmd.show('spheres',cavnbobject)
        cmd.set('sphere_scale',0.8,cavnbobject)
        cmd.color('slate',cavnbobject)
        cmd.orient()

    def get_colors(self,nbins,sat=1.0,value=1.0):
        import colorsys
        self.col=[]
        self.coldesc=[]
        for j in range(nbins):
          self.coldesc.append('col' + str(j))
          # create colors in a gradient from red through magenta to blue
          rgb = [min(1.0, float(j)*2/(nbins-1)), 0.0, min(1.0, float(nbins-j-1)*2/(nbins-1))]

          # convert rgb to hsv,  modify saturation and value and convert back to rgb
          hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
          hsv[1] *= sat
          hsv[2] *= value

          #hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
          #convert to rgb and append to color list
          rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
          self.col.append(rgb)
          
          cmd.set_color(self.coldesc[j],self.col[j])
                
    def colorbar(self,nbins,vmin,vmax):
        my_hex = lambda x:('00'+str(hex(int(255*x)))[2:])[-2:]
        value = vmin
        # remove old stuff
        for s in self.pymol_colorbar.slaves():
            s.destroy()
            
        for i in range(nbins):
            j = nbins-i-1
            hexc = '#'+my_hex(self.col[j][0])+my_hex(self.col[j][1])+my_hex(self.col[j][2])
            label = Label(self.pymol_colorbar,text='%4.2f'%(value),bg=hexc,fg='white')
            label.pack(side=BOTTOM,fill='both',expand=1,ipadx=4,ipady=4)
            value += (vmax-vmin)/nbins

    # rlc color_b.py version 7.0
    # original code 2004 by Robert L. Campbell
    def color_bfactor(self,selection,scale,minimum='',maximum='',mode='hist'):
      """
        based on code from Robert L. Campbell with enhancements from James Stroud
      """
      sat = 1.0
      value = 1.0
      nbins = N_BINS

      # get list of occupancies from selection
      # (PyMOL is unable to read b factors from .vol files)
      m = cmd.get_model(selection)
      sel = []

      if len(m.atom) == 0:
        print "Sorry, no atoms selected"
        b_list = []
      else:
        b_list = [m.atom[i].q for i in range(len(m.atom))]

        max_b = max(b_list)
        min_b = min(b_list)

        if mode == 'ramp':
          # color in bins of equal numbers of atoms
          # subtract 0.1 from the lowest B in order to ensure that the single
          # atom with the lowest B value doesn't get omitted
          b_list.sort()
          bin_num = int(len(b_list)/nbins)
          sel.append(selection + " and (b < %4.4g" % (b_list[bin_num]) + " or b = %4.4g" % (b_list[bin_num]) + ")")
          for j in range(1,nbins):
            sel.append(selection + " and b > %4.4g" % (b_list[j*bin_num]))

        elif mode == 'hist':
          # check if minimum or maximum was specified and use the entered values
          if minimum != '': min_b = float(minimum)
          if maximum != '': max_b = float(maximum)
          # histogram:
          # color in bins of equal B-value ranges
          # subtract 0.1 from the lowest B in order to ensure that the single
          # atom with the lowest B value doesn't get omitted
          bin_width = (max_b - min_b)/nbins
          sel.append(selection + " and (b < %4.4g" % (min_b + bin_width) + " or b = %4.4g" % (min_b + bin_width) + ")")
          for j in range(1,nbins):
            sel.append(selection + " and b > %4.4g" % (min_b + j*bin_width))

        for j in range(nbins): cmd.color(self.coldesc[j],sel[j])
        self.colorbar(nbins,min_b,max_b)


pymol_gui = PymolGUI()

def show_packing(sel='(all)',subset='both'):
    pymol_gui.pymol_sel = re.sub('\(|\)','',sel)
    if subset == 'buried': pymol_gui.subset = 'buried atoms'
    elif subset == 'surface': pymol_gui.subset = 'surface atoms'
    elif subset == 'both': pymol_gui.subset = 'both'
    pymol_gui.pymol_pd_scale()
    
def show_zscore(sel='(all)',subset='both'):
    pymol_gui.pymol_sel = re.sub('\(|\)','',sel)
    if subset == 'buried': pymol_gui.subset = 'buried atoms'
    elif subset == 'surface': pymol_gui.subset = 'surface atoms'
    elif subset == 'both': pymol_gui.subset = 'both'
    pymol_gui.pymol_z_scale()
    
def show_cavities(sel='(all('):
    pymol_gui.pymol_sel = re.sub('\(|\)','',sel)
    pymol_gui.pymol_cav()

def show_cav_neighbors(sel='(all)'):
    pymol_gui.pymol_sel = re.sub('\(|\)','',sel)
    pymol_gui.pymol_cavnb()

