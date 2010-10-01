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

#
# Some of the GUI code was inspired from the APBS_tools by Michael G. Lerner
# whose great contribution to PyMOL is hereby acknowledged!
#

from __future__ import division
from __future__ import generators

import os, string, re, math, tempfile
import Tkinter
from Tkinter import *
import tkFileDialog, tkFont, Pmw

import Voronoia
import VolParser, AveragePacking, ZScore, MoleculeDataset
import Report, ReportDialog


class PackTools:

    def __init__(self,app, pymol_gui=None):
        # option settings
        self.options = Voronoia.get_options()
        self.options['atomtyping'] = 'few'
        self.options['mode'] = 'packing_density'

        # calculate tab elements
        self.input_file = ""
        self.input_file_var = StringVar()
        self.input_file_var.set(self._short_name(self.input_file))

        self.pdb_dataset = None
        self.pdb_dataset_radio = None
        
        self.out_path = os.getcwd()+os.sep
        self.out_path_var = StringVar()
        self.out_path_var.set(self.out_path)

        self.exe_file_var = StringVar()
        self.exe_file_var.set(self._short_name(self.options['executable']))

        self.grid_dist_sel = None

        # options tab elements
        self.subset_radio = None
        self.cavnb_radio  = None
        self.packing_mode_radio = None
        self.atomtyping_radio = None
        self.include_radio = None
        self.radii_radio = None

        # report tab elements
        self.vol_file = ""
        self.vol_file_var = StringVar()
        self.vol_file_var.set(self._short_name(self.vol_file))

        self.vol_dataset = None
        self.vol_dataset_radio = None

        self.ref_file_var = StringVar()
        self.ref_file_var.set(self._short_name(self.options['reference_file']))  


        # now create the dialog    	    
        self.approot = app.root
        dialog = Pmw.Dialog(self.approot,
                                 buttons = ('About','Help', 'Quit'), 
                                 title = 'Voronoia Packing Toolbox',
                                 command = self.execute)

        dialog.withdraw()
        parent = dialog.interior()
        self.parent = parent
        self.dialog = dialog

        self.misc_text = None
        self.misc_html = None
        self.misc_dialog = None
        

	# add logo
        if PIL_FOUND:
            self.voroidlogo = PhotoImage(file = "data"+os.sep+"voronoia.gif") 
            self.pymollogo = PhotoImage(file = "data"+os.sep+"pymol_logo.gif") 
            # ,bg='white'
            can = Canvas(parent,height='96',width='280')
            can.pack()
            self.drawn = can.create_image(0,0,image=self.voroidlogo,anchor=NW)

        # -------------- notebook tabs ----------------------
        notebook = Pmw.NoteBook(parent)
        notebook.pack(fill = 'both', expand = 1, padx = 0, pady = 0)

        # ----- PyMOL stuff --------------------------------------
        if pymol_gui:
            pymol_gui.make_pymol_gui(self,notebook)

        # ----- files --------------------------------------
        page = notebook.add('Calculate Packing')
        notebook.tab('Calculate Packing').focus_set()

        # input file
        frm = Frame(page)
        frm.pack(fill='x',expand=1, padx = 10, pady = 2)
        rb = Button(frm,text='Input PDB file',command=self.in_filedialog)
        rb.pack(side=LEFT)
        rl = Label(frm,textvariable=self.input_file_var,font=tkFont.Font(family="Helvetica",size="9",weight="normal"),relief='sunken', padx = 5, pady = 2)
        rl.pack(side=LEFT,fill='x', expand=1, padx=10, pady=1)

        # checkbox for dataset processing
        frm = Frame(page)
        frm.pack(fill='x',expand=0)
        self.pdb_dataset_radio = Pmw.RadioSelect(frm,buttontype="checkbutton",orient='vertical')
        self.pdb_dataset_radio.pack(side=LEFT,fill='both', expand=0, padx=10, pady=1)
        self.pdb_dataset_radio.add('Process all PDB files in this directory',command=None)
        
        # output path
        frm = Frame(page)
        frm.pack(fill='both',expand=1, padx = 10, pady = 2)
        rb = Button(frm,text='Output directory',command=self.out_pathdialog)
        rb.pack(side=LEFT)
        rl = Label(frm,textvariable=self.out_path_var,justify=LEFT,wraplength=180,font=tkFont.Font(family="Helvetica",size="9",weight="normal"),relief='sunken', padx = 10, pady = 2)
        rl.pack(side=LEFT,fill = 'x', expand = 1,padx=3)

        # start calculation button
        frm = Frame(page)
        frm.pack(fill='both',expand=1, padx = 10, pady = 0)
        button = Button(frm,text='Calculate packing file(s)',command=self.calc_packing)
        button.pack(side=LEFT,fill='x',expand=1,pady=2)

        # calculation parameters        
        group = Pmw.Group(page,tag_text='Calculation parameters')
        group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)

        frm = Frame(group.interior())
        frm.pack(fill='both',expand=1, padx = 0, pady = 0)
        # grid distance selector
        grid_dist = Pmw.OptionMenu(frm,
                    labelpos = 'w',label_text = 'Grid distance [Angstrom]',
                    items = ('0.20','0.10','0.05','0.02','0.01',),
                    )
        grid_dist.pack(side=LEFT,fill='none',expand=0,padx=10)
        self.grid_dist_sel = grid_dist
        
        frm = Frame(group.interior())
        frm.pack(fill='both',expand='1')
        # radii set radiobuttons
        radio = Pmw.RadioSelect(frm,buttontype="radiobutton",orient='vertical')
        radio.pack(side=LEFT)
        radio.add('ProtOr radii (Tsai 1999)')
        radio.add('Stouten radii (Stouten 1993)')
        radio.setvalue('ProtOr radii (Tsai 1999)')
        self.radii_radio = radio
        # hetero/water checkboxes
        radio = Pmw.RadioSelect(frm,buttontype="checkbutton",orient='vertical')
        radio.pack(side=RIGHT) # ,pady=3
        radio.add('Include hetero groups',command=None)
        radio.add('Include water molecules',command=None)
        self.include_radio = radio

        # change executable button
        frm = Frame(group.interior())
        frm.pack(fill='both',expand='1',padx=10)
        rb = Button(frm,text='Change executable',command=self.exe_filedialog)
        rb.pack(side=LEFT,fill='x',expand=0)
        
        # --------------- Options ----------------
        page = notebook.add('Options')

        # buried / surface atoms
        group = Pmw.Group(page,tag_text='Solvent acessibility')
        group.pack(fill='both',expand=1, padx=4, pady=2)
        frm = Frame(group.interior())
        frm.pack(fill='both',expand='1')
        radio = Pmw.RadioSelect(frm,buttontype="checkbutton",orient='vertical')
        radio.pack(side=LEFT)
        radio.add('include buried atoms',command=None)
        radio.add('include solvent accessible atoms',command=None)
        radio.setvalue(['include buried atoms'])
        self.subset_radio = radio

        # cavity-neighbor / non-neighbor atoms
        group = Pmw.Group(page,tag_text='Cavity neighbors')
        group.pack(fill='both',expand=1, padx=4, pady=2)
        frm = Frame(group.interior())
        frm.pack(fill='both',expand='1')
        radio = Pmw.RadioSelect(frm,buttontype="checkbutton",orient='vertical')
        radio.pack(side=LEFT)
        radio.add('include cavity neighbor atoms',command=None)
        radio.add('include all other atoms',command=None)
        radio.setvalue(['include cavity neighbor atoms','include all other atoms'])
        self.cavnb_radio = radio

        # write packing densities or volumes
        group = Pmw.Group(page,tag_text='Output data type')        
        group.pack(side=TOP,fill='both',expand=1, padx=4, pady=2)
        radio = Pmw.RadioSelect(group.interior(),buttontype="radiobutton",orient='vertical')
        radio.pack(side=LEFT)
        radio.add('atomic packing densities')
        radio.add('atomic volumes')
        radio.setvalue('atomic packing densities')
        self.packing_mode_radio = radio

        # concise or long atom types
        group = Pmw.Group(page,tag_text='Atom typing')        
        group.pack(side=TOP,fill='both',expand=1, padx=4, pady=2)
        radio = Pmw.RadioSelect(group.interior(),buttontype="radiobutton",orient='vertical')
        radio.pack(side=LEFT)
        radio.add('19 ProtOr atom types')
        radio.add('173 explicit atom types')
        radio.setvalue('19 ProtOr atom types')
        self.atomtyping_radio = radio

        # surfdist parameters       
        group = Pmw.Group(page,tag_text='Surface distance parameters')
        group.pack(fill = 'both', expand = 1, padx = 0, pady = 0)
        frm = Frame(group.interior())
        frm.pack(fill='both',expand=1, padx = 0, pady = 0)
        n_bins = Pmw.OptionMenu(frm,
                    labelpos = 'w',label_text = 'number of bins',
                    items = ('10','20','50','100',),
                    initialitem='20'
                    )
        max_depth = Pmw.OptionMenu(frm,
                    labelpos = 'w',label_text = 'maximum depth [Angstrom]',
                    items = ('5.0','10.0','15.0','20.0'),
                    initialitem='10.0'
                    )
        n_bins.pack(side=LEFT,fill='none',expand=0,padx=10)
        max_depth.pack(side=LEFT,fill='none',expand=0,padx=10)
        self.n_bins_sel = n_bins
        self.maxdepth_sel = max_depth
        # Pmw.Color.changecolor(page, background = 'white',foreground='#000070')

        # ------------- actions -----------------
        page = notebook.add('Results')

        # input .vol file(s)
        group = Pmw.Group(page,tag_text='Packing files')
        group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)
        frm = Frame(group.interior())
        frm.pack(fill='x',expand=1, padx = 10, pady = 2)
        rb = Button(frm,text='Choose .vol file',command=self.vol_filedialog)
        rb.pack(side=LEFT)
        rl = Label(frm,textvariable=self.vol_file_var,font=tkFont.Font(family="Helvetica",size="9",weight="normal"),relief='sunken', padx = 5, pady = 2)
        rl.pack(side=LEFT,fill='x', expand=1, padx=10, pady=1)
        # checkbox for dataset processing
        frm = Frame(group.interior())
        frm.pack(fill='x',expand=0)
        self.vol_dataset_radio = Pmw.RadioSelect(frm,buttontype="checkbutton",orient='vertical')
        self.vol_dataset_radio.pack(side=LEFT,fill='both', expand=0, padx=10, pady=1)
        self.vol_dataset_radio.add('Process all .vol files in this directory',command=None)

        # Reports
        group = Pmw.Group(page,tag_text='Report')
        group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)
        rb = Button(group.interior(),text='Create packing report',command=self.packing_report)
        rb.pack(fill='x',expand=1, padx = 10, pady = 2)
        rb = Button(group.interior(),text='Calculate average packing',command=self.write_avgpd)
        rb.pack(fill='x',expand=1, padx = 10, pady = 2)

        # compare to reference values
        group = Pmw.Group(page,tag_text='Reference data')
        group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)
        frm = Frame(group.interior())
        frm.pack(fill='x',expand=1, padx = 10, pady = 2)
        rb = Button(frm,text='Choose reference file',command=self.ref_filedialog)
        rb.pack(side=LEFT)
        rl = Label(frm,textvariable=self.ref_file_var,font=tkFont.Font(family="Helvetica",size="9",weight="normal"),relief='sunken', padx = 5, pady = 2)
        rl.pack(side=LEFT,fill='x', expand=1, padx=10, pady=1)
        rb = Button(group.interior(),text='Compare .vol file to reference data',command=self.compare_ref)
        rb.pack(fill='x',expand=1, padx = 10, pady = 2)

        # report surface distances
        group = Pmw.Group(page,tag_text='Surface distance statistics')
        group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)
        rb = Button(group.interior(),text='Report distance to surface',command=self.surfdist)
        rb.pack(fill='x',expand=1, padx = 10, pady = 2)
        
        # write cavity coordinates
        group = Pmw.Group(page,tag_text='PDB output')
        group.pack(fill = 'x', expand = 1, padx = 10, pady = 2)
        rb = Button(group.interior(),text='Write cavities as PDB',command=self.write_cav)
        rb.pack(fill='x',expand=1, padx = 10, pady = 2)
        rb = Button(group.interior(),text='Write cavity neighbors as PDB',command=self.write_cavnb)
        rb.pack(fill='x',expand=1, padx = 10, pady = 2)
                        
        # Pmw.Color.changecolor(page, background = 'white',foreground='#000070')
        notebook.recolorborders()
        notebook.setnaturalsize()

        dialog.show()
    
    # -------------------------------------------------------------
    def filedialog(self,titles,formats=None,var=None,directory=0):
        if not directory:
            name = tkFileDialog.askopenfilename(title=titles,filetypes=formats)
            if os.sep == "\\":
                name = re.sub('\/','\\\\',name)
            shortname = self._short_name(name)
            if var and name: var.set(shortname)
        else:
            name = tkFileDialog.askdirectory(title=titles)
            if os.sep == "\\":
                name = re.sub('\/','\\\\',name)
            if var and name: var.set(name)
            if name[-1]!=os.sep: name += os.sep
        #
        # Here strange bullshit happening on Win is caught.
        # sometimes, Python supplies me with / as file separators,
        # sometimes with \. Crappy.
        #
        return name
        
    def _short_name(self,name):
        """Cuts off all separators from the beginning of a filname."""
        short = name.split(os.sep)[-1]
        return short
    # -------------------------------------------------------------
    def in_filedialog(self):
        myFormats = [('PDB files','*.pdb *.ent')]        
        name = self.filedialog('Please choose an input PDB file',myFormats,self.input_file_var)  
        if name: self.input_file = name
    # -------------------------------------------------------------
    def vol_filedialog(self):
        myFormats = [('Packing files','*.vol')]   
        name = self.filedialog('Please choose an input packing file',myFormats,self.vol_file_var)
        if name: self.vol_file = name
    # -------------------------------------------------------------
    def ref_filedialog(self):
        myFormats = [('Average packing densities','*.avg')]
        name = self.filedialog('Please choose a reference file',myFormats,self.ref_file_var)
        if name: self.ref_file = name
    # -------------------------------------------------------------
    def out_pathdialog(self):
        name = self.filedialog('Please choose an ouput directory',None,self.out_path_var,1)
        if name: self.out_path = name
    # -------------------------------------------------------------
    def exe_filedialog(self):
        myFormats = [('Executable','*.exe, *.bat, *.sh, *.py, *.pl, *.com')]
        name = self.filedialog('Please choose executable for packing calculation',myFormats,self.exe_file_var)
        if name:
            tokens = name.split(os.sep)
            self.exe_file = tokens[-1]
            self.options['executable_path'] = tokens[:-1].join(os.sep)


    # -------------------------------------------------------------
    def execute(self, result):
        if result == 'Help':
            # Create the dialog.
            dialog = Pmw.TextDialog(self.parent, scrolledtext_labelpos = 'n',
                title = 'Help on %s'%APPLICATION_NAME,
                defaultbutton = 0,
                label_text = 'Help on %s'%APPLICATION_NAME)
            dialog.withdraw()
            dialog.insert('end', open(HELP_FILE).read())
            dialog.configure(text_state = 'disabled')
        
            dialog.activate()

        elif result == 'About':
            Pmw.aboutversion('1.0')
            Pmw.aboutcopyright('Copyright Andrean Goede\n and Kristian Rother 2007\nAll rights reserved')
            Pmw.aboutcontact(open(ABOUT_FILE).read())
            about = Pmw.AboutDialog(self.parent,applicationname=APPLICATION_NAME)
            
        elif result == 'Quit' or result == None:
            if __name__ == '__main__':
                self.approot.destroy()
            else:
                self.dialog.withdraw()

    # --------- select vol file first message -----------------
    def select_first_message(self,filetype):
        ms = Pmw.MessageDialog(self.parent,title="Missing file",buttons=['OK'],iconpos='w',icon_bitmap='warning',message_text="Select %s first."%filetype)
        ms.withdraw()
        ms.activate()
        ms.deactivate()

    # --------- text dialog for various purposes ---------------
    def text_dialog(self,titles,labels,text,html):
        self.misc_text = text
        self.misc_html = html

        # Create a dialog window.
        dialog = ReportDialog.ReportDialog(self.parent,titles,labels,text,self.exe_text_dialog)
        self.misc_dialog = dialog
        dialog.withdraw()
        dialog.activate()
        
    def exe_text_dialog(self, result):
        if result == 'Save as HTML': self.write_html()
        elif result == 'Save as ASCII': self.write_ascii()
        else:
	    self.misc_dialog.deactivate()
	    self.misc_dialog = None
	    self.misc_text = None
	    self.misc_html = None
        
    def write_html(self):
        myFormats = [('HTML file','*.html')]        
        name = tkFileDialog.asksaveasfilename(title='Please select file name',filetypes=myFormats)
        if name and len(name)>5:
            if name[-5:]!='.html': name += '.html'
            open(name,'w').write(self.misc_html)
        
    def write_ascii(self):
        myFormats = [('Text file','*.txt')]        
        name = tkFileDialog.asksaveasfilename(title='Please select file name',filetypes=myFormats)
        if name and len(name)>4:
            if name[-4:]!='.txt': name += '.txt'
            open(name,'w').write(self.misc_text)

    #
    # now the action buttons
    #
    def get_vol_files(self):
        # get list of input files
        vol_files = [self.input_file]
        if self.pdb_dataset: vol_files = [self.pdb_dataset + f for f in os.listdir(self.pdb_dataset)]
        return vol_files
        
    # ------------ calculate packing files ---------------------
    def calc_packing(self,show_message=1):
        # check for complete data
        self.get_options()
        if not self.input_file:
            self.select_first_message("PDB file")
            return
        if not self.out_path:
            self.select_first_message("output directory")
            return

        result = []
        for in_file in self.get_vol_files():
            # show calculation message
            ms = Pmw.MessageDialog(title="Calculation running",buttons=[],message_text="""Packing data for the structure
%s
is being calculated. Please wait.

The calculation takes 1-15 minutes per file, depending on the size of the structure
and your hardware."""%in_file)
            ms.show()
            self.approot.update()

            # calculate
            message = Voronoia.process_command('calc',in_file,self.out_path,self.options)
            result.append(message)

            # destroy old message
            ms.withdraw()
            ms.destroy()
            self.approot.update()            

        # set vol file variable
        if self.pdb_dataset:
            text = "%i packing files were calculated:\n"%len(result)
            text += string.join(['\t'+name[1] for name in result],'\n')
            if show_message:
                ms = Pmw.MessageDialog(title=".vol files calculated",message_text=text)
            self.vol_dataset = self.out_path
            self.vol_file = result[0][1]
            self.vol_file_var.set(self.vol_dataset)
            self.vol_dataset_radio.setvalue(['Process all .vol files in this directory'])
            return [name[1] for name in result]
        else:
            message = result[0]
            if show_message:
                ms = Pmw.MessageDialog(title=".vol file calculated",message_text=message[0])
            self.vol_file = message[1]            
            self.vol_file_var.set(self._short_name(message[1])) 
            return message[1]

    # ------------ report packing ------------------------------
    def packing_report(self):
        if not self.vol_file:
            self.select_first_message(".vol file")
            return
        self.get_options()
        self.options['reportmode']='both'

        if self.vol_dataset:
            self.options['dataset'] = self.vol_dataset
            text,html = Voronoia.process_command('report',self.vol_dataset,self.vol_dataset,self.options)
        else:
            text,html = Voronoia.process_command('report',self.vol_file,self.out_path,self.options)
        self.text_dialog('Packing report','Packing report for %s'%(self.vol_file),text,html)

	
    # ------------ calculate z-scores --------------------------    
    def compare_ref(self):
        if not self.vol_file:
            self.select_first_message(".vol file")
            return
        if not self.options['reference_file']:
            self.select_first_message("reference file")
            return
        self.get_options()
        self.options['reportmode']='both'

        if self.vol_dataset:
            self.options['dataset'] = self.vol_dataset
            text,html = Voronoia.process_command('compare',self.vol_dataset,self.vol_dataset,self.options)
        else:
            text,html = Voronoia.process_command('compare',self.vol_file,self.out_path,self.options)        
        self.text_dialog('Packing compared to reference','Packing compared to reference values in\n%s'%(self.options['reference_file']),text,html)

    # ------------ calculate average packing -------------------
    def write_avgpd(self):
        if not self.vol_file:
            self.select_first_message(".vol file")
            return
        self.get_options()
        self.options['reportmode']='both'

        if self.vol_dataset:
            self.options['dataset'] = self.vol_dataset
            text,html = Voronoia.process_command('average',self.vol_dataset,self.vol_dataset,self.options)
        else:
            text,html = Voronoia.process_command('average',self.vol_file,self.out_path,self.options)
        self.text_dialog('Average packing','Average packing values for data from %s'%(self.vol_file),text,html)

    # ------------ write surface distances ----------------------
    def surfdist(self):
        if not self.vol_file:
            self.select_first_message(".vol file")
            return
        self.get_options()
        self.options['reportmode']='both'

        if self.vol_dataset:
            self.options['dataset'] = self.vol_dataset
            text,html = Voronoia.process_command('surfdist',self.vol_dataset,self.vol_dataset,self.options)
        else:
            text,html = Voronoia.process_command('surfdist',self.vol_file,self.out_path,self.options)
        self.text_dialog('Surface distance','Surface distance values for data from %s'%(self.vol_file),text,html)
        
    # ------------ write cavities -------------------------------
    def write_cav(self,response=1):
        if not self.vol_file:
            self.select_first_message(".vol file")
            return

        self.get_options()
        if self.vol_dataset:
            self.options['dataset'] = self.vol_dataset
            message = Voronoia.process_command('cavities',self.vol_dataset,self.vol_dataset,self.options)
        else:    
            message = Voronoia.process_command('cavities',self.vol_file,self.out_path,self.options)
        if response:
            ms = Pmw.MessageDialog(self.parent,title="Cavities written as PDB file",message_text=message[0])
            ms.activate()
            ms.withdraw()
        return message[1]

    # ------------ write cavity neighbors -----------------------
    def write_cavnb(self,response=1):
        if not self.vol_file:
            self.select_first_message(".vol file")
            return

        self.get_options()
        if self.vol_dataset:
            self.options['dataset'] = self.vol_dataset
            message = Voronoia.process_command('cavneighbors',self.vol_dataset,self.vol_dataset,self.options)
        else:
            message = Voronoia.process_command('cavneighbors',self.vol_file,self.out_path,self.options)

        if response:
            ms = Pmw.MessageDialog(self.parent,title="Cavity neighbors written as PDB file",message_text=message[0])
            ms.activate()
            ms.withdraw()
        return message[1]
    
    # ------------ get options from dialog elements -----------------------
    def get_options(self):
        """
        Evaluates values of all control elements.
        """
        # evaluate subset checkboxes
        if self.approot: subset = self.subset_radio.getvalue()
        else: subset = self.subset
        self.options['discard_buried'] = True
        self.options['discard_surface'] = True
        for value in subset:
            if value == 'include buried atoms': self.options['discard_buried'] = False
            elif value == 'include solvent accessible atoms': self.options['discard_surface'] = False
            
        if self.approot: cavnb = self.cavnb_radio.getvalue()
        else: cavnb = self.cavnb
        self.options['discard_cavity_neighbors'] = True
        self.options['discard_non_cavity_neighbors'] = True
        for value in cavnb:
            if value == 'include cavity neighbor atoms': self.options['discard_cavity_neighbors'] = False
            elif value == 'include all other atoms': self.options['discard_non_cavity_neighbors'] = False
        
        # grid distance
        self.options['grid_dist'] = float(self.grid_dist_sel.getvalue())

        # surfdist options
        self.options['surfdist_bins'] = int(self.n_bins_sel.getvalue())
        self.options['surfdist_maxdepth'] = float(self.maxdepth_sel.getvalue())
        
        # hetero and water check boxes
        self.options['include_hetero'] = False
        self.options['include_water'] = False
        for value in self.include_radio.getvalue():
            if value == 'Include hetero groups': self.options['include_hetero'] = True
            if value == 'Include water molecules': self.options['include_water'] = True

        # kind of packing data used
        mode = self.packing_mode_radio.getvalue()
        if mode == 'atomic packing densities': self.options['mode'] = 'packing_density'
        else: self.options['mode'] = 'total_volume'

        # atom tying 
        atomtyping = self.atomtyping_radio.getvalue()
        if atomtyping == '19 ProtOr atom types': self.options['atomtyping'] = 'protor'
        else: self.options['atomtyping'] = 'all'

        # atom radii used
        radii = self.radii_radio.getvalue()
        if radii == 'ProtOr radii (Tsai 1999)': self.options['radii'] = 'protor'
        else: self.options['radii'] = 'stouten'

        # evaluate dataset checkboxes
        self.options['dataset'] = None
        self.pdb_dataset = None
        if self.pdb_dataset_radio.getvalue():
            # extract directory name from input file
            self.pdb_dataset = self.input_file.split(os.sep)[:-1]
            self.pdb_dataset = string.join(self.pdb_dataset,os.sep)
            self.pdb_dataset += os.sep
            
        self.vol_dataset = None
        if self.vol_dataset_radio.getvalue():
            # extract directory name from input file
            self.vol_dataset = self.vol_file.split(os.sep)[:-1]
            self.vol_dataset = string.join(self.vol_dataset,os.sep)
            self.vol_dataset += os.sep
            


#
# startup code
#

# Python backward-compatibility...
True = 1
False = 0

# check for image library
try:
    from PIL import Image
    PIL_FOUND = True
except:
    PIL_FOUND = False


class App:
    def my_show(self,*args,**kwargs):
        pass

if __name__ == '__main__':
    app = App()
    app.root = Tk()
    app.root.withdraw()
    widget = PackTools(app)    

    app.root.mainloop()

