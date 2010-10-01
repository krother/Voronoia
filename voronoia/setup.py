from distutils.core import setup
import py2exe

#
# in die library.zip muss unbedingt das ganze Pmw-Verzeichnis reingeschoben werden.
#
#
file_list = [('data',[
              'data\\voronoia.gif',
              'data\\protor_atomtypes.txt',
              'data\\globular.avg',
              'data\\README.TXT',
              'data\\pymol_logo.gif',
              'data\\get_volume.exe',
              'data\\about.txt',
              ])
    ]

setup(console=['Voronoia.py'],
      data_files=file_list)

setup(windows=['Voronoia_GUI.py'],
      data_files=file_list)
      
#setup(console=['Packing_GUI.py'])


