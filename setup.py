from setuptools import setup, find_packages

try:
    import py2exe
except:
    pass


setup(
    name = 'Voronoia',
    packages = find_packages(),
    setup_requires = ["setuptools_git"],
    install_requires = [],
    version = '1.0.1',
    package_data = {
        'voronoia': ['data/*']
    }
)


if 'py2exe' in locals():
    #
    # in die library.zip muss unbedingt das ganze Pmw-Verzeichnis reingeschoben werden.
    #
    #
    
    data_files = [
        ('data',[
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
          data_files=data_files)
    
    setup(windows=['Voronoia_GUI.py'],
          data_files=data_files)
