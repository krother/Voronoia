
Voronoia V 1.0.1

by Kristian Rother, Peter Hildebrand, Andrean Goede, Bjoern Gruening, and Robert Preissner

------------------------------------------------------------------------
(c) 2008 Kristian Rother (krother@rubor.de)
         and Andrean Goede
         Charite - Medical University of Berlin

This software may be distributed freely, as long as this Copyright
notice is not modified.
------------------------------------------------------------------------

When using this software, please cite:

 Rother K, Preissner R, Goede A, Froemmel C.
 Inhomogeneous molecular density: reference packing densities and distribution of cavities within proteins.
 Bioinformatics. 2003 Nov 1;19(16):2112-21. 


Acknowledgements:
    Fiete Haack
    Cornelius Froemmel

------------------------------------------------------------------------

Usage Instructions

------------------------------------------------------------------------

Voronoia is a software for analyzing the packing of protein structures. It can calculate atomic packing measures, compare them to reference values, and discover atom-sized cavities in the protein interior. Voronoia has been developed for a Windows environment. It runs on Linux machines as well. On both operating systems, Voronoia also works as a plug-in for the PyMOL molecular viewer.

------------------------------------------------------------------------
1. Installation
------------------------------------------------------------------------

1.1) On Windows

Unzip the Voronoia.zip archive. There are two directories: 'standalone' and 'plugin'. In the 'standalone' directory, 'Voronoia_GUI.exe' starts the graphical interface. You can create a link icon on your Desktop for this one. 'Voronoia.exe' provides a command-line-tool, which is best used from a text console (in the Start Menu or by executing the 'cmd'). Just type 'Voronoia.exe' there.
None of the files should be moved out of its directory, but the 'standalone' directory may be moved as a whole without problems.

To use the PyMOL plugin, all contents of the 'plugin\startup' directory should be copied to the 'PYMOL_DIR\modules\pmg_tk\startup' directory. The Voronoia GUI then appears automatically each time PyMOL is started.  Substitute 'PYMOL_DIR' by the location where PyMOL is installed on your machine, most probably 'C:\Program Files\DeLano Scientific\PyMOL\'. We recommend using the 0.99 Version of PyMOL or newer.

If you don't want Voronoia to start each time with PyMOL, change the according switch in the startup/voronoia.py file.


1.2) On Linux

Most important thing first: To calculate packing data, you will need Wine (www.winehq.com). Without Wine, you can still create all kinds of reports from existing .vol packing files. 

To install, Python 2.4 or higher and the Pmw library (1.2 or higher) must be installed. Unpack the Voronoia.tar.gz archive. To start the command-line application, type:

python Voronoia.py

For the graphical interface:

python Voronoia_GUI.py

Wine must be configured to run data/get_volume.exe (depends on your distribution). After you got it running, the 'EXECUTABLE' variable in the Voronoia.py module must be set to the proper Wine call. Of course, you may try to run Voronoia_GUI.exe from the Windows distribution directly by Wine, but we never tried this.

To use the PyMOL plugin, all contents of the 'plugin' directory should be copied to the '$PYMOL\pmg_tk\startup' directory. The Voronoia GUI then appears automatically each time PyMOL is started.  Substitute '$PYMOL' by the path where PyMOL modules are installed on your machine. This depends on your installation. It is *not* /usr/local/bin. Try places like /usr/local/pymol, /var/lib/python-support/python2.4/ or /usr/lib/python2.4/site-packages/ for the pmg_tk directory. We recommend using the 0.99 Version of PyMOL or newer.

If you don't want Voronoia to start each time with PyMOL, change the according switch in the startup/voronoia.py file.

1.3) Peeping at the Source

Unpack the archive. Installing works mostly as described in 1.2). The algorithm in 'get_volume.exe' is precompiled, so you will still require Wine on Linux.


------------------------------------------------------------------------
2. Basics
------------------------------------------------------------------------

This chapter describes the measures calculated by Voronoia.

2.1) Atomic volumes

An atomic volume is the amount of space assigned to a particular atom, delimited by its neighbors.
To calculate it, Voronoi polyhedra with hyperboloid surfaces are constructed around each atom. The Voronoi Cell procedure described in [1] uses a cubic lattice in order to assign the exact values. Different atom radii are taken into account. The procedure distinguishes between the volume inside the atoms' Van-der-Waals radius - the VdW volume V(VdW) - , and a layer of max. 1.4 Angstroem around it - the Solvent Excluded volume V(SE).


2.2) Atomic packing densities

The packing density of an atom is calculated as 

PD = V(VdW) / [V(VdW) + V(SE)],

where V(VdW) and V(SE) are the Van-der-Waals and Solvent Excluded volumes of this atom. Thus, the maximum packing density of an atom is 1.0 (none of the Solvent Excluded volume remains, which nevr occurs for real-wold data). The minimum is 0.0.


2.3) Cavities (packing defects)

These are positions, where a virtual water probe (1.4 Angstrom radius) could fit. They are defined as cavities, when there is no escape path for the probe to the surrounding solvent without collisions with any protein atoms, 

Cavities occupied by water molecules (other hetero atoms?, ANDREAN) are labeled as partially filled. If there are enough water molecules inside, that no additional probe would fit, the cavity is considered filled. Voronoia finds all cavities extensively. 


2.4) Surface and buried atoms

Packing differs principally between atoms on the protein surface and those in the inside. One should not mix these two types without a good reason. The protein surface is determined by rolling a virtual water molecule (1.4 Angstrom radius) over the protein. All atoms touched by the probe are assigned to the surface.

------------------------------------------------------------------------
3. Functions
------------------------------------------------------------------------

This chapter describes what you can do with Voronoia.

3.1) Calculate volume files.

Takes PDB files (.ent or .pdb) and calculates atomic volumes for each atom. It can be adjusted, whether water and other hetero atoms should be removed before the calculation. This affects the volume of surrounding atoms, and eventually resulting cavities. The width of the lattice used for calculation can be tuned for higher accuracy (and computation time). An .vol file with the volumes is written (description in chapter 6).


3.2) Report atomic packing densities or volumes

Creates a report in HTML or ASCII format that summmarizes the packing of a structure, and lists packing values per residue and per atom, and all cavities.


3.3) Calculate average packing values

Calculates average packing densities (or volumes) for all atoms of a structure or structure dataset. The data is grouped either by 'native' atom types (residue + atom name), or by 17 concise ProtOr types from [2] grouping similarly packed atom types together. The ProtOr atom types are labeled XmHnz, where X is the element, m the number of total bonds, n the number of bonded hydrogens and z a letter distinguishing some sub-types. Example: the C3H1u type contains sp2-carbon atoms with 3 bonds (2 single and 1 double), 1 of which is to a hydrogen.

The output is a table where for each atom type, the number of atoms, average packing density, the absolute and the relative standard deviation are listed. The output file can be used as a reference file for 3.4.


3.4) Compare a volume file to reference values.

Given a file with reference average packing densities, it can be calculated how well the packing in a .vol file fits to them. A file with reference data for buried atoms in high-resolution structures is provided, but all sets calculated with 3.3 can be used, as long as the atom typing scheme (ProtOr or native) is the same.

For the comparison, the z-score-rms is calculated:

z-score-rms = sqrt( sum( (PD(i,k) - <PD(k)>)^2 / StdDev(k) ) ),

where PD(i,k) is the packing density of the i'th atom having the type k, <PD(k)> is the average packing and StdDev(k) the standard deviation for that type. Typical values for the z-score-rms in proteins are:

   0.0       : perfect match, this is definitely an artifact.
   0.8 - 1.2 : nominal for buried atoms in proteins
   1.3 - 1.4 : slightly different from the reference, maybe a structure with low resolution.
   1.5+      : clearly different data; membrane proteins compared to soluble; loosely packed structures; erroneous structures. Reference set contains water but examined structure does not.
   3.0+      : some error in the data or in program usage, e.g. surface atoms compared with buried ones.


3.5) Write cavity positions

Creates a PDB file, where the center points of all cavities in a structure are stored as pseudoatoms. The points are calculated as the centers of mass from the cavity neighbor atoms (all atoms touched by a 1.4 Angstrom probe inside the cavity). 


3.6) Write cavity neighbors

Creates a PDB file, where the neighbor atoms of a cavity (all atoms touched by a 1.4 Angstrom probe inside the cavity) are listed. The residue number of the atoms is relabeled to contain the cavity number.


3.7) Defining which atoms are to be used for calculation

Voronoia provides several switches to include or exclude general subsets of atoms [with defaults]:
- surface atoms [off]
- buried atoms [on]
- cavity neighbors [on]
- non-cavity neighbors [on]

Additionally, specific regions can be defined by adding SELECT lines to the .vol files right after the header.
Lines starting with SELECT will determine the parts of the file that will be selected for calculation. The selections are additive. Examples;
       SELECT chain A AND resi 1-100
       SELECT name CA AND resn TYR
       SELECT resn HOH
       SELECT resi 65
When no SELECT commands are specified, all atoms are selected automatically according to the more general restrictions like surface, buried etc.



------------------------------------------------------------------------
4. Using the graphical Interface
------------------------------------------------------------------------

Most of the GUI elements for the functions in 3. can be used intuitively. Some points should be clarified here:

 * Prior to all operations, volume files need to be calculated. Alternatively, pre-calculated .vol files can be provided. Upon calculation with the 'Calculate structures' button, the name of the resulting file is used for the other functions automatically.
 
 * To calculate packing data for a whole dataset, choose a directory with .ent/.vol files. All of them will be processed subsequently.
 
 * The PyMOL tab only appears if Voronoia is installed as a PyMOL plugin.

 
------------------------------------------------------------------------
5. Command Line Reference
------------------------------------------------------------------------

Generally, Voronoia is called from a command line console with:

    C:\Python2.4\Python.exe Voronoia.py <command> [options] <input file> (Windows)
    or
    Voronoia.py <command> [options] <input file> (Linux),

where <command> is one of:
    calc         - calculate packing .vol file(s)
    report       - write concise packing report from .vol file(s)
    average      - calculate average packing densities from .vol file(s)
    compare      - compare .vol file(s) to average packing densities
    cavities     - write cavities from .vol file(s) as PDB file(s)
    cavneighbors - write cavity neighbor atoms from .vol file(s) as PDB file(s)

where [options] for 'calc' are none to many of:
    -g <float> - grid distance for 'calc' (default 0.1)
    -stouten   - use Stouten radii instead of ProtOr radii set
    -e <file>  - location of executable (get_volume.exe)
    
where [options] for other commands are none to many of:
    -d             - process all files in directory <input file>
    -o             - output directory
    -surface       - use surface atoms only
    -buried        - use non-surface atoms only
    -cavnb         - use cavity neighbor atoms only
    -no_cavnb      - discard atoms neighboring cavities     
    -h             - include hetero groups 
    -w             - include water     
    -r <file>      - specify the reference file for 'compare'
    -volume        - calculate atomic volumes instead of packing densities
    -protor        - use concise ProtOr atom typing scheme with 19 types
    -html          - HTML output of reports.


------------------------------------------------------------------------
6. Volume Files (.vol)
------------------------------------------------------------------------

Most parts of Voronoia work with a special PDB structure format, the .vol files. It corresponds mostly to the PDB format, with a few extensions:
- After the header, a LENGTH line gives a precise count of residues and atoms.

- The following NRHOLE line gives the numbers of internal cavities that are partially or not filled, before and after hetero atoms (including water???) are removed.
ANDREAN!?

- A number of HOLE NUMBER lines follow each containing a list of the atom indices for a particular cavity.

- The ATOM lines contain the follwing modified data: 
    col.61-67: b_factor substituted by total atomic volume
    col.68-73: Van-der-Waals atomic volume
    col.74-80: Solvent Excluded atomic volume
    col.82   : bit indicating solvent accessibility (0-surface, 1-buried)

        
-------------------------------------------------------------------------
7. Appendix
-------------------------------------------------------------------------

Technical Notes:

Voronoia consists of a core algorithm compiled from Delphi (in get_volume.exe) and a wrapper software written in Python. It uses the Pmw and PIL libraries. All front-ends (Windows/Linux/GUI/PyMOL/commandline) use the same modules. The executables for Windows have been compiled using Py2exe and distutils.

To see the logo, the PIL library is necessary (therefore we made it optional).

-------------------------------------------------------------------------

ProtOr atom radii:

The atom radii below are used depending on the chemical configuration of an atom, as published by Tsai et al [2]:



-------------------------------------------------------------------------

References:

[1]   1.Goede, A., Preissner, R., and Frömmel, C., Voronoi Cell: New Method for Allocation of Space among Atoms: Elimination of Avoidable Errors in Calculation of Atomic Volume and Density. 1997.

[2]   Tsai J, Taylor R, Chothia C, Gerstein M. The packing density in proteins: standard radii and volumes. J Mol Biol. 1999 Jul 2;290(1):253-66.
   PMID: 10388571

[3]   Hildebrand PW, Rother K, Goede A, Preissner R, Frommel C. Molecular Packing and Packing Defects in Helical Membrane Proteins. Biophys J. 2004 Nov 19; 
   PMID: 15556989 

[4]   Rother K, Preissner R, Goede A, Frommel C. Inhomogeneous molecular density: reference packing densities and distribution of cavities within proteins. Bioinformatics. 2003 Nov 1;19(16):2112-21.
      PMID: 14594717 

