import sys
"""
Dies ist Knut das Problemloeseprogramm.

Liest ein .vol file und macht ein normal aussehendes PDB draus.

Schnipselt aus allen zeilen mit ATOM das Ende nach dem B-Faktor ab.

Starten mit python knut.py <volfile> <pdbfile>
"""
in_file = sys.argv[1]
out_file = sys.argv[2]

out = []
for line in open(in_file):
    if line.startswith('ATOM'):
        line = line[:60]+line[61:67]+'      VORO'+line[22:26]+'\n'
        # evtl stattdessen
        # line = line[:60]+line[60:67]+'     VORO'+line[22:26]+'\n'

    out.append(line)

open(out_file,'w').writelines(out)

