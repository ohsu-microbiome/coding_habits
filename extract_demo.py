#!/usr/bin/env python3

from Bio import AlignIO
import sys
import os
import argparse
from datetime import datetime

#--------------------------
#       argparse!
#--------------------------

parsing=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
    extract_demo.py  is used to extract primer amplicons or variable regions from 
    a single MSA of the ribosomal 16S gene sequence (rrn). Because the position of these features are based 
    on the E. coli genome, the coordinates refer to the position on the E. coli sequence. 
    Either CP016007.2543965.2545520 or CP016007.3589827.3591382 
    from the Silva database (v132) are considered close enough, but the first is favored. 
    
    The start and end coordinates on the reference sequence are passed through 
    the --excise option. The start coordinate should be smaller than the end 
    coordinate, but I'm using <list name>.sort() to force the order from low to high.
    
    Coordinates are as used by UGENE. The first nucleotide of the sequence that starts 
    the MSA is counted as 1. Gaps are not counted in individual sequences, but 
    are counted in the MSA. Coordinates values are inclusive.
             
        seq3 A T C - - C G T A -->
             1 2 3 . . 4 5 6 7 
             
        seq4 - - - C G - - T A -->
             . . . 1 2 . . 3 4
             
        MSA  1 2 3 4 5 6 7 8 9 -->

    Just pass the start and end values seperated by a comma:
    
        >extract_demo.py -i infile -e 2,5
        
    If seq3 is the E. coli reference, the MSA coordinates are determined and extracted
    
        seq3 T C - - C G
             2 3 . . 4 5
             
        seq4 - - C G - - 
             . . 1 2 . .
             
        MSA  2 3 4 5 6 7
        
    and written to a new multi-record fasta file, without the "-" characters.
    
        >seq3
        TCCG
        >seq4
        CG 
""", epilog="This script has no other dependencies besides BioPython")

parsing.add_argument("--infile", "-i", required=True, help="relative or absolute path to the MSA")
parsing.add_argument("--outfile", "-o", help="new filename, automatically written to ../processed_data")
parsing.add_argument("--excise", "-e", help="the start and stop coordinates on the reference E. coli sequence.")
parsing.add_argument("--verbose", "-v", action="store_true", help="spell out what is happening")
  
args=parsing.parse_args()

#--------------------------
#   get to work
#--------------------------

runtime=datetime.now().strftime('%Y-%m-%d_%H%M')

# *** validate data args ***

# check that the infile is a clustal msa
try:
    msa = AlignIO.read(args.infile, "clustal")
    if args.verbose:
        print("Format of infile is fine.")
except ValueError or AssertionError:
    print("\n\tThe format of the file \"{}\" doesn't meet the Clustal file specifications:\n\n\t1) The first line in the file must start with the words \"CLUSTAL W\" or \"CLUSTALW\". Other information in the first line is ignored.\n\t2) One or more empty lines.\n\t3) One or more blocks of sequence data. Each block consists of:\n\t\t* One line for each sequence in the alignment. Each line consists of:\n\t\t\t* the sequence name\n\t\t\t* white space\n\t\t\t* up to 60 sequence symbols.\n\t\t\t* optional - white space followed by a cumulative count of residues for the sequences\n\t\t* A line showing the degree of conservation for the columns of the alignment in this block\n\t\t(this can be left as blank. \" \" indicates no match, but for this application the conservation isn't important)\n\t\t* One or more empty lines.\n\n\tQuitting now.".format(args.infile))
    sys.exit(1)
    
# check that there's an E. coli rrn sequence CP016007.2543965.2545520 or CP016007.3589827.3591382

standards = ["CP016007.2543965.2545520", "CP016007.3589827.3591382"]
names=[record.id for record in msa]

if set(standards).intersection(set(names)):
    if standards[0] in names:
        map_reference=standards[0]
    else:
        map_reference=standards[1] 
    if args.verbose:
        print("Found at least one of the standards: {}".format(set(standards).intersection(set(names))))
else:
    print("\n\tThe MSA doesn't appear to have either of the E. coli rrn sequences needed for coordinates.\n\tCheck that the MSN includes one of these two sequences:\n\t\tCP016007.2543965.2545520\n\t\tCP016007.3589827.3591382\n\n\tQuitting now.")
    sys.exit(1)
    
# new folder to hold output
final_dir="processed_data/extracted_16s_regions_{}".format(datetime.now().strftime('%Y-%m-%d'))
if not os.path.exists(final_dir):
    os.makedirs(final_dir)
    if args.verbose==True:
        print("creating "+final_dir)

if args.outfile:
    outfile="{0}/extracted_{1}_{2}.fna".format(final_dir, args.outfile, runtime)
else:
    outfile="{0}/extracted_16s_{1}.fna".format(final_dir, runtime)
    
coords=[int(x) for x in args.excise.split(',')]
coords.sort()
# add 1 to the last value, AlignIO format is [inclusive:exclusive]
coords[-1]+=1

# make a dict out of the MSA   
msa_dict={m.id:m.seq for m in msa}

# list is coord of ref seq, index is position in MSA 
# probably the kind of thing that's already solved, don't care  
counter=0
ref_map=list()
for x in msa_dict[map_reference]:
    if x!="-":
        ref_map.append(counter)
        counter+=1
    else:
        ref_map.append(None)

# write as dictionary
extract=msa[:, ref_map.index(coords[0]):ref_map.index(coords[1])]

# I can't figure out how to do this directly through Biopython
# so I'm going to do it the dumb way:
# parse the fasta file to a dictionary while removing the gap characters
# then write the dictionary to a new fasta file

dumb={u.id:u.seq.ungap('-') for u in extract}
    
print("writing to {}".format(outfile))
with open(outfile, "w") as f:
    for k,v in dumb.items():
        blocks=[v[i:i+60] for i in range(0, len(v), 60)]
        f.write(">{}\n".format(k))
        for b in blocks:
            f.write(str(b)+"\n")






    
