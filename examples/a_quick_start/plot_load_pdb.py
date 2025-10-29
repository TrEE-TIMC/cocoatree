
"""
=========================
Load a PDB structure file
=========================

This small example shows how to use :func:`cocoatree.io.load_pdb` to import
your own PDB structure file.
"""

# %%
# Import necessary packages
from cocoatree.io import load_pdb

# Provide path for thumbnail image
# sphinx_gallery_thumbnail_path = '../../doc/images/3TGI.png'

# %%
# We will use the PDB structure of rat trypsin as an example. The file can be
# downloaded at: https://www.rcsb.org/structure/3TGI. It is also the structure
# that is included in :func:`cocoatree.datasets.load_S1A_serine_proteases`.

# %%
pdb_seq, pdb_pos = load_pdb('data/3TGI.pdb', pdb_id='3TGI', chain='E')

# %%
# The ``pdb_id`` argument is the ID that will be used for the structure, in
# this case, we choose `3TGI`, but it could be `TRYPSIN` or any name that the
# user finds fitting.
#
# In order to understand the function's ``chain`` argument, you need to open
# the PDB in an editor to check the information included within your file:

"""HEADER    COMPLEX (SERINE PROTEASE/INHIBITOR)     15-JUL-98   3TGI              
TITLE     WILD-TYPE RAT ANIONIC TRYPSIN COMPLEXED WITH BOVINE                   
TITLE    2 PANCREATIC TRYPSIN INHIBITOR (BPTI)                                  
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: TRYPSIN;                                                   
COMPND   3 CHAIN: E;                                                            
COMPND   4 EC: 3.4.21.4;                                                        
COMPND   5 MOL_ID: 2;                                                           
COMPND   6 MOLECULE: BOVINE PANCREATIC TRYPSIN INHIBITOR;                       
COMPND   7 CHAIN: I;                                                            
COMPND   8 SYNONYM: BPTI                                                        
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: RATTUS NORVEGICUS;                              
SOURCE   3 ORGANISM_COMMON: NORWAY RAT;                                         
SOURCE   4 ORGANISM_TAXID: 10116;                                               
SOURCE   5 ORGAN: PANCREATIC;                                                   
SOURCE   6 MOL_ID: 2;                                                           
SOURCE   7 ORGANISM_SCIENTIFIC: BOS TAURUS;                                     
SOURCE   8 ORGANISM_COMMON: CATTLE;                                             
SOURCE   9 ORGANISM_TAXID: 9913                          """

# %%
# As you can see, there are actually two molecules in the PDB:
#   - rat trypsin
#   - bovine pancreatic trypsin inhibitor
#
# It is necessary to specify which molecule you wish to load by using the
# ``chain`` argument. In this case, it is ``chain='E'`` for trypsin, and
# ``chain='I'`` for the trypsin inhibitor.
#
# For a comparison, here are the first lines of E. coli dihydrofolate
# reductase's PDB file (which is the one included in cocoatree's
# :func:`cocoatree.datasets.load_DHFR`):

"""HEADER    OXIDOREDUCTASE                          02-FEB-11   3QL0              
TITLE     CRYSTAL STRUCTURE OF N23PP/S148A MUTANT OF E. COLI DIHYDROFOLATE      
TITLE    2 REDUCTASE                                                            
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: DIHYDROFOLATE REDUCTASE;                                   
COMPND   3 CHAIN: A;                                                            
COMPND   4 EC: 1.5.1.3;                                                         
COMPND   5 ENGINEERED: YES;                                                     
COMPND   6 MUTATION: YES                                                        
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: ESCHERICHIA COLI;                               
SOURCE   3 ORGANISM_TAXID: 364106;                                              
SOURCE   4 STRAIN: UTI89 / UPEC;                                                
SOURCE   5 GENE: FOLA, UTI89_C0054;                                             
SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 
SOURCE   7 EXPRESSION_SYSTEM_TAXID: 562                                         """

# %%
# In this case, there is only one molecule, which is specified by
# ``chain='A'``.
#
# You can now access the pdb sequence and the residue numbering:
print(pdb_seq)
print(pdb_pos)
