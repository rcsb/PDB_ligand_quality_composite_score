# PDB_ligand_quality_composite_score
(1) Description of "PCA_PDB_references.csv"\
Reference data for multivariate data analysis in constructing PDB ligand quality composite ranking score\
Column description:\
"pdb_ligand": PDB ID and CCD ID of the ligand instance. Multiple instances of the same compound in the same structure are averaged in calculating the values of quality indicators\
"rsr": Real space R factor\
"rscc": Real space correlation coefficient\
"mogul_bonds_rmsz": Root-Mean-Squared deviation Z-score of all bond lengths\
"mogul_angles_rmsz": Root-Mean-Squared deviation Z-score of all bond angles\
"clash_per_atom": Inter-molecular clashes on a ligand structure divided by the number of non-hydrogen atoms of the ligand\
"formula_weight": Formula weight\
"refine.ls_d_res_high": High resolution limit of the PDB structure with the ligand\
"fit_pc1": 1st principal component of RSR and RSCC\
"geo_pc1": 1st principal component of mogul_bonds_rmsz and mogul_angles_rmsz\
Usage: This reference data set sets up a scale of quality scores from worst to best for fit_pc1 and geo_pc1. Any particular ligand's fit_pc1 and geo_pc1 scores are to be compared with this reference to determin its quality standing percentile, i.e. composite ranking score.\
\
(2) Description of "PDB_ligands_raw_part_#.csv"\
The complete raw data (all ligand validation data of PDB archive as of Sept 2020). It is too big for GitHub, so it is split into three pieces.\
Column description: (in addition to the above)\
"pdb_id": PDB ID\
"ligand_id": Ligand CCD ID\
"model": Model number\
"chain": Author's chain ID\
"resnum": Authhor's residue number\
"icode": Insertion code is applicable, empty if none\
"altcode" : Alternate conformer ID if applicable, empty if none\
"non_Hydrogen_atoms": Number of non-hydrogen atoms based on chemical definition\
"NatomsEDS": Number of modeled atoms for EDS calculation\
"clash": Number of intermolecular clashes\
"chiral-outlier": Number of stereochemical outliers\
\
(3) Description of "calculateCompositeScore.py"\
The python script to calculate composite ranking score of ligand quality. The script takes input of one ligand instance’s quality data (under “main” function), with an example embedded in it. User can modify the script by provide the validation data for a particular ligand instance. The ligand validation data for input can be obtained from any PDB structure’s publicly available PDF or XML formatted validation report. 
