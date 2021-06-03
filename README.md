# PDB_ligand_quality_composite_score
Data for multivariate data analysis in constructing PDB ligand quality composite score\
Column description:\
"pdb_ligand": PDB ID and CCD ID of the ligand instance. Multiple instanes of the same compound in the same structure are averaged in calculating the values of quality indicators\
"rsr": Real space R factor\
"rscc": Real space correlation coefficient\
"mogul_bonds_rmsz": Root-Mean-Squared deviation Z-score of all bond lengths\
"mogul_angles_rmsz": Root-Mean-Squared deviation Z-score of all bond angles\
"clash_per_atom": Inter-molecular clashes on a ligand structure divided by the number of non-hydrogen atoms of the ligand\
"formula_weight": Formula weight\
"refine.ls_d_res_high": High resolution limit of the PDB structure with the ligand\
"fit_pc1": 1st princinpal component of RSR and RSCC\
"geo_pc1": 1st princinpal component of mogul_bonds_rmsz and mogul_angles_rmsz
