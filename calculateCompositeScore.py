##
# File: calculateCompositeScore.py
# Date: 2020-11-29 Chenghua Shao
##
""" calcuate fitting and geometry composite scores for any ligand, referencing the entire PDB archive
The composite score is the ranking percentile of a ligand quality metrics of either experimental fitting or geometry.
The smaller the percentile, the worst the quality, 0% is the worst, 100% is the best.
There are two steps in calculating the percentile composite scores:
(1) Calculate the 1st Principal Component (PC1) of multivariate data:
(1a) Calculate PC1 of rsr and rscc, namely, fit_pc1 by the following formula:
((rsr-mean_rsr)/std_rsr)*loading_rsr + ((rscc-mean_rscc)/std_rscc)*loading_rscc
where mean_rsr is the mean of all rsr values in the reference data, and std_rsr is the standard deviation of all rsr values
loading_rsr is sqrt(2)/2.0, and loading_rscc is -sqrt(2)/2.0
(1b) Calculate PC1 of mogul_bonds_rmsz and mogul_angles_rmsz, namly, geo_pc1 by the following formula:
((mogul_bonds_rmsz-mean_mogul_bonds_rmsz)/std_mogul_bonds_rmsz)*loading_mogul_bonds_rmsz +
((mogul_angles_rmsz-mean_mogul_angles_rmsz)/std_mogul_angles_rmsz)*loading_mogul_angles_rmsz
where both loading_mogul_bonds_rmsz and loading_mogul_angles_rmsz is sqrt(2)/2.0
(2) Once the fit_pc1 and geo_pc1 is calculated, find the ranking percentile of the value in the PDB archive

input: reference PDB archive data, and rsr/rscc/mogul_bonds_rmsz/mogul_angles_rmsz values of any ligand
The reference PDB data have fit_pc1 and geo_pc1 pre-calculated for each averaged PDB ligand instance.
The ligands of the same type within one particular structure is combined to avoid overwelming impact on the ranking.
e.g. if there are 20 ATP ligands in PDB entry 1ABC, they are combined into one row of 1ABC-ATP
output: the composite ranking scores of fitting and geometry
"""

import csv
import numpy

class PCA:
    """ class to load reference PDB data, and then calculate PCA parameters for subsequent PCA analysis.
    This class should only loaded once and then used repeatly as parameters for calculation on multiple ligands    
    input: reference PDB archive data
    The reference PDB data have fit_pc1 and geo_pc1 pre-calculated for each PDB ligand instance.
    The ligands of the same type within one particular structure is combined to avoid overwelming impact on the ranking.
    e.g. if there are 20 ATP ligands in PDB entry 1ABC, they are combined into one row of 1ABC-ATP
    output: the references in organized format, d_ref, and the PCA parameters, d_par
    """
    def __init__(self):
        self.d_ref = {}
        self.d_par = {}

    def loadReference(self, filepath_ref = "Data/PCA_PDB_references.csv"): #load reference data
        self._readPcaReferences(filepath_ref) #read PCA reference data
        self._calculatePcaParameters() #calculate PCA parameters
        
    def _readPcaReferences(self,filepath, sep=","):
        ## load reference data as a dictionary, with each column of the reference as the key, and list of data as value
        ## e.g. self.d_ref["rsr"] is the list of rsr values for all archived PDB ligand instances 
        try:
            f = open(filepath)
            file = csv.reader(f, delimiter=sep)
        except IOError as msg: #minimal exception handling because the data is always packed in the proper format
            print(msg)
        else:
            l_header = next(file, None)
            for item in l_header[1:]:
                self.d_ref[item] = []
            for l_line in file:
                for (item, value) in zip(l_header[1:], l_line[1:]):
                    self.d_ref[item].append(float(value))
        finally:
            f.close()

    def _calculatePcaParameters(self):
        ## calculate PCA parameters to use, e.g. the mean, std and loading of each variable used in the formula
        for var in ["rsr","rscc","mogul_bonds_rmsz","mogul_angles_rmsz"]:
            self.d_par[var] = {}
            l_value = self.d_ref[var]
            self.d_par[var]["mean"] = numpy.mean(l_value)
            self.d_par[var]["std"] = numpy.std(l_value)
            self.d_par[var]["loading"] = numpy.sqrt(2)/2.0
        self.d_par["rscc"]["loading"] = -self.d_par[var]["loading"] #must change sign for rscc's loading, because it's opposite to rsr
            
    def _getRankingPercentile(self, score, l_ref):
        ## general function to calculate the rankign percentile of a new value in an existing list of values.
        l_ref.append(score) #add the new value to the list, this won't change the list self.d_ref.
        percentile = sorted(l_ref).index(score)/float(len(l_ref)-1) #using the ranking and then divided by length, -1 because index starts at 0
        return 1-percentile #use the opposite percentile so that the worst is the smallest, for PDB convention.

    def calculateCompositeScores(self, d_lig):
        """calculate the ranking composite scores of fitting and geometry of a ligand
        There are two steps in calculating the percentile composite scores:
        (1) Calculate the 1st Principal Component (PC1) of multivariate data:
        (1a) Calculate PC1 of rsr and rscc, namely, fit_pc1 by the following formula:
        ((rsr-mean_rsr)/std_rsr)*loading_rsr + ((rscc-mean_rscc)/std_rscc)*loading_rscc
        where mean_rsr is the mean of all rsr values in the reference data, and std_rsr is the standard deviation of all rsr values
        loading_rsr is sqrt(2)/2.0, and loading_rscc is -sqrt(2)/2.0
        (1b) Calculate PC1 of mogul_bonds_rmsz and mogul_angles_rmsz, namly, geo_pc1 by the following formula:
        ((mogul_bonds_rmsz-mean_mogul_bonds_rmsz)/std_mogul_bonds_rmsz)*loading_mogul_bonds_rmsz +
        ((mogul_angles_rmsz-mean_mogul_angles_rmsz)/std_mogul_angles_rmsz)*loading_mogul_angles_rmsz
        where both loading_mogul_bonds_rmsz and loading_mogul_angles_rmsz is sqrt(2)/2.0
        (2) Once the fit_pc1 is calculated, find the ranking percentile of the value in the PDB archive
        """
        try:
            if d_lig["num_atoms"]>0 and d_lig["num_modeled_atoms"]>0 and (d_lig["num_atoms"]-d_lig["num_modeled_atoms"])>1:
                ## adjustment for incomplete ligands with more than 1 non-hydrogen atoms missing
                incompleteness = float(d_lig["num_atoms"]-d_lig["num_modeled_atoms"])/float(d_lig["num_atoms"])
                rsr = float(d_lig["rsr"]) + 0.08235*incompleteness
                rscc = float(d_lig["rscc"]) - 0.09652*incompleteness
            else:
                rsr = float(d_lig["rsr"])
                rscc = float(d_lig["rscc"])
                
            fit_pc1 = ((rsr-self.d_par["rsr"]["mean"])/self.d_par["rsr"]["std"]) * \
                      self.d_par["rsr"]["loading"] + \
                      ((rscc-self.d_par["rscc"]["mean"])/self.d_par["rscc"]["std"]) * \
                      self.d_par["rscc"]["loading"]
            d_lig["fit_pc1"] = fit_pc1
            d_lig["fit_p"] = self._getRankingPercentile(fit_pc1, self.d_ref["fit_pc1"])
        except ValueError:
            d_lig["fit_pc1"] = "NA"
            d_lig["fit_p"] = "NA"
        except TypeError:
            d_lig["fit_pc1"] = "NA"
            d_lig["fit_p"] = "NA"

        try:
            mogul_bonds_rmsz = float(d_lig["mogul_bonds_rmsz"])
            mogul_angles_rmsz = float(d_lig["mogul_angles_rmsz"])


            geo_pc1 = ((mogul_bonds_rmsz-self.d_par["mogul_bonds_rmsz"]["mean"])/self.d_par["mogul_bonds_rmsz"]["std"]) * \
                      self.d_par["mogul_bonds_rmsz"]["loading"] + \
                      ((mogul_angles_rmsz-self.d_par["mogul_angles_rmsz"]["mean"])/self.d_par["mogul_angles_rmsz"]["std"]) * \
                      self.d_par["mogul_angles_rmsz"]["loading"]
            d_lig["geo_pc1"] = geo_pc1
            d_lig["geo_p"] = self._getRankingPercentile(geo_pc1, self.d_ref["geo_pc1"])
        except ValueError:
            d_lig["geo_pc1"] = "NA"
            d_lig["geo_p"] = "NA"
            
        return d_lig

def main():
    pca = PCA()
    pca.loadReference("PCA_PDB_references.csv")
    ## print(pca.d_par["rsr"])
    ## print(pca.d_par["rscc"])
    ## print(pca.d_par["mogul_bonds_rmsz"])
    ## print(pca.d_par["mogul_angles_rmsz"])

    ## use PDB entry 6WJC's Y01 ligand in Figure 5 as example to caculate
    d_lig = {}
    d_lig["rsr"] = 0.227 
    d_lig["rscc"] = 0.892
    d_lig["mogul_bonds_rmsz"] =  1.31
    d_lig["mogul_angles_rmsz"] = 1.34
    d_lig["num_atoms"] = 35 #number of non-hydrogen atoms based on chemical definition
    d_lig["num_modeled_atoms"] = 35 #number of modeled non-hydrogen atoms

    d_lig = pca.calculateCompositeScores(d_lig)
    print(d_lig["fit_p"]) #composite ranking score for fitting, as percentile of standing
    print(d_lig["geo_p"]) #composite ranking score for geometry, as percentile of standing
    
if __name__ == "__main__":
    main()

## 2020 archive parameters from the reference set
## rsr: {'std': 0.08816137336094877, 'loading': 0.7071067811865476, 'mean': 0.17678949511295347}
## rscc: {'std': 0.09031199482827988, 'loading': -0.7071067811865476, 'mean': 0.8937142765679803}
## mogul_bonds_rmsz: {'std': 1.1218352583745501, 'loading': 0.7071067811865476, 'mean': 1.1699254841050957}
## mogul_angles_rmsz: {'std': 1.0221227753881565, 'loading': 0.7071067811865476, 'mean': 1.1854950479577229}
    
                                               
                                                 
        
