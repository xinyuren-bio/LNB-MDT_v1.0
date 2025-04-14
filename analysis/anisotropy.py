import warnings
warnings.filterwarnings('ignore')

import pandas as pd
from scipy.linalg import eigh
import numpy as np

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['Anisotropy']


class Anisotropy(AnalysisBase):
    """
    A class for calculating the anisotropy of lipid molecules in a bilayer.
    
    This class analyzes the orientational order of lipid molecules by calculating
    the order parameter tensor and its eigenvalues.
    """

    def __init__(self, universe, residueGroup: dict, file_path: str = None):
        """
        Initialize the Anisotropy analysis class.
        
        Parameters
        ----------
        universe : MDAnalysis.Universe
            The MDAnalysis Universe object containing the molecular dynamics trajectory.
            This should include both the structure file (e.g., .gro) and trajectory file (e.g., .xtc).
            
        residueGroup : dict
            A dictionary specifying the head atoms for each lipid type.
            Format: {'lipid_name': ['head_atom_names']}
            Example: {
                'DPPC': ['PO4'],
                'CHOL': ['ROH']
            }
            
        file_path : str, optional
            The path where the analysis results will be saved as a CSV file.
            If None, results will not be saved to disk.
            
        Attributes
        ----------
        headAtoms : MDAnalysis.AtomGroup
            The selected head atoms of all lipid molecules.
            
        _n_residues : int
            The total number of lipid molecules in the system.
            
        resids : numpy.ndarray
            The residue IDs of all lipid molecules.
            
        resnames : numpy.ndarray
            The residue names of all lipid molecules.
            
        results.Anisotropy : numpy.ndarray
            A 2D array storing the anisotropy values for all lipid molecules across all frames.
            Shape: (n_residues, n_frames)
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residueGroup)
        self.file_path = file_path

        # Convert head atom names to space-separated strings
        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        # Initialize atom selection
        self.headAtoms = self.u.atoms[[]]

        # Select head atoms for all specified lipid types
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (self.residues[i], self.headSp[self.residues[i]]), updating=False)

        # Set basic attributes
        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Anisotropy = None

        # Record analysis parameters
        self.parameters = str(residueGroup)

    @property
    def Anisotropy(self):
        return self.results.Anisotropy

    def _prepare(self):
        self.results.Anisotropy = np.full([self.n_frames], fill_value=np.nan)

    def _single_frame(self):
        atomsPos = self.headAtoms.positions
        ovariance_matrix = np.cov(atomsPos.T)
        # 计算特征值
        eigenvalues, _ = eigh(ovariance_matrix)
        # print(q, eigenvalues)
        # 计算椭圆偏心率
        anisotropy = np.sqrt(1 - eigenvalues[0] ** 2 / eigenvalues[2] ** 2)
        self.results.Anisotropy[self._frame_index] = anisotropy

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': [i for i in range(self.start, self.stop, self.step)]
                , 'results': self.results.Anisotropy
                , 'file_path': self.file_path
                , 'description': 'Anisotropy'
                , 'parameters': self.parameters
            }
            WriteExcelBubble(**dict_parameter).run()


if __name__ == "__main__":

    import MDAnalysis as mda

    gro_file = "../cases/lnb.gro"
    xtc_file = "../cases/md.xtc"
    csv_file = "../cases/csv/area_step5_lnb.csv"
    u = mda.Universe(gro_file, xtc_file)
    cls1 = Anisotropy(u, {'DPPC':['PO4'], 'DUPC':['PO4'], 'CHOL':['ROH']}, file_path=csv_file)
    cls1.run(start=10, step=5, verbose=True)

