import warnings
warnings.filterwarnings('ignore')

from MDAnalysis.lib.distances import capped_distance
import numpy as np
import scipy.sparse
import scipy.stats

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['Cluster']


# suppress some MDAnalysis warnings when writing PDB/GRO files


class Cluster(AnalysisBase):
    """
    A class for analyzing lipid clustering in a bilayer system.
    
    This class identifies and analyzes clusters of lipid molecules based on their
    spatial proximity. It uses a distance-based criterion to determine which lipids
    belong to the same cluster.
    """
    
    def __init__(self, universe, residueGroup: dict, cutoff: float = None, file_path: str = None):
        """
        Initialize the Cluster analysis class.
        
        Parameters
        ----------
        universe : MDAnalysis.Universe
            The MDAnalysis Universe object containing the molecular dynamics trajectory.
            This should include both the structure file (e.g., .gro) and trajectory file (e.g., .xtc).
            
        residueGroup : dict
            A dictionary specifying the atoms to use for clustering analysis.
            Format: {'lipid_name': ['atom_names']}
            Example: {
                'DPPC': ['PO4'],
                'CHOL': ['ROH']
            }
            
        cutoff : float, optional
            The distance cutoff (in Å) used to determine if two lipids belong to the same cluster.
            If two lipids are closer than this distance, they are considered to be in the same cluster.
            If None, a default value will be used.
            
        file_path : str, optional
            The path where the analysis results will be saved as a CSV file.
            If None, results will not be saved to disk.
            
        Attributes
        ----------
        headAtoms : MDAnalysis.AtomGroup
            The selected atoms of all lipid molecules used for clustering.
            
        _n_residues : int
            The total number of lipid molecules in the system.
            
        resids : numpy.ndarray
            The residue IDs of all lipid molecules.
            
        resnames : numpy.ndarray
            The residue names of all lipid molecules.
            
        results.Cluster : numpy.ndarray
            A 2D array storing the cluster assignments for all lipid molecules across all frames.
            Shape: (n_residues, n_frames)
            
        cutoff : float
            The distance cutoff used for clustering.
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residueGroup)
        self.cutoff = cutoff
        self.file_path = file_path

        # Convert atom names to space-separated strings
        self.atomSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Atoms for clustering:", self.atomSp)

        # Initialize atom selection
        self.headAtoms = self.u.atoms[[]]

        # Select atoms for all specified lipid types
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (self.residues[i], self.atomSp[self.residues[i]]), updating=False)

        # Set basic attributes
        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Cluster = None

        # Record analysis parameters
        self.parameters = str(residueGroup) + 'Cutoff:' + str(self.cutoff)

    @property
    def Cluster(self):
        return self.results.Cluster

    @property
    def Resindices(self):
        return self.results.Resindices

    def _prepare(self):
        self.results.Cluster = np.zeros([self.n_frames])
        self.results.Resindices = np.full(self.n_frames, fill_value=0, dtype=object)

    def _single_frame(self):
        positions_single_frame = self.headAtoms.positions
        pairs = capped_distance(
            positions_single_frame
            , positions_single_frame
            , max_cutoff=self.cutoff
            , box=self._ts.dimensions
            , return_distances=False
        )

        ref, nei = np.unique(self.resids[pairs], axis=0).T

        index_needs = ref != nei

        ref = ref[index_needs]
        nei = nei[index_needs]

        data = np.ones_like(ref)

        neighbours_frame = scipy.sparse.csr_matrix(
            (data, (ref, nei))
            , dtype=np.int8
            , shape=(self._n_residues, self._n_residues)
        )

        _, com_labels = scipy.sparse.csgraph.connected_components(neighbours_frame)

        unique_com_labels, counts = np.unique(com_labels, return_counts=True)
        largest_label = unique_com_labels[np.argmax(counts)]
        self.results.Cluster[self._frame_index] = max(counts)
        self.results.Resindices[self._frame_index] = self.resnames[com_labels == largest_label]

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': [i for i in range(self.start, self.stop, self.step)]
                , 'results': self.results.Cluster
                , 'file_path': self.file_path
                , 'description': 'Cluster'
                , 'parameters': self.parameters
            }
            WriteExcelBubble(**dict_parameter).run()


def make_data(residues_group):
    ls1 = []
    for residue, atoms in residues_group.items():
        str_atom = ' '.join(atoms)
        ls1.append(f'resname {residue} and name {str_atom}')
    return ls1


if __name__ == "__main__":
    import MDAnalysis as mda

    gro_file = "../cases/lnb.gro"
    xtc_file = "../cases/md.xtc"
    csv_file = "../cases/csv/area_step5_lnb.csv"
    u = mda.Universe(gro_file, xtc_file)
    cls1 = Cluster(u, {'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}, cutoff=20, file_path='E:/untitled1.csv')
    cls1.run(0, 100)
