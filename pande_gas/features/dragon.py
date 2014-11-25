"""
Dragon descriptors.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import numpy as np

from pande_gas.features import Featurizer
from pande_gas.utils.dragon_utils import Dragon


class DragonDescriptors(Featurizer):
    """
    Calculate Dragon descriptors.

    Parameters
    ----------
    assign_stereo_from_3d : bool, optional (default False)
        Assign stereochemistry from 3D coordinates. This will overwrite any
        existing stereochemistry information on molecules.
    """
    name = 'dragon'

    def __init__(self, assign_stereo_from_3d=False):
        self.engine = Dragon(assign_stereo_from_3d=assign_stereo_from_3d)

    def featurize(self, mols, parallel=False, client_kwargs=None,
                  view_flags=None):
        """
        Calculate features for molecules.

        Parameters
        ----------
        mols : iterable
            RDKit Mol objects.
        parallel : bool, optional (default False)
            Train subtrainers in parallel using IPython.parallel.
        client_kwargs : dict, optional
            Keyword arguments for IPython.parallel Client.
        view_flags : dict, optional
            Flags for IPython.parallel LoadBalancedView.
        """
        if parallel:
            from IPython.parallel import Client

            if client_kwargs is None:
                client_kwargs = {}
            if view_flags is None:
                view_flags = {}
            client = Client(**client_kwargs)
            client.direct_view().use_dill()  # use dill
            view = client.load_balanced_view()
            view.set_flags(**view_flags)
            call = view.map(
                self._featurize,
                np.array_split(mols, len(client.direct_view())), block=False)
            features = call.get()
            features = np.vstack(features)

            # get output from engines
            call.display_outputs()

        else:
            features = self._featurize(mols)

        return np.asarray(features)

    def _featurize(self, mols):
        """
        Calculate Dragon descriptors.

        Parameters
        ----------
        mols : array_like
            Molecules.
        """
        return self.engine.get_descriptors(mols)
