"""
IPython.parallel utilities.
"""
import subprocess
import time
import uuid


class LocalCluster(object):
    """
    Start an IPython.parallel cluster on localhost.

    Parameters
    ----------
    n_engines : int
        Number of engines to initialize.
    """
    def __init__(self, n_engines):
        self.n_engines = n_engines

        # placeholders
        self.cluster_id = None
        self.controller = None
        self.engines = []

        # initialize the cluster
        self.start()

    def __del__(self):
        """
        Shut down the cluster.
        """
        self.stop()

    def start(self):
        """
        Start the cluster by running ipcontroller and ipengine.
        """
        self.cluster_id = uuid.uuid4()
        self.controller = subprocess.Popen(
            ['ipcontroller', '--ip=*',
             '--cluster-id={}'.format(self.cluster_id)])
        time.sleep(1)  # wait for controller to initialize
        for i in xrange(self.n_engines):
            engine = subprocess.Popen(
                ['ipengine', '--cluster-id={}'.format(self.cluster_id)])
            self.engines.append(engine)
        time.sleep(10)  # wait for engines to initialize

    def stop(self):
        """
        Shut down the cluster.
        """
        for engine in self.engines:
            engine.terminate()
        self.controller.terminate()
