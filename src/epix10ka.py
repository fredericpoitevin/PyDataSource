from __future__ import absolute_import
from . import PyDataSource

class Epix10ka(PyDataSource.Detector):
    """Epix10ka Detector Class.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)

    @property
    def size(self):
        """
        size of data array
        """
        return self.configData.numberOfRows*self.configData.numberOfColumns
    
    @property
    def shape(self):
        """
        shape of data array
        """
        return (self.configData.numberOfRows, self.configData.numberOfColumns)

    @property
    def calibData(self):
        """
        Not yet implemented
        """
        return None

