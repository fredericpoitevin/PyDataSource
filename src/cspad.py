from __future__ import absolute_import
from . import PyDataSource

class Cspad(PyDataSource.Detector):
    """Cspad Detector Class.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)

        try:
            self.add.count('corr', name='corr_count')
        except:
            pass
#        self.add.projection('corr', 'r', name='corr_r')

        try:
            self._init_streak(**kwargs)
        except:
            pass

    def _init_streak(self, beam_x=None, beam_y=None, **kwargs):
        """
        Output cspad jet streak information 
        """
        # try getting from 'CXI:SC1:DIFFRACT:BEAM_X' 
        if not beam_x:
            beam_x = 2094.93
        if not beam_y:
            beam_y = -1796.57
        
        cy, cx = self.get_center(beam_x_pv, beam_y_pv)
        j_map_1, j_map_2 = find_proj_mapping(cy, cx)
        detector.add.parameter(proj_map_1 = j_map_1, proj_map_2 = j_map_2)
        detector.add.property(streak_angle)
        detector.add.property(streak_intensity)
        detector.add.property(streak_width)
        detector.add.property(streak_present)


    def add_max_plot(self):
        """
        Add max count plot of CsPad
        """
        if not self._det_config['stats'].get('corr_stats'):
            self.add.stats('corr')
            next(self)
        self.add.property(img_max)
        self.add.psplot('img_max')

    @property
    def streak_angle_raw(self):
        """
        Jet streak calculation
        Returns: jet angle, jet intensity (as standard deviations from the mean),
        jet width
        
        """
        import numpy as np
        from scipy.signal import peak_widths

        im1 = self.corr[1,-100:,:100]
        im2 = self.corr[17,-100:,:100]
        proj1 = np.zeros((100,80))
        proj2 = np.zeros((100,80))
        for a in range(-40,40):
            for i in range(im1.shape[0]):
                proj1[i,a+40]=im1[i,self.proj_map_1[i,a+40]]
                proj2[i,a+40]=im2[i,self.proj_map_2[i,a+40]]
        s = proj1.sum(axis=0)+proj2.sum(axis=0)
        s -= s.mean()
        s /= np.roll(s,10-s.argmax())[20:].std()
        peak = s[1:-1].argmax()+1
        try:
            peakwidth = peak_widths(s, [peak])[0][0]
        except Exception as e:
            peakwidth = 5
        return (np.pi*(peak-40)/360.0, s.max(), peakwidth)
    
    def get_center(self, x0_pv, y0_pv):
        center = (y0_pv.get(), x0_pv.get())
        cy, cx = get_center_coords(self, center)
        cy -= 185
        return cy, cx

    def to_pad_coord(self, point, i):
        '''Point: (y,x)'''
        import numpy as np

        pad = [1,9,17,25][i]
        origin = np.asarray((self.calibData.coords_x[pad,0,0], self.calibData.coords_y[pad,0,0]))
        unit_y = ((self.calibData.coords_x[pad,1,0]-self.calibData.coords_x[pad,0,0]),
                  (self.calibData.coords_y[pad,1,0]-self.calibData.coords_y[pad,0,0]))
        unit_x = ((self.calibData.coords_x[pad,0,1]-self.calibData.coords_x[pad,0,0]),
                  (self.calibData.coords_y[pad,0,1]-self.calibData.coords_y[pad,0,0]))
        matrix = np.asarray([[unit_y[0], unit_x[0]],[unit_y[1], unit_x[1]]])
        pos = np.linalg.solve(matrix, np.asarray(point)-origin)
        return pos

    def get_center_coords(self, center):
        import numpy as np

        cy = np.zeros(4)
        cx = np.zeros(4)
        for i in range(4):
            pos = to_pad_coord(self, center, i)
            cy[i], cx[i] = pos[0], pos[1]
        return cy, cx



def streak_angle(self):
    """
    Streak angle calculation
    """
    return self.streak_angle_raw[0]

def streak_intensity(self):
    """
    Streak intensity calculation
    """
    return self.streak_angle_raw[1]

def streak_width(self):
    """
    Streak width calculation
    """
    return self.streak_angle_raw[2]


def img_max(self):
    """
    Method to make image from corr_stats
    """
    try:
        import numpy as np
        stat = 'max'
        code = self.sourceData.eventCode
        attr = 'corr_stats'
        dat = getattr(self, attr).sel(stat=stat).sel(codes=code)[-1]
        xx = self.calibData.indexes_x
        yy = self.calibData.indexes_y
        a = np.zeros([int(xx.max())+1,int(yy.max())+1])
        a[xx,yy] = dat.data
        return a
    except:
        return None


def streak_present(self):
    im1 = self.asic[0]
    im2 = self.asic[2]
    return streak_present_im(im1) and streak_present_im(im2)

def streak_present_im(im):
    '''im is 2D np-array'''
    s = im[-10:].sum(axis=0)
    s -= s.mean()
    s /= np.roll(s,10-s.argmax())[20:].std()
    return s.max()>5

def find_proj_mapping(cy,cx):
    import numpy as np
    sq = 0

    j_index_1 = np.zeros((100,80), dtype=np.int64)
    j_index_2 = np.zeros((100,80), dtype=np.int64)
    for a in range(-40,40):
        ang = np.radians(float(a)/2)
        for i in range(100):
            j = int(np.tan(ang)*(100-i + cy[sq]) + cx[sq]) % 100
            j_index_1[i,a+40]=j
            j = int(np.tan(ang)*(100-i + cy[(sq+2)%4]) + cx[(sq+2)%4]) % 100
            j_index_2[i,a+40]=j
    return j_index_1, j_index_2



# methods from /reg/neh/operator/cxiopr/experiments/cxilr6716/src/to_hdf5.py

#def asic(self, attr='calib', 
#        asics=[1,9,17,25], ka=0,kb=194): 
#    """
#    Select inner asics
#    """
#    return getattr(self, attr)[asics,:,ka:kb]
#
#def streak_angle_raw(self):
#    """
#    Jet streak calculation
#    Returns: jet angle, jet intensity (as standard deviations from the mean),
#    jet width
#    
#    """
#    import numpy as np
#    from scipy.signal import peak_widths
#
#    sq = 0
#
#    asic = self.asic
#    im1 = asic[sq][-100:,:100]
#    im2 = asic[(sq+2)%4][-100:,:100]
#    proj1 = np.zeros((100,80))
#    proj2 = np.zeros((100,80))
#    for a in range(-40,40):
#        for i in range(im1.shape[0]):
#            proj1[i,a+40]=im1[i,self.proj_map_1[i,a+40]]
#            proj2[i,a+40]=im2[i,self.proj_map_2[i,a+40]]
#    s = proj1.sum(axis=0)+proj2.sum(axis=0)
#    s -= s.mean()
#    s /= np.roll(s,10-s.argmax())[20:].std()
#    
#    return (float(s.argmax()-40)/2, s.max(), peak_widths(s,[s.argmax()])[0][0]/2)
#
#def find_proj_mapping(cy,cx):
#    import numpy as np
#    sq = 0
#
#    j_index_1 = np.zeros((100,80), dtype=np.int64)
#    j_index_2 = np.zeros((100,80), dtype=np.int64)
#    for a in range(-40,40):
#        ang = np.radians(float(a)/2)
#        for i in range(100):
#            j = int(np.tan(ang)*(100-i + cy[sq]) + cx[sq]) % 100
#            j_index_1[i,a+40]=j
#            j = int(np.tan(ang)*(100-i + cy[(sq+2)%4]) + cx[(sq+2)%4]) % 100
#            j_index_2[i,a+40]=j
#    return j_index_1, j_index_2   
#

