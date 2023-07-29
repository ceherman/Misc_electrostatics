import numpy as np
import pandas as pd

import sys
sys.path.append("/home/chase/codes/python_functions/")
import akta

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

# Column info___________________________________________________________________
class column():
    def __init__(self, length, diameter, v_extra=0, sup_vel=300, name='some_column'):
        self.name       = name
        self.length     = length    # [cm]
        self.diameter   = diameter  # [cm]
        self.v_extra    = v_extra   # [ml]
        self.sup_vel    = sup_vel   # [cm h-1]
        self.get_area()             # [cm2]
        self.get_volume()           # [ml]
        self.get_flow()             # [ml min-1]
        self.get_cv_extra()         # column volumes
        return

    def get_area(self):
        self.area = np.pi/4.0 * self.diameter**2    # [cm2]
        return

    def get_volume(self):
        self.volume = self.area * self.length       # [ml]
        return

    def get_flow(self):
        self.flow = round(self.sup_vel/60.0 * self.area, 2) # [ml min-1]
        return

    def get_cv_extra(self):
        self.cv_extra = self.v_extra/self.volume
        return

# Run data (read)_______________________________________________________________

class run_data():
    def __init__(self, folder, file, buffer_ion_str=0.0):
        self.folder         = folder
        self.file           = file
        self.df             = pd.read_csv(self.folder + self.file)
        self.buffer_ion_str = buffer_ion_str
        self.column         = column(10.19, 0.5, 0.327)
        self.get_file_info()
        self.get_ion_str()
        self.clean_data()
        self.get_first_moment()     # [CV]
        return

    def get_file_info(self):
        """Extract info from the file name"""
        self.protein    = 'lys'
        self.resin      = 'sep'
        self.pH         = float(self.file[3])
        self.nacl       = float(self.file[5:find_nth(self.file, '_', 3)]) # [mM]
        self.replicate  = int(self.file[-1])
        return

    def get_ion_str(self):
        """Estimate ionic strength in the solution"""
        self.ion_str = self.buffer_ion_str + self.nacl
        return

    def clean_data(self):
        self.df['vol_ml'] = self.df.time_min * self.column.flow
        self.df['cv'] = self.df.vol_ml/self.column.volume # convert to CV
        self.df.absorb_mau -= self.df.at[0, 'absorb_mau'] # zero absorbance
        self.df['absorb_corrected'] = akta.correct_baseline_drift(self.df.vol_ml,
                                        self.df.absorb_mau) # correct for baseline drift
        return

    def get_first_moment(self):
        moment = akta.get_first_moment(self.df.cv, self.df.absorb_corrected)
        self.cv_moment = moment - self.column.cv_extra
        return

    def get_kprime(self):
        """Compute kprime"""
        self.kprime = (self.cv_moment - self.eps_t)/self.eps_t
        return

    def get_keq(self):
        """Compute dimensionless Keq"""
        self.keq = (self.cv_moment - self.eps_t)/(1.0 - self.eps_t)
        return






















#
