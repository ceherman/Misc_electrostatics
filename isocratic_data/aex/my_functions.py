import numpy as np
import pandas as pd

import sys
sys.path.append("/home/chase/codes/python_functions/")
import akta

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

# Hard-coded data_______________________________________________________________
columns = {}
columns['phq']   = column(4.2, 0.5, 0.335)
columns['pxq_1'] = column(4.3, 0.5, 0.335)
columns['pxq_2'] = column(4.3, 0.5, 0.335, sup_vel=150)
columns['caq_1'] = column(5.1, 0.5, 0.335)
columns['caq_2'] = column(5.35, 0.5, 0.335)
columns['caq_3'] = column(3.7, 1.0, 1.085)
for name in columns.keys():
    columns[name].name = name


# Run data (read)_______________________________________________________________
porosities = pd.read_csv('./data/porosities.csv')

class run_data():
    def __init__(self, folder, file, buffer_ion_str=0.0):
        self.folder         = folder
        self.file           = file
        self.df             = pd.read_csv(self.folder + self.file)
        self.buffer_ion_str = buffer_ion_str
        self.get_file_info()
        self.get_ion_str()
        self.get_column()
        self.clean_data()
        self.get_first_moment()     # [CV]
        return

    def get_file_info(self):
        """Extract info from the file name"""
        self.protein    = self.file[:3]
        self.resin      = self.file[4:7]
        self.nacl       = float(self.file[8:8+self.file[8:].find('_')]) # [mM]
        self.replicate  = int(self.file[-1])
        return

    def get_ion_str(self):
        """Estimate ionic strength in the solution"""
        self.ion_str = self.buffer_ion_str + self.nacl
        return

    def get_column(self):
        if self.resin == 'phq':
            self.column = columns['phq']
        elif self.resin == 'pxq':
            if self.protein == 'adh' and self.nacl == 50.0 and self.replicate == 2:
                self.column = columns['pxq_2']
            else:
                self.column = columns['pxq_1']
        else:
            if self.protein == 'blg' or self.protein == 'cat':
                self.column = columns['caq_1']
            elif self.protein == 'adh':
                self.column = columns['caq_3']
            else:
                if self.nacl in [150.0, 100.0]:
                    self.column = columns['caq_2']
                else:
                    self.column = columns['caq_1']
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

    def get_eps_t(self):
        """Read total porosity (from the porosities file)"""
        self.eps_t = porosities.loc[porosities[(porosities.protein==self.protein) &\
                                (porosities.resin==self.resin)].index, 'eps_t'].iloc[0]
        return

    def get_phi(self):
        """Read the dimensionless phase ratio (from the porosities file)"""
        self.phi = porosities.loc[porosities[(porosities.protein==self.protein) &\
                                (porosities.resin==self.resin)].index, 'phi'].iloc[0]
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
