#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

class AtomicLevel:
# Just default values. These get reassigned in element.add_level
    def __init__(self,i,n,E,g):
        self.list_index = i
        self.n_quantum = n 
        self.excitation_energy_ev = E
        self.statistical_weight = g
        self.ionization_frequency_threshold = 0.
        self.pop_fraction = 0.
        self.alpha = 0.









