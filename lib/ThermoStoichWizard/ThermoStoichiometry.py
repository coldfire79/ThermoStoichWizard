import re

class ThermoStoichiometry(object):
    """extract thermo stoichiometry from a chemical formula"""
    
    CHEMICAL_ELEMENTS = ["C","H","N","O","S","P"]

    def __init__(self, chemical_formula, ignore_isotopes={'C13':1}):
        '''
            chemical_formula: chemical formula (e.g. C27H15O6NC13SP)
            ignore_isotopes: ignore if you have this in the '>index' of the 
                chemical_formula. it should follow the name convention
        '''
        super(ThermoStoichiometry, self).__init__()
        self.chemical_formula = chemical_formula.upper()
        if ignore_isotopes is not None:
            self.ignore_isotopes = ignore_isotopes
            self.__ignore_isotopes()

    def __ignore_isotopes(self):
        '''remove the isotope tags in the end of the chemical formula
        '''
        for isotope_tag in self.ignore_isotopes:
            if self.chemical_formula.index(isotope_tag)>=self.ignore_isotopes[isotope_tag]:
                self.chemical_formula = self.chemical_formula.replace(isotope_tag, '')

    def extract_composition(self):
        '''extract the chemimcal composition from the chemical formula
        '''
        # search number of elements
        chem_comp = {}
        for e in self.CHEMICAL_ELEMENTS:
            composition = re.search('{}(\d*)'.format(e), self.chemical_formula)

            if composition:
                n_element = composition.group(1)
                if n_element=='': chem_comp[e] = 1
                else: chem_comp[e] = int(n_element)
            else: chem_comp[e] = 0
        return chem_comp
                
    def get_thermo_stoich(self):
        chem_comp = self.extract_composition()
        # "C","H","N","O","S","P"
        a = chem_comp['C']
        b = chem_comp['H']
        c = chem_comp['N']
        d = chem_comp['O']
        e = chem_comp['S']
        f = chem_comp['P']
        z = 0
        #######################################################################
        # Step 1a) stoichiometries for an electron donor
        #######################################################################
        stoich_electron_donor = np.zeros(10)
        stoich_electron_donor[0] = -1  # ySource
        stoich_electron_donor[1] = -(3*a+4*e-d)  # yH2o
        stoich_electron_donor[2] = a  # yHco3
        stoich_electron_donor[3] = c  # yNh4
        stoich_electron_donor[4] = e  # yHpo4
        stoich_electron_donor[5] = f  # yHs
        stoich_electron_donor[6] = 5*a+b-4*c-2*d+7*e-f  # yH
        stoich_electron_donor[7] = -z+4*a+b-3*c-2*d+5*e-2*f  # yE
        stoich_electron_donor[8] = 0 # add additional components: e-acceptor
        stoich_electron_donor[9] = 0 # add additional components: biomass

        #######################################################################
        # Step 1b) stoichiometries for an electron acceptor (i.e., oxygen)
        #######################################################################
        stoich_electron_acceptor = np.zeros(10)
        stoich_electron_acceptor[8] = -1  # oxygen
        stoich_electron_acceptor[6] = -4  #  h+
        stoich_electron_acceptor[7] = -4  #  e-
        stoich_electron_acceptor[1] = 2  #  h2o

        #######################################################################
        # Step 1c) stoichCat: stoichiometries for catabolic reaciton
        #######################################################################
        yEd = stoich_electron_donor[7]
        yEa = stoich_electron_acceptor[7]
        stoich_cat_rxns = stoich_electron_donor-(yEd/yEa)*stoich_electron_acceptor

        #######################################################################
        # Step 2a) stoichAnStar: stoichiometries for anabolic reaciton 
        #          (N source = NH4+)
        #######################################################################
        yEd = stoich_electron_donor[7]
        yEa = stoich_electron_acceptor[7]
        stoich_cat_rxns = stoich_electron_donor-(yEd/yEa)*stoich_electron_acceptor
        