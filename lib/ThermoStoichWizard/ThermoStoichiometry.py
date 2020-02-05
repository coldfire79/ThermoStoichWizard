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
        chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]  # C H_1.8 N_0.2 O_0.5
        aB = chemFormBiom[0]
        bB = chemFormBiom[1]
        cB = chemFormBiom[2]
        dB = chemFormBiom[3]
        eB = chemFormBiom[4]
        fB = chemFormBiom[5]
        zB = chemFormBiom[6]

        ySource = -1
        yH2o = -(3*aB+4*eB-dB)
        yHco3 = aB
        yNh4 = cB
        yHpo4 = eB
        yHs = fB
        yH = 5*aB+bB-4*cB-2*dB+7*eB-fB
        yE = -zB+4*aB+bB-3*cB-2*dB+5*eB-2*fB
        # add additional components: e-acceptor and biomass in the end
        stoichAnStarB = np.array([ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE,0,0])
        stoichAnStarB = -stoichAnStarB
        stoichAnStarB[-1] = stoichAnStarB[0]
        stoichAnStarB[0] = 0

        # Step 2b) "overall" anabolic reaction
        eA4Anabolic = [ # electron acceptor for anabolic reaction
            'O2',    # Kleerebezem and Van Loosdrecht (2010)
            'HCO3-' # % McCarty (year?)
        ]

        for i in eA4Anabolic:
            eA4Ana = eA4Anabolic[i]
            if eA4Ana == 'O2':
                stoichAnStar_O2 = stoichAnStarB+(1/a)*stoich_electron_donor
                yEana = stoichAnStar_O2[7]
                if yEana > 0:
                    stoichAn_O2 = stoichAnStar_O2-yEana/yEa*stoich_electron_acceptor
                elif yEana < 0:
                    stoichAn_O2 = stoichAnStar_O2-yEana/yEd*stoich_electron_donor
                else:
                    stoichAn_O2 = stoichAnStar_O2
            elif eA4Ana == 'HCO3-':
                yEd = stoich_electron_donor[7]
                yEa = stoichAnStarB[7]
                stoichAn_HCO3 = stoich_electron_donor-(yEd/yEa)*stoichAnStarB
                stoichAn_HCO3 = stoichAn_HCO3/stoichAn_HCO3[9]

        # Step 3: get lambda
  
        # - estimate delGcox0 using LaRowe and Van Cappellen (2011)
        ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D 
        nosc = -ne/a+4  # nominal oxidataion state of carbon 
        delGcox0PerE = 60.3-28.5*nosc  # kJ/C-mol
        delGcox0 = delGcox0PerE*a*abs(stoich_electron_donor[0])  # kJ/rxn

        # - estimate delGf0 for electron donor
        delGf0_D_zero = 0
        delGf0_zero = [delGf0_D_zero, -237.2, -586.8, -79.3, -1096.1, 12.1, 0, 0, 16.4, -67]
        # delGcox0_zero = drop(delGf0_zero %*% stoichD)
        # delGf0_D_est = (delGcox0-delGcox0_zero)/stoichD[1]
        # # - finally, delGf0
        # delGf0 = delGf0_zero
        # delGf0[1] = delGf0_D_est

        # # - standard delG at pH=0
        # delGcat0 = drop(delGf0 %*% stoichCat)
        # delGan0_O2 = drop(delGf0 %*% stoichAn_O2)
        # delGan0_HCO3 = drop(delGf0 %*% stoichAn_HCO3)