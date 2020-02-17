import numpy as np
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
        self.chemical_composition = self.extract_composition()

    def __ignore_isotopes(self):
        '''remove the isotope tags in the end of the chemical formula
        '''
        for isotope_tag in self.ignore_isotopes:
            try:
                if self.chemical_formula.index(isotope_tag)>=self.ignore_isotopes[isotope_tag]:
                    self.chemical_formula = self.chemical_formula.replace(isotope_tag, '')
            except ValueError as e: # if not isotope_tag in the formular
                continue

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

    def get_stoich_electron_donor(self):
        a = self.chemical_composition['C']
        b = self.chemical_composition['H']
        c = self.chemical_composition['N']
        d = self.chemical_composition['O']
        e = self.chemical_composition['S']
        f = self.chemical_composition['P']
        z = 0

        return np.array([
            -1,
            -(3*a+4*e-d),
            a,
            c,
            e,
            f,
            5*a+b-4*c-2*d+7*e-f,
            -z+4*a+b-3*c-2*d+5*e-2*f,
            0,
            0
        ])

    def get_stoich_electron_acceptor(self):
        stoich_electron_acceptor = np.zeros(10)
        stoich_electron_acceptor[8] = -1  # oxygen
        stoich_electron_acceptor[6] = -4  #  h+
        stoich_electron_acceptor[7] = -4  #  e-
        stoich_electron_acceptor[1] = 2  #  h2o
        return stoich_electron_acceptor

    def get_stoich_catabolic_reaciton(self, stoich_electron_donor, stoich_electron_acceptor):
        yEd = stoich_electron_donor[7]
        yEa = stoich_electron_acceptor[7]
        return stoich_electron_donor-(yEd/yEa)*stoich_electron_acceptor

    def get_stoich_anabolic_reaction(self, 
                                     chemical_composition, 
                                     stoich_electron_donor, 
                                     stoich_electron_acceptor):
        # "C","H","N","O","S","P"
        a = self.chemical_composition['C']
        # b = self.chemical_composition['H']
        # c = self.chemical_composition['N']
        # d = self.chemical_composition['O']
        # e = self.chemical_composition['S']
        # f = self.chemical_composition['P']
        # z = 0

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

        for eA4Ana in eA4Anabolic:
            if eA4Ana == 'O2':
                stoichAnStar_O2 = stoichAnStarB+(1/a)*stoich_electron_donor
                yEana = stoichAnStar_O2[7]
                if yEana > 0:
                    yEa = stoich_electron_acceptor[7]
                    stoichAn_O2 = stoichAnStar_O2-yEana/yEa*stoich_electron_acceptor
                elif yEana < 0:
                    yEd = stoich_electron_donor[7]
                    stoichAn_O2 = stoichAnStar_O2-yEana/yEd*stoich_electron_donor
                else:
                    stoichAn_O2 = stoichAnStar_O2
            elif eA4Ana == 'HCO3-':
                yEd = stoich_electron_donor[7]
                yEa = stoichAnStarB[7]
                stoichAn_HCO3 = stoich_electron_donor-(yEd/yEa)*stoichAnStarB
                stoichAn_HCO3 = stoichAn_HCO3/stoichAn_HCO3[9]

        return stoichAn_O2, stoichAn_HCO3

    def get_lambda(self,
                   chemical_composition, 
                   stoich_electron_donor, 
                   stoich_cat_rxns, 
                   stoich_anabolic_O2, 
                   stoich_anabolic_HCO3):
        a = chemical_composition['C']
        b = chemical_composition['H']
        c = chemical_composition['N']
        d = chemical_composition['O']
        e = chemical_composition['S']
        f = chemical_composition['P']
        z = 0

        ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D 
        nosc = -ne/a+4  # nominal oxidataion state of carbon 
        delGcox0PerE = 60.3-28.5*nosc  # kJ/C-mol
        delGcox0 = delGcox0PerE*a*np.abs(stoich_electron_donor[0])  # kJ/rxn

        # - estimate delGf0 for electron donor
        delGf0_D_zero = 0
        delGf0_zero = np.array([delGf0_D_zero, -237.2, -586.8, -79.3, -1096.1, 12.1, 0, 0, 16.4, -67])
        delGcox0_zero = np.dot(delGf0_zero, stoich_electron_donor)
        delGf0_D_est = (delGcox0-delGcox0_zero)/stoich_electron_donor[0]
        # - finally, delGf0
        delGf0 = delGf0_zero
        delGf0[0] = delGf0_D_est

        # - standard delG at pH=0
        delGcat0 = np.dot(delGf0, stoich_cat_rxns)
        delGan0_O2 = np.dot(delGf0, stoich_anabolic_O2)
        delGan0_HCO3 = np.dot(delGf0, stoich_anabolic_HCO3)

        # - stadard delG at pH=7
        R = 0.008314  # kJ/(K.mol)
        T = 298  # K
        iProton = 6  # [eD,h2o,hco3-,nh4+,hpo4**2-,hs-,h+,e-,eA,biom]
        delGcox = delGcox0+R*T*stoich_electron_donor[iProton]*np.log(1e-7)
        delGcat = delGcat0+R*T*stoich_cat_rxns[iProton]*np.log(1e-7)
        delGan_O2 = delGan0_O2+R*T*stoich_anabolic_O2[iProton]*np.log(1e-7)
        delGan_HCO3 = delGan0_HCO3+R*T*stoich_anabolic_HCO3[iProton]*np.log(1e-7)

        # The Thermodynamic Electron Equivalents Model (TEEM)
        # --------
        eta = 0.43
        delGsyn = 200  # kJ/(mol.X)
        if delGan_O2 < 0:
            m_O2 = 1
        else:
            m_O2 = -1

        if delGan_HCO3 < 0:
            m_HCO3 = 1
        else:
            m_HCO3 = -1

        lambda_O2 = (delGan_O2*eta**m_O2+delGsyn)/(-delGcat*eta)
        lambda_HCO3 = (delGan_HCO3*eta**m_HCO3+delGsyn)/(-delGcat*eta)

        if lambda_O2 > 0:
            stoichMet_O2 = lambda_O2*stoich_cat_rxns+stoich_anabolic_O2
        else:
            stoichMet_O2 = stoich_anabolic_O2

        if lambda_HCO3 > 0:
            stoichMet_HCO3 = lambda_HCO3*stoich_cat_rxns+stoich_anabolic_HCO3
        else:
            stoichMet_HCO3 = stoich_anabolic_HCO3

        delGdis_O2 = np.dot(delGf0, stoichMet_O2) + R*T*stoichMet_O2[iProton]*np.log(1e-7)
        delGdis_HCO3 = np.dot(delGf0, stoichMet_HCO3) + R*T*stoichMet_HCO3[iProton]*np.log(1e-7)
        # delGdis = 200+18*(6-a)**1.8 + np.exp(((-0.2+nosc)**2)**0.16*(3.6+0.4*a))

        return \
            [lambda_O2,lambda_HCO3], \
            [delGcox0PerE,delGcox0,delGcox,delGcat0,delGcat,delGan0_O2,delGan0_HCO3,\
                delGan_O2,delGan_HCO3,delGdis_O2,delGdis_HCO3],\
            stoichMet_O2,\
            stoichMet_HCO3

    def get_all_thermo_stoich(self):
        #######################################################################
        # Step 1a) stoichiometries for an electron donor
        #######################################################################
        self.stoich_electron_donor = self.get_stoich_electron_donor()

        #######################################################################
        # Step 1b) stoichiometries for an electron acceptor (i.e., oxygen)
        #######################################################################
        self.stoich_electron_acceptor = self.get_stoich_electron_acceptor()
        
        #######################################################################
        # Step 1c) stoichCat: stoichiometries for catabolic reaciton
        #######################################################################
        self.stoich_cat_rxns = \
            self.get_stoich_catabolic_reaciton(self.stoich_electron_donor,
                                               self.stoich_electron_acceptor)

        #######################################################################
        # Step 2a) stoichAnStar: stoichiometries for anabolic reaciton 
        #          (N source = NH4+)
        #######################################################################
        self.stoich_anabolic_O2, self.stoich_anabolic_HCO3 = \
            self.get_stoich_anabolic_reaction(self.chemical_composition,
                                              self.stoich_electron_donor,
                                              self.stoich_electron_acceptor)

        # Step 3: get lambda
        # - estimate delGcox0 using LaRowe and Van Cappellen (2011)
        self.th_lambda, self.delta_gibbs_energy, self.stoich_metabolic_O2, self.stoich_metabolic_HCO3 = \
            self.get_lambda(self.chemical_composition,
                            self.stoich_electron_donor,
                            self.stoich_cat_rxns,
                            self.stoich_anabolic_O2, 
                            self.stoich_anabolic_HCO3)

        return self.delta_gibbs_energy + self.th_lambda + \
            list(self.stoich_electron_donor) + \
            list(self.stoich_electron_acceptor) + \
            list(self.stoich_cat_rxns) + \
            list(self.stoich_anabolic_O2) + \
            list(self.stoich_anabolic_HCO3) + \
            list(self.stoich_metabolic_O2) + \
            list(self.stoich_metabolic_HCO3)

