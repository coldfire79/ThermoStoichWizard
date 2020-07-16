'''
Compute the thermodynamic stoichiometries for chemical compositions
'''

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import re
import os

# CHNOPS chemical elements
CHEMICAL_ELEMENTS = ["C","H","N","O","P","S"]
# TODO: how to use Candidates
REQUIRED_COLUMNS = CHEMICAL_ELEMENTS#+['Candidates']

class FTICRResult(object):
    """FTICR Result"""
    def __init__(self, tbl, dtype=np.int):
        super(FTICRResult, self).__init__()
        if self.isvalid(tbl):
            # drop the peaks with the same compositions
            self.tbl = tbl.drop_duplicates(subset=CHEMICAL_ELEMENTS+['Na','C13'])
            self._assigned_tbl = self._filter(tbl, dtype=dtype)

            # mapping table: cpd id and molecular formula (unique)
            self.id2mf = self._assigned_tbl.mf.to_dict()
            self.mf2id = pd.Series(self._assigned_tbl.index.values, index=self._assigned_tbl.mf.values).to_dict()
            
            self._num_peaks = tbl.shape[0]
            self._num_cpds = self._assigned_tbl.shape[0]
        else:
            print('[Error] Input table requires these columns (output format in Formularity)\n{}'
                .format(REQUIRED_COLUMNS))
            raise(Exception('Input table format'))
        
    @property
    def num_peaks(self):
        return self._num_peaks

    @property
    def num_cpds(self):
        return self._num_cpds

    def isvalid(self, tbl):
        '''
            validate if the input table contains the essential columns, 
            which aligns the output format of "Formularity".
        '''
        isvalid = np.sum([c not in tbl.columns for c in REQUIRED_COLUMNS])==0
        return isvalid

    def _filter(self, tbl, dtype=np.int):
        '''filter out unassigned peaks and assign formulas
            TODO: how to deal with C13 and Na
            TODO: how to deal with the duplicated mf
        '''
        # assign formulas
        def assign_formula(row):
            mf = ''
            for ele in CHEMICAL_ELEMENTS:
                if (ele in row)&(row[ele]>0):
                    if row[ele]==1: mf += ele
                    else: mf += ele+str(row[ele])
            return mf
        tbl[CHEMICAL_ELEMENTS] = tbl[CHEMICAL_ELEMENTS].astype(dtype)
        tbl['mf'] = tbl.apply(assign_formula, axis=1)
        tbl['cpd_id'] = ['xcpd__{}'.format(i) for i in range(tbl.shape[0])]
        tbl = tbl.set_index('cpd_id')
        
        # filter out unassigned peaks
        filter_condition = tbl[CHEMICAL_ELEMENTS].sum(axis=1)>0
        if 'C13' in tbl.columns:
            tbl['C13'] = tbl['C13'].astype(np.int)
            filter_condition &= tbl['C13']==0
        if 'Na' in tbl.columns:
            tbl['Na'] = tbl['Na'].astype(np.int)
            filter_condition &= tbl['Na']==0
        return tbl[filter_condition]

    def to_csv(self, fout):
        self._assigned_tbl.to_csv(fout)

    def _batch_stoichiometries(self):
        '''extract all stoichiometries
        '''
        def get_stoichiometry(row):
            chem_comp = {e:row[e] for e in CHEMICAL_ELEMENTS}
            thermo_stoich = ThermoStoichiometry(chem_comp)
            thermo_stoich.get_all_thermo_stoich()
            return thermo_stoich
        return self._assigned_tbl.apply(get_stoichiometry, axis=1).to_dict()

    def run(self):
        self.all_stoich = self._batch_stoichiometries()
        # "stoichD","stoichA","stoichCat","stoichAn_O2","stoichAn_HCO3","stoichMet_O2","stoichMet_HCO3"
        stoich_colnames = ["donor","h2o","hco3","nh4","hpo4","hs","h","e","acceptor","biom"]
        self.stoichD = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_electron_donor for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        self.stoichA = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_electron_acceptor for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        self.stoichCat = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_cat_rxns for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        self.stoichAn_O2 = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_anabolic_O2 for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        self.stoichAn_HCO3 = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_anabolic_HCO3 for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        self.stoichMet_O2 = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_metabolic_O2 for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        self.stoichMet_HCO3 = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].stoich_metabolic_HCO3 for cpd in self.all_stoich},
            orient='index', columns=stoich_colnames)
        
        thermo_colnames = ["delGcox0PerC","delGcox0","delGcox","delGcat0","delGcat","delGan0_O2","delGan0_HCO3",
                "delGan_O2","delGan_HCO3","delGdis_O2","delGdis_HCO3","lambda_O2","lambda_HCO3"]
        self.thermo = pd.DataFrame.from_dict({self.id2mf[cpd]:self.all_stoich[cpd].delta_gibbs_energy+self.all_stoich[cpd].th_lambda for cpd in self.all_stoich},
            orient='index', columns=thermo_colnames)
        
    def save_result_files(self, folder):
        # save to csv files
        # self.stoichD.to_csv(folder+'/stoichD.csv')
        # self.stoichA.to_csv(folder+'/stoichA.csv')
        # self.stoichCat.to_csv(folder+'/stoichCat.csv')
        # self.stoichAn_O2.to_csv(folder+'/stoichAn_O2.csv')
        # self.stoichAn_HCO3.to_csv(folder+'/stoichAn_HCO3.csv')
        self.stoichMet_O2.to_csv(folder+'/stoichMet_O2.csv')
        # self.stoichMet_HCO3.to_csv(folder+'/stoichMet_HCO3.csv')
        self.thermo.to_csv(folder+'/thermodynamic_props.csv')
    
    def create_fba_model_files(self, folder, prefix='temp'):
        compounds_file = os.path.join(folder, "{}_comps.tsv".format(prefix))
        self.create_cpd_file_fba_model(compounds_file)
        # self.create_rxn_file_fba_model(self.stoichD, os.path.join(folder, "temp_stoichD.tsv"))
        # self.create_rxn_file_fba_model(self.stoichA, os.path.join(folder, "temp_stoichA.tsv"))
        # self.create_rxn_file_fba_model(self.stoichCat, os.path.join(folder, "temp_stoichCat.tsv"))
        # self.create_rxn_file_fba_model(self.stoichAn_O2, os.path.join(folder, "temp_stoichAn_O2.tsv"))
        # self.create_rxn_file_fba_model(self.stoichAn_HCO3, os.path.join(folder, "temp_stoichAn_HCO3.tsv"))
        self.create_rxn_file_fba_model(self.stoichMet_O2, os.path.join(folder, "{}_stoichMet_O2.tsv".format(prefix)))
        # self.create_rxn_file_fba_model(self.stoichMet_HCO3, os.path.join(folder, "temp_stoichMet_HCO3.tsv"))

    def create_cpd_file_fba_model(self, fout):
        comp_cols = ['id','name','formula','charge','inchikey','smiles','deltag','kegg id','ms id']
        compounds = [{'id':'{}_c0'.format(cid),'formula':self.id2mf[cid]} for cid in self.id2mf]
        compounds.append({'id':'h2o_c0','formula':'H2O'})
        compounds.append({'id':'hco3_c0','formula':'HCO3'})
        compounds.append({'id':'nh4_c0','formula':'NH4'})
        compounds.append({'id':'hpo4_c0','formula':'HPO4'})
        compounds.append({'id':'hs_c0','formula':'HS'})
        compounds.append({'id':'h_c0','formula':'H'})
        compounds.append({'id':'e_c0','formula':'e-'})
        compounds.append({'id':'acceptor_c0','formula':'O2'})
        compounds.append({'id':'biom_c0','formula':'CH1.8N0.2O0.5'})
        comp_df = pd.DataFrame(compounds, columns=comp_cols)
        comp_df.to_csv(fout, sep='\t', index=False)

    def create_rxn_file_fba_model(self, stoich_mat, fout):
        rxn_cols = ['id','direction','compartment','gpr','name','enzyme','deltag','reference','equation',
            'definition','ms id','bigg id','kegg id','kegg pathways','metacyc pathways']
        
        stoich_colnames = ["donor","h2o","hco3","nh4","hpo4","hs","h","e","acceptor","biom"]

        def generate_equation(r):
            reactants = []
            products = []
            for col in stoich_colnames:
                if col=='donor': name = self.mf2id[r.name]
                else: name = col

                if r[col] == 0: continue
                elif r[col] > 0:
                    products.append('({0})  {1}[c0]'.format(r[col], name))
                else:
                    reactants.append('({0})  {1}[c0]'.format(-r[col], name))
            return '{} <=> {}'.format(' + '.join(reactants), ' + '.join(products))
        
        reactions = []
        for i, eq in enumerate(stoich_mat.apply(generate_equation, axis=1).tolist()):
            reactions.append({'id':'xrxn{}_c0'.format(i+1),'equation':eq})
        
        rxn_df = pd.DataFrame(reactions,columns=rxn_cols)
        rxn_df.to_csv(fout, sep='\t', index=False)

    def create_media_file(self, media_file):
        media_cols = ['compounds','name','formula','minFlux','maxFlux','concentration']
        media_compounds = [{'compounds':_id,'formula':self.id2mf[_id],'name':self.id2mf[_id],'minFlux':-1000,'maxFlux':1000,'concentration':1} for _id in self.id2mf]
        media_df = pd.DataFrame(media_compounds, columns=media_cols)
        media_df.to_csv(media_file, sep='\t', index=False)

    def plot_lambda_dist(self, fout='lambda_dist.png'):
        if self.thermo is not None:
            plt.close('all')
            g = sns.distplot(self.thermo.lambda_O2, label=r'$\lambda$')
            # g = sns.distplot(self.thermo.lambda_HCO3, label='HCO3')
            plt.xlabel(r'$\lambda$', fontsize=15)
            plt.ylabel('Distribution', fontsize=15)
            plt.xlim([0,0.3])
            plt.legend(fontsize=15)
            plt.tight_layout()
            if fout: plt.savefig(fout)
        else:
            print('[Warning] "plot_lambda_dist" requires self.thermo. Please use run().')

    def plot_delta_gibb_dist(self, colname, label, fout='dist.png'):
        if self.thermo is not None:
            plt.close('all')
            g = sns.distplot(self.thermo[colname], label=label)
            g.set_xlabel(label+"[kJ/C-mol]", fontsize=15)
            g.set_ylabel('Distribution', fontsize=15)
            plt.legend(fontsize=15)
            plt.tight_layout()
            if fout: plt.savefig(fout)
        else:
            print('[Warning] "plot_lambda_dist" requires self.thermo. Please use run().')

    def get_summary(self, colname):
        if self.thermo is not None:
            return (self.thermo[colname].mean(),
                    self.thermo[colname].std(ddof=1),
                    self.thermo[colname].median())
        else:
            print('[Warning] "plot_lambda_dist" requires self.thermo. Please use run().')
            return (np.nan, np.nan, np.nan)

    def plot_van_krevelen(self, fout):
        df = self._assigned_tbl.copy()
        plt.figure(figsize=(10,8))
        df["H:C"] = df.H / df.C
        df["O:C"] = df.O / df.C

        g = sns.scatterplot("O:C", "H:C", hue="Class", alpha=1, s=15, data=df)
        g.set_xlabel("O:C", fontsize=15)
        g.set_ylabel("H:C", fontsize=15)
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=15)
        plt.tight_layout()
        
        plt.savefig(fout)
        
    def average_by_lambda_bins(self, n_bins=10, cutoff=5):
        '''
        average compositions per each bin in the lambda distribution after 
        filtering out the two-side tails by a cutoff percent (%)
        '''
        assert 0 <= cutoff < 100, "cutoff must be 0 <= cutoff < 100"
        assert 0 < n_bins, "n_bins must be 0 < n_bins"

        # data
        lambda_dist = self.thermo.lambda_O2.values
        comp_df = self._assigned_tbl[REQUIRED_COLUMNS].copy()
        # th_lambda[0] --> lambda_O2
        comp_df['lambda'] = comp_df.apply(lambda x: self.all_stoich[x.name].th_lambda[0], axis=1)

        # get the boundary
        if cutoff > 0:
            lambda_min = np.percentile(lambda_dist, cutoff)
            lambda_max = np.percentile(lambda_dist, 100-cutoff)
        elif cutoff == 0:
            lambda_min = 0
            lambda_max = np.amax(lambda_dist)
        print('lambda_min', lambda_min, 'lambda_max', lambda_max)

        # binning
        bins = np.linspace(lambda_min, lambda_max, n_bins+1)
        labels = ['Bin{}'.format(i+1) for i in range(n_bins)]
        comp_df['Class'] = pd.cut(comp_df['lambda'], bins=bins, labels=labels)
        tdf = comp_df[comp_df['Class'].notnull()]

        new_comp = tdf.groupby('Class').mean()
        new_comp['Na'] = 0
        new_comp['C13'] = 0

        return new_comp.reset_index()



class ThermoStoichiometry(object):
    """extract thermo stoichiometry from a chemical formula"""
    
    # def __init__(self, chemical_formula, ignore_isotopes={'C13':1}):
    #     '''
    #         chemical_formula: chemical formula (e.g. C27H15O6NC13SP)
    #         ignore_isotopes: ignore if you have this in the '>index' of the 
    #             chemical_formula. it should follow the name convention
    #     '''
    #     super(ThermoStoichiometry, self).__init__()
    #     self.chemical_formula = chemical_formula.upper()
    #     if ignore_isotopes is not None:
    #         self.ignore_isotopes = ignore_isotopes
    #         self.__ignore_isotopes()
    #     self.chemical_composition = self.extract_composition()

    def __init__(self, chemical_composition):
        """extract thermo stoichiometry from a given chemical composition"""
        super(ThermoStoichiometry, self).__init__()
        self.chemical_composition = chemical_composition

    # def __ignore_isotopes(self):
    #     '''remove the isotope tags in the end of the chemical formula
    #     '''
    #     for isotope_tag in self.ignore_isotopes:
    #         try:
    #             if self.chemical_formula.index(isotope_tag)>=self.ignore_isotopes[isotope_tag]:
    #                 self.chemical_formula = self.chemical_formula.replace(isotope_tag, '')
    #         except ValueError as e: # if not isotope_tag in the formular
    #             continue

    # def extract_composition(self):
    #     '''extract the chemimcal composition from the chemical formula
    #     '''
    #     # search number of elements
    #     chem_comp = {}
    #     for e in CHEMICAL_ELEMENTS:
    #         composition = re.search('{}(\d*)'.format(e), self.chemical_formula)

    #         if composition:
    #             n_element = composition.group(1)
    #             if n_element=='': chem_comp[e] = 1
    #             else: chem_comp[e] = int(n_element)
    #         else: chem_comp[e] = 0
    #     return chem_comp

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
        stoich_electron_acceptor[6] = -4  # h+
        stoich_electron_acceptor[7] = -4  # e-
        stoich_electron_acceptor[1] =  2  # h2o
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
                stoichAnStar_O2 = stoichAnStarB+(1/float(a))*stoich_electron_donor
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
        nosc = -ne/float(a)+4  # nominal oxidataion state of carbon 
        delGcox0PerC = 60.3-28.5*nosc  # kJ/C-mol
        delGcox0 = delGcox0PerC*a*np.abs(stoich_electron_donor[0])  # kJ/rxn

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
            [delGcox0PerC,delGcox0,delGcox,delGcat0,delGcat,delGan0_O2,delGan0_HCO3,\
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

