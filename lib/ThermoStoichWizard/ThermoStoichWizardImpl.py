# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
import pandas as pd

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.fba_toolsClient import fba_tools

from ThermoStoichWizard.ThermoStoichiometry import ThermoStoichiometry

#END_HEADER


class ThermoStoichWizard:
    '''
    Module Name:
    ThermoStoichWizard

    Module Description:
    A KBase module: ThermoStoichWizard
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/coldfire79/ThermoStoichWizard.git"
    GIT_COMMIT_HASH = "4d5870ae2d8b1d6baef2cb17e8e3f0869b596cd2"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = os.path.abspath(config['scratch'])
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_ThermoStoichWizard(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_ThermoStoichWizard
        uuid_string = str(uuid.uuid4())

        #######################################################################
        #  check out the input table
        #######################################################################
        print ("Input parameter", params['input_tbl'])
        dfu = DataFileUtil(self.callback_url)
        input_tbl = dfu.get_objects({'object_refs': [params['input_tbl']]})['data'][0]

        print(input_tbl['data']['attributes'])
        for peak in input_tbl['data']['instances']:
            print(peak, input_tbl['data']['instances'][peak])
            break
        print(input_tbl['info'])

        #######################################################################
        #  compute thermo stoichiometry
        #######################################################################

        mf = 'C35H32O6C13S2'
        therm = ThermoStoichiometry(mf)
        print(therm.extract_composition())

        #######################################################################
        #  create the tsv files for fba
        #######################################################################
        comp_cols = ['id','name','formula','charge','inchikey','smiles','deltag','kegg id','ms id']
        rxn_cols = ['id','direction','compartment','gpr','name','enzyme','deltag','reference','equation','definition','ms id','bigg id','kegg id','kegg pathways','metacyc pathways']
        
        compounds = []
        compounds.append({'id':'comp1_c0','formula':'C35H32O6C13S2'})
        compounds.append({'id':'h2o_c0','formula':'H2O'})
        compounds.append({'id':'hco3_c0','formula':'HCO3'})
        compounds.append({'id':'nh4_c0','formula':'NH4'})
        compounds.append({'id':'hpo4_c0','formula':'HPO4'})
        compounds.append({'id':'hs_c0','formula':'HS'})
        compounds.append({'id':'h_c0','formula':'H'})


        comp_df = pd.DataFrame(compounds, columns=comp_cols)

        reactions = []
        reactions.append({'id':'biomass_c0','equation':'(1)  comp1[c0] <=> (1)  h2o[c0] + (1)  hco3[c0] + (1)  nh4[c0] + (1)  hpo4[c0] + (1)  hs[c0] + (1)  h[c0]'})

        rxn_df = pd.DataFrame(reactions,columns=rxn_cols)

        compounds_file = os.path.join(self.shared_folder, "temp_comps.tsv")
        reactions_file = os.path.join(self.shared_folder, "temp_rxns.tsv")

        comp_df.to_csv(compounds_file, sep='\t', index=False)
        rxn_df.to_csv(reactions_file, sep='\t', index=False)

        #######################################################################
        #  generate fbamodel
        #######################################################################
        fba_param = {
            # 'model_name':'model' + uuid_string,
            'file_type':'tsv',
            'compounds_file':{'path': compounds_file},
            'model_file':{'path': reactions_file},
            'biomass':['biomass_c0'],
            'model_name': "ThermoStoic_model",
            'workspace_name': params['workspace_name']
        }
        fbaobj = fba_tools(self.callback_url)
        fba_model_wref = fbaobj.tsv_file_to_model(p=fba_param)
        print('fba_model:', fba_model_wref)

        #######################################################################
        #  create the tsv files for media
        #######################################################################
        media_cols = ['compounds','name','formula','minFlux','maxFlux','concentration']

        media_compounds = []
        media_compounds.append({'id':'comp1_c0','formula':'C35H32O6C13S2','name':'comp1','minFlux':-1000,'maxFlux':1000,'concentration':0})
        media_compounds.append({'id':'comp2_c0','formula':'C35H32O6S2','name':'comp2','minFlux':-1000,'maxFlux':1000,'concentration':0})
        media_compounds.append({'id':'comp3_c0','formula':'C35H13S2','name':'comp3','minFlux':-1000,'maxFlux':1000,'concentration':0})

        media_df = pd.DataFrame(media_compounds, columns=media_cols)

        media_tsv_file = os.path.join(self.shared_folder, "temp_media.tsv")

        media_df.to_csv(media_tsv_file, sep='\t', index=False)
        

        #######################################################################
        #  generate media
        #######################################################################
        media_param = {
            'file_type':'tsv',
            'media_file':{'path': media_tsv_file},
            'media_name': "ThermoStoic_media",
            'workspace_name': params['workspace_name']
        }
        media_wref = fbaobj.tsv_file_to_media(p=media_param)
        print('media:',media_wref)

        #######################################################################
        # 
        #######################################################################
        html_folder = os.path.join(self.shared_folder, 'html')
        os.mkdir(html_folder)

        html_str = "<html><head>Thermo Stoich Wizard Report</head><body><br><br>{}:{}</body></html>"
        html_str = html_str.format(mf, therm.extract_composition())

        with open(os.path.join(html_folder, "index.html"), 'w') as index_file:
            index_file.write(html_str)

        report = KBaseReport(self.callback_url)
        html_dir = {
            'path': html_folder,
            'name': 'index.html',  # MUST match the filename of your main html page
            'description': 'Thermo Stoich Wizard Report'
        }
        report_info = report.create_extended_report({
            'html_links': [html_dir],
            'direct_html_link_index': 0,
            'report_object_name': 'thermo_stoich_wizard_report_' + uuid_string,
            'workspace_name': params['workspace_name']
        })
        # report_info = report.create({'report': {'objects_created':[],
        #                                         'text_message': "OK"},
        #                                         'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_ThermoStoichWizard

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_ThermoStoichWizard return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
