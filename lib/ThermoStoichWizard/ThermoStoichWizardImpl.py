# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
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
        #  run
        #######################################################################
        print ("Input parameter", params['input_tbl'])
        dfu = DataFileUtil(self.callback_url)
        input_tbl = dfu.get_objects({'object_refs': [params['input_tbl']]})['data'][0]

        print(input_tbl['data']['attributes'])
        for peak in input_tbl['data']['instances']:
            print(peak, input_tbl['data']['instances'][peak])
            break
        print(input_tbl['info'])

        from ThermoStoichWizard.ThermoStoichiometry import ThermoStoichiometry

        mf = 'C35H32O6C13S2'
        therm = ThermoStoichiometry(mf)
        print(therm.extract_composition())
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
