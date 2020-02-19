# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
import pandas as pd

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.fba_toolsClient import fba_tools

from ThermoStoichWizard.ThermoStoichiometry import ThermoStoichiometry, FTICRResult

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
        objects_created = []
        output_files = []
        
        #######################################################################
        #  check out the input table
        #######################################################################
        print ("Input parameter", params['input_tbl'])
        dfu = DataFileUtil(self.callback_url)
        input_tbl = dfu.get_objects({'object_refs': [params['input_tbl']]})['data'][0]

        # # investigate input_tbl
        # for peak in input_tbl['data']['instances']:
        #     print(peak, input_tbl['data']['instances'][peak])
        #     break
        # print(input_tbl['info'])

        tbl_cols = [info['attribute'] for info in input_tbl['data']['attributes']]
        tbl_df = pd.DataFrame.from_dict(input_tbl['data']['instances'], orient='index', columns=tbl_cols)

        #######################################################################
        #  compute thermo stoichiometry
        #######################################################################
        fticr = FTICRResult(tbl_df)
        fticr.to_csv(os.path.join(self.shared_folder, "input_compounds.csv"))

        fticr.run()
        output_folder = os.path.join(self.shared_folder, 'csv')
        os.mkdir(output_folder)
        fticr.save_result_files(output_folder)
        output_filenames = ["stoichD","stoichA","stoichCat","stoichAn_O2","stoichAn_HCO3","stoichMet_O2","stoichMet_HCO3","thermodynamic_props"]
        output_files = [{'path': output_folder+'/{}.csv'.format(n), 'name': n, 'label': n, 'description': n} for n in output_filenames]

        # filter out the unassigned peaks
        num_peaks = fticr.num_peaks
        num_cpds = fticr.num_cpds

        print('num_peaks:{}, num_cpds:{}'.format(num_peaks, num_cpds))

        # #######################################################################
        # #  create the tsv files for fba
        # #######################################################################
        fticr.create_fba_model_files(self.shared_folder)

        #######################################################################
        #  generate fbamodel
        #######################################################################
        fbaobj = fba_tools(self.callback_url)
        def generate_fbamodel(fbaobj, fticr, model_prefix, workspace_name, compounds_file, reactions_file):
            fba_param = {
                # 'model_name':'model' + uuid_string,
                'file_type':'tsv',
                'compounds_file':{'path': compounds_file},
                'model_file':{'path': reactions_file},
                'biomass':['rxn1_c0'],  # TODO: how to define a biomass reaction
                'model_name': "Thermo_{}_{}".format(model_prefix, uuid_string),
                'workspace_name': workspace_name # 
            }
            fba_model_wref = fbaobj.tsv_file_to_model(p=fba_param)
            print('fba_model:', fba_model_wref)
            return fba_model_wref

        stoichiometries = ["stoichD","stoichA","stoichCat","stoichAn_O2","stoichAn_HCO3","stoichMet_O2","stoichMet_HCO3"]
        for stoich in stoichiometries:
            fba_model_wref = generate_fbamodel(fbaobj, fticr, model_prefix=stoich,
                workspace_name=params['workspace_name'],
                compounds_file=os.path.join(self.shared_folder, "temp_comps.tsv"),
                reactions_file=os.path.join(self.shared_folder, "temp_{}.tsv".format(stoich)))
            objects_created.append({'ref': fba_model_wref,
                'description': "FBA model for ".format(stoich)})
        #######################################################################
        #  create the tsv files for media
        #######################################################################
        media_tsv_file = os.path.join(self.shared_folder, "temp_media.tsv")
        fticr.create_media_file(media_tsv_file)
        #######################################################################
        #  generate media
        #######################################################################
        media_param = {
            'file_type':'tsv',
            'media_file':{'path': media_tsv_file},
            'media_name': "ThermoStoic_media_{}".format(uuid_string),
            'workspace_name': params['workspace_name']
        }
        media_wref = fbaobj.tsv_file_to_media(p=media_param)
        print('media:', media_wref)
        objects_created.append({'ref': media_wref,
            'description': "Media object contains the initial condition."})

        #######################################################################
        # figures
        #######################################################################
        fig_folder = os.path.join(self.shared_folder, 'fig')
        os.mkdir(fig_folder)
        lambda_dist_path = os.path.join(fig_folder, "lambda_dist.png")
        fticr.plot_lambda_dist(fout=lambda_dist_path)
        delGcox_dist_path = os.path.join(fig_folder, "delGcox_dist.png")
        fticr.plot_delta_gibb_dist('delGcox', r'$\Delta G^{Cox}$', delGcox_dist_path)
        delGcat_dist_path = os.path.join(fig_folder, "delGcat_dist.png")
        fticr.plot_delta_gibb_dist('delGcat', r'$\Delta G^{Cat}$', delGcat_dist_path)

        #######################################################################
        # html report
        #######################################################################
        html_folder = os.path.join(self.shared_folder, 'html')
        os.mkdir(html_folder)

        # html_str = '\
        # <html>\
        #   <head><title>Thermo Stoich Wizard Report</title>\
        #   <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/jquery.dataTables.min.css">\
        #   </head>\
        #   <body><br><br>{}:{}<br>{}:{}<div id="cpd_tbl">{}</div>\
        #   </body>\
        #   <script src="https://code.jquery.com/jquery-3.3.1.js"></script> \
        #   <script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js"></script>\
        #   <script type="text/javascript">\
        #     $(document).ready( function () {{\
        #       $("#cpd_tbl table").DataTable();\
        #     }});\
        #   </script>\
        # </html>'

        # html_str = '\
        # <html>\
        #   <head><title>Thermo Stoich Wizard Report</title>\
        #   <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.12/css/jquery.dataTables.min.css">\
        #   </head>\
        #   <body><br><br>{}:{}<br>{}:{}\
        #   <div><img src="{}" alt="Lambda distribution"></div><div><img src="{}" alt="delGcat distribution"></div><div><img src="{}" alt="delGcox distribution"></div>\
        #   </body>\
        #   <script src="https://code.jquery.com/jquery-3.3.1.js"></script> \
        #   <script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js"></script>\
        #   </script>\
        # </html>'
        html_str = '<!doctype html>'
        html_str += '<html lang="en">'
        html_str += '<head>'
        html_str += '<meta charset="utf-8">'
        html_str += '<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">'
        html_str += '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">'
        html_str += '<title>Thermo Stoich Wizard Report</title>'
        html_str += '</head>'
        html_str += '<body>'
        html_str += '<h1>Thermo Stoich Wizard Report</h1>'
        html_str += '<h4>{}:{}</h4>'
        html_str += '<h4>{}:{}</h4>'
        html_str += '<div class="container">'
        html_str += '<div class="row">'
        html_str += '<div class="col-md-4">'
        html_str += '<div class="card mb-4 box-shadow">'
        html_str += '<img class="card-img-top" alt="lambda_dist" src="../fig/lambda_dist.png" style="height: 225px; width: 100%; display: block;">'
        html_str += '<div class="card-body">'
        html_str += '<p class="card-text">Energy coupling thermodynamic parameter</p>'
        html_str += '<div class="d-flex justify-content-between align-items-center">'
        html_str += '<div class="btn-group">'
        html_str += '<button type="button" class="btn btn-sm btn-outline-secondary">Save</button>'
        html_str += '</div>'
        html_str += '<small class="text-muted">9 mins</small>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '<div class="col-md-4">'
        html_str += '<div class="card mb-4 box-shadow">'
        html_str += '<img class="card-img-top" alt="delGcat_dist" src="../fig/delGcat_dist.png" style="height: 225px; width: 100%; display: block;">'
        html_str += '<div class="card-body">'
        html_str += '<p class="card-text">Gibbs free energy change for catabolic reaction</p>'
        html_str += '<div class="d-flex justify-content-between align-items-center">'
        html_str += '<div class="btn-group">'
        html_str += '<button type="button" class="btn btn-sm btn-outline-secondary">Save</button>'
        html_str += '</div>'
        html_str += '<small class="text-muted">9 mins</small>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '<div class="col-md-4">'
        html_str += '<div class="card mb-4 box-shadow">'
        html_str += '<img class="card-img-top" alt="delGcox_dist" src="../fig/delGcox_dist.png" style="height: 225px; width: 100%; display: block;">'
        html_str += '<div class="card-body">'
        html_str += '<p class="card-text">Gibbs free energy change for an electron donor half reaction</p>'
        html_str += '<div class="d-flex justify-content-between align-items-center">'
        html_str += '<div class="btn-group">'
        html_str += '<button type="button" class="btn btn-sm btn-outline-secondary">Save</button>'
        html_str += '</div>'
        html_str += '<small class="text-muted">9 mins</small>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '</div>'
        html_str += '<script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>'
        html_str += '<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>'
        html_str += '<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>'
        html_str += '</body>'
        html_str += '</html>'
        
        html_str = html_str.format("Number of peaks", num_peaks,
                                   "Number of compounds", num_cpds)

        with open(os.path.join(html_folder, "index.html"), 'w') as index_file:
            index_file.write(html_str)

        report = KBaseReport(self.callback_url)
        html_dir = {
            'path': html_folder,
            'name': 'index.html',  # MUST match the filename of your main html page
            'description': 'Thermo Stoich Wizard Report'
        }
        report_info = report.create_extended_report({
            'objects_created': objects_created,
            'file_links': output_files,
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
