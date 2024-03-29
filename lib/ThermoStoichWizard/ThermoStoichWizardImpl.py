# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.fba_toolsClient import fba_tools

from ThermoStoichWizard.ThermoStoichiometry import ThermoStoichiometry, FTICRResult
from ThermoStoichWizard.LambdaAnalysis import LambdaAnalysis

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
    VERSION = "0.0.2"
    GIT_URL = "https://github.com/coldfire79/ThermoStoichWizard.git"
    GIT_COMMIT_HASH = "1c3776e6e16251db8a131ca5a9397b0c746588c2"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = os.path.abspath(config['scratch'])
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        self.lambda_analysis = LambdaAnalysis(self.config)
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

        n_lambda_bins = int(params['n_lambda_bins'])
        lambda_cutoff = float(params['lambda_cutoff'])

        
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
        # output_filenames = ["stoichD","stoichA","stoichCat","stoichAn_O2","stoichAn_HCO3","stoichMet_O2","stoichMet_HCO3","thermodynamic_props"]
        output_filenames = ["stoichMet_O2","thermodynamic_props"]
        output_files = [{
                'path': output_folder+'/{}.csv'.format(n),
                'name': '{}.csv'.format(n),
                'label': n, 'description': n
            } for n in output_filenames]

        # filter out the unassigned peaks
        num_peaks = fticr.num_peaks
        num_cpds = fticr.num_cpds

        print('num_peaks:{}, num_cpds:{}'.format(num_peaks, num_cpds))

        #######################################################################
        #  average compositions by lambda bins
        #######################################################################
        if params['bin_method'] == "cumulative":
            new_comp = fticr.average_by_lambda_bins(n_bins=n_lambda_bins, cutoff=lambda_cutoff)
        elif params['bin_method'] == "uniform":
            new_comp = fticr.average_by_lambda_bins_uniform(n_bins=n_lambda_bins, cutoff=lambda_cutoff)
        else:
            raise("bin_method was wrong:", params['bin_method'])
        average_comp_path = os.path.join(output_folder, "avg_comp_from_lambda_bins.csv")
        new_comp.to_csv(average_comp_path)

        output_files.append({'path': average_comp_path,
                            'name': 'avg_comp_from_lambda_bins.csv',
                            'label': 'average compositions for each lambda bin',
                            'description': 'average compositions for each lambda bin'})

        # compute the reactions for bin averaged compositions
        new_fticr = FTICRResult(new_comp, dtype=np.float)

        new_fticr.run()
        selected_folder = os.path.join(self.shared_folder, 'bin_avg')
        os.mkdir(selected_folder)
        new_fticr.save_result_files(selected_folder)
        output_filenames = ["stoichMet_O2"]
        output_files += [{
            'path': selected_folder+'/{}.csv'.format(n),
            'name': '{}_from_lambda_bins.csv'.format(n),
            'label': '{}_from_lambda_bins'.format(n),
            'description': '{}_from_lambda_bins'.format(n),
        } for n in output_filenames]

        #######################################################################
        #  create the tsv files for fba
        #######################################################################
        fticr.create_fba_model_files(self.shared_folder)
        new_fticr.create_fba_model_files(self.shared_folder, prefix='bin_avg')

        #######################################################################
        #  generate fbamodel
        #######################################################################
        fbaobj = fba_tools(self.callback_url)
        def generate_fbamodel(fbaobj, model_prefix, workspace_name, compounds_file, reactions_file):
            fba_param = {
                # 'model_name':'model' + params['output_surfix'],
                'file_type':'tsv',
                'compounds_file':{'path': compounds_file},
                'model_file':{'path': reactions_file},
                'biomass':['xrxn1_c0'],  # TODO: how to define a biomass reaction
                'model_name': "{}_{}".format(model_prefix, params['output_surfix']),
                'workspace_name': workspace_name # 
            }
            fba_model_wref = fbaobj.tsv_file_to_model(p=fba_param)
            print('fba_model:', fba_model_wref)
            return fba_model_wref

        # stoichiometries = ["stoichD","stoichA","stoichCat","stoichAn_O2","stoichAn_HCO3","stoichMet_O2","stoichMet_HCO3"]
        stoichiometries = ["stoichMet_O2"]
        for stoich in stoichiometries:
            fba_model_wref = generate_fbamodel(fbaobj, model_prefix=stoich,
                workspace_name=params['workspace_name'],
                compounds_file=os.path.join(self.shared_folder, "temp_comps.tsv"),
                reactions_file=os.path.join(self.shared_folder, "temp_{}.tsv".format(stoich)))
            objects_created.append({'ref': fba_model_wref['ref'],
                'description': "FBA model for {}".format(stoich)})

            fba_model_wref = generate_fbamodel(fbaobj, model_prefix="Bin_Averaged_"+stoich,
                workspace_name=params['workspace_name'],
                compounds_file=os.path.join(self.shared_folder, "bin_avg_comps.tsv"),
                reactions_file=os.path.join(self.shared_folder, "bin_avg_{}.tsv".format(stoich)))
            objects_created.append({'ref': fba_model_wref['ref'],
                'description': "FBA model for {}".format(stoich)})
        #######################################################################
        #  create the tsv files for media
        #######################################################################
        # media_tsv_file = os.path.join(self.shared_folder, "temp_media.tsv")
        # fticr.create_media_file(media_tsv_file)

        #######################################################################
        #  generate media
        #######################################################################
        # media_param = {
        #     'file_type':'tsv',
        #     'media_file':{'path': media_tsv_file},
        #     'media_name': "ThermoStoic_media_{}".format(params['output_surfix']),
        #     'workspace_name': params['workspace_name']
        # }
        # media_wref = fbaobj.tsv_file_to_media(p=media_param)
        # print('media:', media_wref)
        # objects_created.append({'ref': media_wref['ref'],
        #     'description': "Media object contains the initial condition."})

        #######################################################################
        # html report
        #######################################################################
        html_folder = os.path.join(self.shared_folder, 'html')
        os.mkdir(html_folder)

        #######################################################################
        # figures
        #######################################################################
        if "Class" in fticr._assigned_tbl.columns:
            van_krevelen_available = True
        else:
            van_krevelen_available = False

        if van_krevelen_available:
            van_krevelen_path = os.path.join(html_folder, "van_krevelen.png")
            # fticr.plot_van_krevelen(fout=van_krevelen_path)
            
            # fig, ax = plt.subplots(1,2,figsize=(12,6),sharex=True,sharey=True)
            plt.figure(figsize=(7,5))
            df = fticr._assigned_tbl.copy()
            
            df["H:C"] = df.H / df.C
            df["O:C"] = df.O / df.C

            g1 = sns.scatterplot("O:C", "H:C", hue="Class", alpha=1, s=15, data=df)
            g1.set_xlabel("O:C", fontsize=15)
            g1.set_ylabel("H:C", fontsize=15)
            plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=10)
            plt.tight_layout()
            plt.savefig(van_krevelen_path)

            van_krevelen_lambda_bins_path = os.path.join(html_folder, "van_krevelen_by_lambda_bins.png")
            plt.figure(figsize=(7,5))
            new_comp["H:C"] = new_comp.H / new_comp.C
            new_comp["O:C"] = new_comp.O / new_comp.C
            g = sns.scatterplot("O:C", "H:C", hue="Class", s=100, data=new_comp)
            g.set_xlabel("O:C", fontsize=15)
            g.set_ylabel("H:C", fontsize=15)
            g.set_xlim(g1.get_xlim())
            g.set_ylim(g1.get_ylim())
            plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=12)
            plt.tight_layout()
            # plt.savefig(van_krevelen_path)
            plt.savefig(van_krevelen_lambda_bins_path)

        lambda_dist_path = os.path.join(html_folder, "lambda_dist.png")
        fticr.plot_lambda_dist(fout=lambda_dist_path)
        # delGcat0_dist_path = os.path.join(html_folder, "delGcat0_dist.png")
        # fticr.plot_delta_gibb_dist('delGcat0', r'$\Delta G_{Cox}^0$', delGcat0_dist_path)
        # delGcat_dist_path = os.path.join(html_folder, "delGcat_dist.png")
        # fticr.plot_delta_gibb_dist('delGcat', r'$\Delta G_{Cox}$', delGcat_dist_path)
        delGcox0_dist_path = os.path.join(html_folder, "delGcox0_dist.png")
        fticr.plot_delta_gibb_dist('delGcox0PerC', r'$\Delta G_{Cox}^0$', delGcox0_dist_path)
        # delGcox_dist_path = os.path.join(html_folder, "delGcox_dist.png")
        # fticr.plot_delta_gibb_dist('delGcox', r'$\Delta G_{Cox}$', delGcox_dist_path)

        if van_krevelen_available:
            output_files.append({'path': van_krevelen_path, 'name': 'van_krevelen.png',
                'label': 'van Krevelen diagram for compounds', 'description': 'van Krevelen diagram for compounds'})
            output_files.append({'path': van_krevelen_lambda_bins_path, 'name': 'van_krevelen_by_lambda_bins.png',
                'label': 'van Krevelen diagram for each lambda bin', 'description': 'van Krevelen diagram for each lambda bin'})

        output_files.append({'path': lambda_dist_path, 'name': 'lambda_dist.png',
            'label': 'lambda distribution', 'description': 'lambda distribution'})
        # output_files.append({'path': delGcat0_dist_path, 'name': 'delGcat0_dist.png',
        #     'label': 'delGcat0 distribution',
        #     'description': 'Gibbs free energy change for an electron donor half reaction'})
        # output_files.append({'path': delGcat_dist_path, 'name': 'delGcat_dist.png',
        #     'label': 'delGcat distribution', 'description': 'Gibbs free energy change for catabolic reaction'})
        output_files.append({'path': delGcox0_dist_path, 'name': 'delGcox0_dist.png',
            'label': 'delGcox0 distribution',
            'description': 'Gibbs energies for the oxidation half reactions'})
        # output_files.append({'path': delGcox_dist_path, 'name': 'delGcox_dist.png',
        #     'label': 'delGcox distribution',
        #     'description': 'Gibbs energies for the oxidation half reactions'})

        summary_str = '<ul class="list-group list-group-flush">'
        summary_str += '<li class="list-group-item">Average: {:.3f}</li>'
        summary_str += '<li class="list-group-item">Standard deviation: {:.3f}</li>'
        summary_str += '<li class="list-group-item">Median: {:.3f}</li>'
        summary_str += '</ul>'

        html_str = '<div class="col-md-6">'
        html_str += '<div class="card mb-6 box-shadow">'
        html_str += '<img class="card-img-top" alt="lambda_dist" src="lambda_dist.png" style="width: 100%; display: block;">'
        html_str += '<div class="card-body">'
        html_str += '<p class="card-text">Energy coupling thermodynamic parameter</p>'
        html_str += '</div>'
        html_str += summary_str.format(*fticr.get_summary('lambda_O2'))
        html_str += '</div>'
        html_str += '</div>'

        html_str += '<div class="col-md-6">'
        html_str += '<div class="card mb-6 box-shadow">'
        html_str += '<img class="card-img-top" alt="delGcox0_dist" src="delGcox0_dist.png" style="width: 100%; display: block;">'
        html_str += '<div class="card-body">'
        html_str += '<p class="card-text">Gibbs free energy change for catabolic reaction</p>'
        html_str += '</div>'
        html_str += summary_str.format(*fticr.get_summary('delGcox0PerC'))
        html_str += '</div>'
        html_str += '</div>'
        
        # html_str += '<div class="col-md-4">'
        # html_str += '<div class="card mb-4 box-shadow">'
        # html_str += '<img class="card-img-top" alt="delGcat_dist" src="delGcat_dist.png" style="width: 100%; display: block;">'
        # html_str += '<div class="card-body">'
        # html_str += '<p class="card-text">Gibbs free energy change for an electron donor half reaction</p>'
        # html_str += '</div>'
        # html_str += summary_str.format(*fticr.get_summary('delGcat'))
        # html_str += '</div>'
        # html_str += '</div>'
        
        if van_krevelen_available:
            html_str += '<div class="col-md-6">'
            html_str += '<div class="card mb-6 box-shadow">'
            html_str += '<img class="card-img-top" alt="van_krevelen" src="van_krevelen.png" style="width: 100%; display: block;">'
            html_str += '<div class="card-body">'
            html_str += '<p class="card-text">Van Krevelen diagram for all compositions</p>'
            html_str += '</div>'
            html_str += '</div>'
            html_str += '</div>'

            html_str += '<div class="col-md-6">'
            html_str += '<div class="card mb-6 box-shadow">'
            html_str += '<img class="card-img-top" alt="van_krevelen_by_lambda_bins" src="van_krevelen_by_lambda_bins.png" style="width: 100%; display: block;">'
            html_str += '<div class="card-body">'
            html_str += '<p class="card-text">Van Krevelen Diagram for average compositions of lambda bins</p>'
            html_str += '</div>'
            html_str += '</div>'
            html_str += '</div>'

        with open(os.path.join(os.path.dirname(__file__), 'templates', 'template.html'),
                  'r') as template_file:
            report_html = template_file.read()
            report_html = report_html.replace('Number of peaks:', 'Number of peaks: {}'.format(num_peaks))
            report_html = report_html.replace('Number of compounds:', 'Number of compounds: {}'.format(num_cpds))
            report_html = report_html.replace('<!--[Results]-->', html_str)
            
        with open(os.path.join(html_folder, "index.html"), 'w') as index_file:
            index_file.write(report_html)

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
            'report_object_name': 'thermo_stoich_wizard_report_' + params['output_surfix'],
            'workspace_name': params['workspace_name']
        })
        
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

    def run_lambda_analysis(self, ctx, params):
        """
        run_lambda_analysis: perform lambda analysis
        :param params: instance of type "LambdaParams" -> structure:
           parameter "lambda_tbl" of type "obj_ref" (An X/Y/Z style
           reference), parameter "stoich_tbl" of type "obj_ref" (An X/Y/Z
           style reference), parameter "vh_cs" of String, parameter "vh_o2"
           of String, parameter "workspace_name" of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_lambda_analysis
        output = self.lambda_analysis.run(params)
        #END run_lambda_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_lambda_analysis return value ' +
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
