import logging
import os
import uuid
import pandas as pd
import numpy as np
import itertools
from scipy.stats import pearsonr

import matplotlib.pyplot as plt
import seaborn as sns

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

class LambdaAnalysis(object):
    """docstring for LambdaAnalysis"""
    def __init__(self, config):
        super(LambdaAnalysis, self).__init__()
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

    def run(self, params):
        print("run lambda analysis")
        vh_cs = [float(i) for i in params["vh_cs"].split(",")]
        vh_o2 = [float(i) for i in params["vh_o2"].split(",")]
        
        df = self._fetch_df_from_refs([params["lambda_tbl"], params["stoich_tbl"]])

        mu_max = 1

        for vhcs, vho2 in itertools.product(vh_cs, vh_o2):
            df['r_biom'] = mu_max * np.exp(-np.abs(df['donor'])/vhcs)* np.exp(-np.abs(df['acceptor'])/vho2)
            df['r_o2'] = np.abs(df['acceptor']) * df['r_biom']
            df['r_hco3'] = np.abs(df['hco3']) * df['r_biom']

            r_lambda_rbiom = pearsonr(df['lambda_O2'], df['r_biom'])
            r_lambda_ro2 = pearsonr(df['lambda_O2'], df['r_o2'])
            r_lambda_rhco3 = pearsonr(df['lambda_O2'], df['r_hco3'])

            print(vhcs, vho2)
            print(r_lambda_rbiom, r_lambda_ro2, r_lambda_rhco3)



        # self.plot(df)
        
        report = KBaseReport(self.callback_url)
        report_info = report.create_extended_report({
            # # 'objects_created': objects_created,
            # 'file_links': output_files,
            # 'html_links': [html_dir],
            # 'direct_html_link_index': 0,
            # 'report_object_name': 'miia_report_' + params['output_suffix'],
            'workspace_name': params['workspace_name']
        })

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }

        return output

    def _fetch_df_from_refs(self, object_refs):
        dfu = DataFileUtil(self.callback_url)
        tables = dfu.get_objects({'object_refs': object_refs})['data']
        lambda_df = self._fetch_df_from_json(tables[0])
        stoich_df = self._fetch_df_from_json(tables[1])

        df = stoich_df.merge(lambda_df["lambda_O2"], left_index=True, right_index=True)
        return df

    def _fetch_df_from_json(self, json_data):
        tbl_cols = [info['attribute'] for info in json_data['data']['attributes']]
        return pd.DataFrame.from_dict(json_data['data']['instances'], orient='index', columns=tbl_cols, dtype=np.float)
