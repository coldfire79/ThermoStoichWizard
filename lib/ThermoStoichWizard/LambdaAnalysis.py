import logging
import os
import uuid
import pandas as pd
import numpy as np
import itertools
from scipy.stats import pearsonr

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

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

        html_folder = os.path.join(self.shared_folder, 'html')
        os.mkdir(html_folder)

        vis_content = ''
        X, Y = np.meshgrid(vh_cs, vh_o2)
        corr_mats = {i:np.zeros(X.shape) for i in ['r_lambda_rbiom','r_lambda_ro2','r_lambda_rhco3']}
        print(corr_mats)
        for i in range(len(vh_cs)):
            for j in range(len(vh_o2)):
                vhcs = X[i,j]
                vho2 = Y[i,j]

                df['r_biom'] = mu_max * np.exp(-np.abs(df['donor'])/vhcs)* np.exp(-np.abs(df['acceptor'])/vho2)
                df['r_o2'] = np.abs(df['acceptor']) * df['r_biom']
                df['r_hco3'] = np.abs(df['hco3']) * df['r_biom']

                print(vhcs, vho2)

                fout = 'correlation_plot_vhcs={:.2f}_vho2={:.2f}.png'.format(vhcs, vho2)
                fpath = os.path.join(html_folder, fout)
                r_lambda_rbiom, r_lambda_ro2, r_lambda_rhco3 = self._plot_correlation(df, fout=fpath)

                corr_mats['r_lambda_rbiom'][i,j] = r_lambda_rbiom[0]
                corr_mats['r_lambda_ro2'][i,j] = r_lambda_ro2[0]
                corr_mats['r_lambda_rhco3'][i,j] = r_lambda_rhco3[0]

                vis_content += '<div>'
                vis_content += '<h3>Vh[Cs]={:.2f}, Vh[O2]={:.2f}<h4>'.format(vhcs, vho2)
                vis_content += '<img alt="{0}" src="{0}" style="width: 100%; display: block;">'.format(fout)
                vis_content += '</div>'

        plt.close('all')
        fig = plt.figure(figsize=(12,3))
        for i, (corr, label) in enumerate(zip(['r_lambda_rbiom','r_lambda_ro2','r_lambda_rhco3'],
                                     ["$r_{Biom}$","$r_{O_2}$","$r_{HCO_3^-}$"])):
            
            ax = fig.add_subplot(1, 3, i+1)
            tdf = pd.DataFrame(corr_mats[corr], columns=vh_cs)
            tdf.index = vh_o2
            sns.heatmap(tdf, center=0, cmap="RdBu_r", annot=True, fmt=".2f", ax=ax)
            ax.set_xlabel(r"$V_h[C_s]$")
            ax.set_ylabel(r"$V_h[O_2]$")

            # print(i)
            # Z = corr_mats[corr]
            # # set up the axes for the first plot
            # ax = fig.add_subplot(1, 3, i+1, projection='3d')

            # # plot a 3D surface like in the example mplot3d/surface3d_demo
            # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
            #                        linewidth=0, antialiased=False)
            
            # ax.dist = 10
            # ax.azim = -115
            # ax.set_zlim(-1.01, 1.01)
            # ax.set_xlabel(r"$V_h[C_s]$")
            # ax.set_ylabel(r"$V_h[O_2]$")
            # ax.set_zlabel(r"$\rho$", rotation=180)
            # ax.set_title(label)

            # fig.colorbar(surf, shrink=0.5, aspect=10)

        plt.tight_layout()
        plt.savefig(os.path.join(html_folder, "correlation_3d.png"))
        
        corr_content = '<img alt="{0}" src="{0}" style="width: 100%; display: block;">'.format("correlation_3d.png")

        with open(os.path.join(os.path.dirname(__file__), 'templates', 'lambda_template.html'),
                  'r') as template_file:
            report_html = template_file.read()
            report_html = report_html.replace('<p>Visualization_Content</p>', vis_content)
            report_html = report_html.replace('<p>3D Correlation Plot</p>', corr_content)
            # report_html = report_html.replace('Number of compounds:', 'Number of compounds: {}'.format(num_cpds))
            # report_html = report_html.replace('<!--[Results]-->', html_str)
            
        with open(os.path.join(html_folder, "index.html"), 'w') as index_file:
            index_file.write(report_html)

        report = KBaseReport(self.callback_url)
        html_dir = {
            'path': html_folder,
            'name': 'index.html',
            'description': 'Lambda Analysis Report'
        }
        report_info = report.create_extended_report({
            # # 'objects_created': objects_created,
            # 'file_links': output_files,
            'html_links': [html_dir],
            'direct_html_link_index': 0,
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

    def _plot_correlation(self, df, fout):
        r_lambda_rbiom = pearsonr(df['lambda_O2'], df['r_biom'])
        r_lambda_ro2 = pearsonr(df['lambda_O2'], df['r_o2'])
        r_lambda_rhco3 = pearsonr(df['lambda_O2'], df['r_hco3'])

        # print(r_lambda_rbiom, r_lambda_ro2, r_lambda_rhco3)

        fig, ax = plt.subplots(2,3, figsize=(11,5))
        sns.distplot(df['r_biom'], ax=ax[0,0])
        sns.distplot(df['r_o2'], ax=ax[0,1])
        sns.distplot(df['r_hco3'], ax=ax[0,2])
        ax[0,0].set_xlabel(r"$r_{Biom}$")
        ax[0,1].set_xlabel(r"$r_{O_2}$")
        ax[0,2].set_xlabel(r"$r_{HCO_3^-}$")

        sns.scatterplot(df['lambda_O2'], df['r_biom'], alpha=0.6, ax=ax[1,0])
        sns.scatterplot(df['lambda_O2'], df['r_o2'], alpha=0.6, ax=ax[1,1])
        sns.scatterplot(df['lambda_O2'], df['r_hco3'], alpha=0.6, ax=ax[1,2])
        ax[1,0].set_ylabel(r"$r_{Biom}$")
        ax[1,1].set_ylabel(r"$r_{O_2}$")
        ax[1,2].set_ylabel(r"$r_{HCO_3^-}$")

        ax[1,0].set_title(r"$\rho$={:.3f}".format(r_lambda_rbiom[0]))
        ax[1,1].set_title(r"$\rho$={:.3f}".format(r_lambda_ro2[0]))
        ax[1,2].set_title(r"$\rho$={:.3f}".format(r_lambda_rhco3[0]))

        for i in range(3):
            ax[1,i].set_xlabel(r"$\lambda$")
            
        if fout:
            plt.tight_layout()
            plt.savefig(fout)
        
        return r_lambda_rbiom, r_lambda_ro2, r_lambda_rhco3