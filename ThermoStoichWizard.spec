/*
A KBase module: ThermoStoichWizard
*/

module ThermoStoichWizard {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_ThermoStoichWizard(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

    /* An X/Y/Z style reference */
    typedef string obj_ref;

    typedef structure {
        obj_ref lambda_tbl;
        obj_ref stoich_tbl;
        string vh_cs;
        string vh_o2;
        string workspace_name;
      } LambdaParams;

    /* run_lambda_analysis: perform lambda analysis*/
    funcdef run_lambda_analysis(LambdaParams params) returns (ReportResults output) authentication required;

};
