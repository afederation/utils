import pandas as pd


def load_datasets(data_dir='/Users/afederation/data/paris/', error=False):
    df_redcap = pd.read_csv(data_dir + 'redcap_data.csv', sep='\t',
                            error_bad_lines=error, warn_bad_lines=error)
    df_drugs = pd.read_csv(data_dir + 'drug_hts_data.csv', sep='\t',
                           error_bad_lines=error, warn_bad_lines=error)
    df_biomarkers = pd.read_csv(data_dir + 'patient_biomarker_table.csv',
                                error_bad_lines=error, warn_bad_lines=error)

    return df_redcap, df_drugs, df_biomarkers


def add_spm3_column(data):
    data.loc[(data['standard_zscore'] <= -1.25), 'pesronalization'] = 'high'
    data.loc[(data['standard_zscore'] <= -0.75) & (data['standard_zscore'] > -1.25), 'pesronalization'] = 'good'
    data.loc[(data['standard_zscore'] <= -0.5) & (data['standard_zscore'] > -0.75), 'pesronalization'] = 'moderate'
    data.loc[(data['standard_zscore'] <= 0) & (data['standard_zscore'] > -0.5), 'pesronalization'] = 'low'
    data.loc[(data['standard_zscore'] > 0), 'pesronalization'] = 'none'

    data.loc[(data['norm_auc'] <= 0.3), 'sensitivity'] = 'high'
    data.loc[(data['norm_auc'] <= 0.5) & (data['norm_auc'] > 0.3), 'sensitivity'] = 'good'
    data.loc[(data['norm_auc'] <= 0.75) & (data['norm_auc'] > 0.5), 'sensitivity'] = 'moderate'
    data.loc[(data['norm_auc'] <= 0.85) & (data['norm_auc'] > 0.75), 'sensitivity'] = 'low'
    data.loc[(data['norm_auc'] > 0.85), 'sensitivity'] = 'none'

    data.loc[(data['pesronalization'] == 'high' ) & (data['sensitivity'] == 'high'), 'spm3'] = 16
    data.loc[(data['pesronalization'] == 'high' ) & (data['sensitivity'] == 'good'), 'spm3'] = 15
    data.loc[(data['pesronalization'] == 'high' ) & (data['sensitivity'] == 'moderate'), 'spm3'] = 14
    data.loc[(data['pesronalization'] == 'high' ) & (data['sensitivity'] == 'low'), 'spm3'] = 13
    data.loc[(data['pesronalization'] == 'high' ) & (data['sensitivity'] == 'none'), 'spm3'] = 1

    data.loc[(data['pesronalization'] == 'good' ) & (data['sensitivity'] == 'high'), 'spm3'] = 14
    data.loc[(data['pesronalization'] == 'good' ) & (data['sensitivity'] == 'good'), 'spm3'] = 13
    data.loc[(data['pesronalization'] == 'good' ) & (data['sensitivity'] == 'moderate'), 'spm3'] = 12    
    data.loc[(data['pesronalization'] == 'good' ) & (data['sensitivity'] == 'low'), 'spm3'] = 11
    data.loc[(data['pesronalization'] == 'good' ) & (data['sensitivity'] == 'none'), 'spm3'] = 1

    data.loc[(data['pesronalization'] == 'moderate' ) & (data['sensitivity'] == 'high'), 'spm3'] = 13
    data.loc[(data['pesronalization'] == 'moderate' ) & (data['sensitivity'] == 'good'), 'spm3'] = 12
    data.loc[(data['pesronalization'] == 'moderate' ) & (data['sensitivity'] == 'moderate'), 'spm3'] = 11    
    data.loc[(data['pesronalization'] == 'moderate' ) & (data['sensitivity'] == 'low'), 'spm3'] = 10
    data.loc[(data['pesronalization'] == 'moderate' ) & (data['sensitivity'] == 'none'), 'spm3'] = 1  

    data.loc[(data['pesronalization'] == 'low' ) & (data['sensitivity'] == 'high'), 'spm3'] = 12
    data.loc[(data['pesronalization'] == 'low' ) & (data['sensitivity'] == 'good'), 'spm3'] = 11
    data.loc[(data['pesronalization'] == 'low' ) & (data['sensitivity'] == 'moderate'), 'spm3'] = 10
    data.loc[(data['pesronalization'] == 'low' ) & (data['sensitivity'] == 'low'), 'spm3'] = 9
    data.loc[(data['pesronalization'] == 'low' ) & (data['sensitivity'] == 'none'), 'spm3'] = 1 

    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'high') & (data['standard_zscore'] >= 1.25), 'spm3'] = 5
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'good') & (data['standard_zscore'] >= 1.25), 'spm3'] = 4
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'moderate') & (data['standard_zscore'] >= 1.25), 'spm3'] = 3
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'low') & (data['standard_zscore'] >= 1.25), 'spm3'] = 2
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'none') & (data['standard_zscore'] >= 1.25), 'spm3'] = 1  

    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'high') & (data['standard_zscore'] >= 0.75), 'spm3'] = 6
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'good') & (data['standard_zscore'] >= 0.75), 'spm3'] = 5
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'moderate') & (data['standard_zscore'] >= 0.75), 'spm3'] = 4
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'low') & (data['standard_zscore'] >= 0.75), 'spm3'] = 3
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'none') & (data['standard_zscore'] >= 0.75), 'spm3'] = 1   

    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'high') & (data['standard_zscore'] >= 0.5), 'spm3'] = 7
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'good') & (data['standard_zscore'] >= 0.5), 'spm3'] = 6
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'moderate') & (data['standard_zscore'] >= 0.5), 'spm3'] = 5
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'low') & (data['standard_zscore'] >= 0.5), 'spm3'] = 4
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'none') & (data['standard_zscore'] >= 0.5), 'spm3'] = 1    

    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'high') & (data['standard_zscore'] >= 0), 'spm3'] = 8
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'good') & (data['standard_zscore'] >= 0), 'spm3'] = 7
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'moderate') & (data['standard_zscore'] >= 0), 'spm3'] = 6
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'low') & (data['standard_zscore'] >= 0), 'spm3'] = 5
    data.loc[(data['pesronalization'] == 'none' ) & (data['sensitivity'] == 'none') & (data['standard_zscore'] >= 0), 'spm3'] = 1   

    return data


def create_drug_name_to_id_dict(fn='/Users/afederation/data/paris/annotated_complete_drug_library_latest.xlsx'):
    df_drug_annot = pd.read_excel(fn)
    drug_name_to_id_dict = {}

    for index, row in df_drug_annot.iterrows():
        drug_id = row.id_drug_sengine
        drug_names = row.drug_name_alternates.split(', ')
        drug_names.append(row.drug_name_commercial)
        drug_names.append(row.drug_name_display)
        drug_names.append(row.drug_name_biomedtracker)
        if pd.notna(row.drug_name_alternates_sengine):
            for x in row.drug_name_alternates_sengine.split(', '):
                drug_names.append(x)

        drug_names_lower = []
        for drug in list(set(drug_names)):
            if pd.notna(drug):
                drug_names_lower.append(drug.lower())

        for name in drug_names_lower:
            drug_name_to_id_dict[name] = drug_id

    return drug_name_to_id_dict