import json
import os
import shutil
from typing import Dict, List
import logging
import requests
from Parser1 import Parser1
from Helperfuncs import Helperfuncs as helpers
import pandas as pd
import csv


info_file='templets/weekly_run_parameters_table.csv'
data = {}
with open(info_file, mode='r') as file:
    csv_reader = csv.reader(file, delimiter=',')
    for row in csv_reader:
        if row:  # Avoid empty rows
            key, value = row
            data[key] = value
print(data)

output_folder = data['outfolder']
batch_name = data['batch_id']
genome_build = data['genome_build']
logging.basicConfig(
    level=logging.ERROR,   # Set the minimum level to DEBUG
    format='%(message)s',  # Custom format
    filename=f'{output_folder}/weekly_run.log',    # Output to a file
)

logger = logging.getLogger(__name__)
api_geneyx = {"apiUserKey": "PfiiZhaPOjywmersy9uI6Zma/XmSZGGYQwRP8+LhJrmm52aecPdhzA==", 
              "apiUserId": "HdpXjtnwUF0xG+EcNXuGberA752xUDKVLkER39QKcm8=" }
protocol_dict={}
with open("templets/protocols_dict_add.tsv", "r") as f:
    csv_reader = csv.DictReader(f, delimiter='\t')
    for row in csv_reader:
        protocol_key = row['protocol_key'].lower()
        print(protocol_key)
        print(row[genome_build])
        protocol_dict[protocol_key] = row[genome_build]
    

print(protocol_dict)
batch_patient_dict={}
batch_post_pattern = "https://analysis.geneyx.com/api/BatchVcfSamples?batchName="

batch_call_post_url = batch_post_pattern + batch_name
batch_data_response = requests.post(batch_call_post_url, data = api_geneyx)

all_patients_json = batch_data_response.json()["Data"]

for i in all_patients_json:
    flag = "YES"
    batch_patient_dict[i["Patient"]] = [i["SerialNumber"],flag]
    

def main(project_table_data, file_data_parser, file_log):
    parser = Parser1(file_data_parser, file_log, batch_patient_dict,protocol_dict,g_build=genome_build)
    cases_api_df = pd.DataFrame()
    problem_cases=pd.DataFrame()
     # Convert specific columns to lowercase
    cols_to_lower = ['geneyx protocol', 'single/duo/trio/quatro/panel/mtDNA', 'Examinee']
    project_table_data[cols_to_lower] = project_table_data[cols_to_lower].apply(lambda col: col.str.lower())

    # Parse single 
    single_df = project_table_data.loc[project_table_data['single/duo/trio/quatro/panel/mtDNA'] == 'single']
    for _, proband_row in single_df.iterrows():
        if not proband_row.empty:
            api_data=parser.parse_row(proband_row)
            if not api_data: 
                logger.error(f"Error parsing row: {ex_number}")
                problem_cases=pd.concat([problem_cases, pd.DataFrame([proband_row])], ignore_index=True)
            else:
                for case in api_data:
                    cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)
            


    print(len(cases_api_df))
    # Parse duo
    duo_df = project_table_data.loc[project_table_data['single/duo/trio/quatro/panel/mtDNA'] == 'duo']
    uniq_hpo_duo = set(duo_df.EX)

    for ex_number in uniq_hpo_duo:
        ex_number = str(ex_number)
        ex_specific = duo_df.loc[duo_df['EX'] == ex_number]
        mother_row=ex_specific.loc[ex_specific['Examinee'] == "spouse#1"].squeeze()
        mother_Sn=helpers.get_full_sampleSn(mother_row['DNA number'],batch_patient_dict,genome_build)
        if not mother_Sn:
            logger.error(f"not found sampleSn for {mother_row['DNA number']}")
            continue
        father_row=ex_specific.loc[ex_specific['Examinee'] == "spouse#2"].squeeze()
        father_Sn=helpers.get_full_sampleSn(father_row['DNA number'],batch_patient_dict,genome_build)
        if not father_Sn:
            logger.error(f"not found sampleSn for {father_row['DNA number']}")
            problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
        if not mother_row.empty and not father_row.empty:
            AssociatedSamples= [
            {"SampleId": mother_Sn, "Relation": "mother", "Affected": "Unaffected"},
            {"SampleId": father_Sn, "Relation": "father", "Affected": "Unaffected"}]
            api_data=parser.parse_row(mother_row, AssociatedSamples)
            if not api_data: 
                logger.error(f"Error parsing row: {ex_number}")
                problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
                continue
            for case in api_data:
                cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)

        
    # Parse trio
    trio_df = project_table_data.loc[(project_table_data['single/duo/trio/quatro/panel/mtDNA'] == 'trio') | (project_table_data['single/duo/trio/quatro/panel/mtDNA'] == 'quatro')]
    uniq_hpo_trio = set(trio_df.EX)
    for ex_number in uniq_hpo_trio:
        AssociatedSamples, AssociatedSamples2 = [],[]
        
        proband_row=None
        ex_number = str(ex_number)
        ex_specific = trio_df.loc[trio_df['EX'].astype(str) == ex_number]
        # if ex_specific.empty:
        #     ex_specific = trio_df.loc[trio_df['EX'].astype(str) == ex_number]
        print(ex_specific[['Examinee', 'DNA number']])
        
        proband_row=ex_specific.loc[(ex_specific['Examinee'] == "proband") | (ex_specific['Examinee'] == 'proband1')].squeeze()
        if proband_row.empty:
            logger.error(f"not proband found for {ex_number} exome number")
            problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
            continue
        proband1_Sn = helpers.get_full_sampleSn(proband_row['DNA number'],batch_patient_dict,genome_build)

        mother_row=ex_specific.loc[ex_specific['Examinee'] == "mother"].squeeze()

        father_row=ex_specific.loc[ex_specific['Examinee'] == "father"].squeeze()

        
        bool_father, bool_mother = False, False
        if not mother_row.empty:
            bool_mother=True
            mother_Sn=helpers.get_full_sampleSn(mother_row['DNA number'],batch_patient_dict,genome_build)
            AssociatedSamples.append({"SampleId": mother_Sn, "Relation": "mother", "Affected": "Unaffected"})
        
        if not father_row.empty:
            bool_father=True
            father_Sn=helpers.get_full_sampleSn((father_row['DNA number']),batch_patient_dict,genome_build)
            AssociatedSamples.append({"SampleId": father_Sn, "Relation": "father", "Affected": "Unaffected"})
        if not any([bool_father, bool_mother]):
            logger.error(f"not found mother or father for {ex_number} exome number")
            continue
        proband2_row=ex_specific.loc[ex_specific['Examinee'] == "proband2"].squeeze()

        if not proband2_row.empty:
            sibling_Sn=helpers.get_full_sampleSn((proband2_row['DNA number']),batch_patient_dict,genome_build)
            AssociatedSamples.append({"SampleId": sibling_Sn, "Relation": "sibling", "Affected": "Unknown"})
            AssociatedSamples2.append({"SampleId": proband1_Sn, "Relation": "sibling", "Affected": "Unknown"})
            proband2_r2 = proband2_row.copy()
            proband2_r2['EX'] = f"{ex_number}_3"
            api_data=parser.parse_row(proband2_r2, AssociatedSamples2)
            if not api_data: 
                logger.error(f"Error parsing row: {ex_number}")
                problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
                continue
            for case in api_data:
                cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)
            
        if AssociatedSamples:
            api_data=parser.parse_row(proband_row, AssociatedSamples)
            if not api_data: 
                logger.error(f"Error parsing row: {ex_number}")
                problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
                continue
            for case in api_data:
                cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)
        else: 
            file_log.write(f"{ex_number} havent mother or father")

        ########### add shred cariers here
        if bool_father and bool_mother:
            mother_r2 = mother_row.copy()

            mother_r2[['EX', 'single/duo/trio/quatro/panel/mtDNA', 'geneyx protocol', 'Phenotype']] = [f"{ex_number}_2", "duo", 'shared carrier screening', "shared carrier screening"]

            
            AssociatedSamples= [
            {"SampleId": mother_Sn, "Relation": "mother", "Affected": "Unaffected"},
            {"SampleId": father_Sn, "Relation": "father", "Affected": "Unaffected"}]
            api_data=parser.parse_row(mother_r2, AssociatedSamples)
            if not api_data: 
                logger.error(f"Error parsing row: {ex_number}")
                problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
            for case in api_data:
                cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)
        elif bool_father:
            father_r2 = mother_row.copy()
            father_r2['EX'] = f"{ex_number}_2"
            father_r2['single/duo/trio/quatro/panel/mtDNA'] = "single"
            father_r2['geneyx protocol'] = 'clinical'
            father_r2['Phenotype'] = "father carrier screening"
            api_data=parser.parse_row(father_r2)
            if not api_data: 
                logger.error(f"Error parsing row: {ex_number}")
                problem_cases=pd.concat([problem_cases, ex_specific], ignore_index=True)
            for case in api_data:
                cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)
        elif bool_mother:
            mother_r2 = mother_row.copy()
            mother_r2['EX'] = f"{ex_number}_2"
            mother_r2['single/duo/trio/quatro/panel/mtDNA'] = "single"
            mother_r2['geneyx protocol'] = 'clinical'
            mother_r2['Phenotype'] = "mother carrier screening"
            api_data=parser.parse_row(mother_r2)
            for case in api_data:
                cases_api_df = pd.concat([cases_api_df, pd.DataFrame([case])], ignore_index=True)
    print(len(cases_api_df))
    col = cases_api_df.pop('AssociatedSamples')
    cases_api_df.insert(4, 'AssociatedSamples', col)
    file_path = output_folder+"/api_commands.csv"
    # if api_commands file exists, move it to a new name and create new one
    if os.path.exists(file_path):
        new_file_path = output_folder+"/api_commands_1.csv"
        shutil.move(file_path, new_file_path)

    cases_api_df.to_csv(output_folder+"/api_commands.csv", index=False)
    problem_cases.to_csv(output_folder+"/problem_cases.csv", index=False)
        
        
if __name__ == "__main__":
    # Load project table data
    project_table_data = pd.read_excel(output_folder+"/hadas.xlsx")

    # Open files
    file_data_parser = open("data_parser.json", "w")
    file_log = open(output_folder+"/log.txt", "w")

    main(project_table_data, file_data_parser, file_log)

    # Close files
    file_data_parser.close()
    file_log.close()
