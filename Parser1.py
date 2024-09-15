import csv
import json
from typing import Dict, List
from Helperfuncs import Helperfuncs as helpers
import pandas as pd
import logging
logger = logging.getLogger(__name__)

# file_log = open("log.txt", "a")


# with open("templets/protocols_dict.csv", "r") as f:
#     reader = csv.reader(f)
#     protocol_dict = {row[0]: row[1] for row in reader}

api_geneyx = {"apiUserKey": "PfiiZhaPOjywmersy9uI6Zma/XmSZGGYQwRP8+LhJrmm52aecPdhzA==", 
              "apiUserId": "HdpXjtnwUF0xG+EcNXuGberA752xUDKVLkER39QKcm8=" }
class Parser1:
    def __init__(self, file_data_parser, file_log,batch_dict,protocol_dict,g_build):
        self.file_data_parser = file_data_parser
        self.file_log = file_log
        self.batch_dict = batch_dict
        self.protocol_dict = protocol_dict
        self.g_bluid = g_build

    def parse_row(self, row, assocateds=None):
        if row.empty: return
        print(type(row))
        protocol_key=row['geneyx protocol']
        if protocol_key=='clinical':
            protocol_key='clinical'+f"_{row['single/duo/trio/quatro/panel/mtDNA']}"
        if type(protocol_key)!=str:
            protocol_key=protocol_key.values[0]
        
        protocol_key=protocol_key.lower()
        if protocol_key.endswith(' '):
            protocol_key = protocol_key[:-1]
        
        if not protocol_key in self.protocol_dict.keys():

            logger.warning(f"{row['DNA number']} {row['geneyx protocol']} havent protocol")
            return
        sampleSn=helpers.get_full_sampleSn(row['DNA number'],self.batch_dict,self.g_bluid)
        if not sampleSn: 
            logger.error(f"not found sampleSn for {row['DNA number']}")
            return
        data_dict = {
        "SerialNumber": str(row['EX']),
        "SubjectId": str(row['DNA number']),
        "ProbandSampleId": sampleSn,
        "ProtocolId": self.protocol_dict[protocol_key],
        "Phenotypes": helpers.check_phenotype(row["HPO (terms)"]),
        "Description": helpers.check_phenotype(row["Phenotype"]),
         "Owner": row["Sender"],
        }
        if assocateds: data_dict["AssociatedSamples"]=assocateds
            

        protocols_lst=data_dict['ProtocolId'].split(',')
        if len(protocols_lst)>1:
            for i, protocol in enumerate(protocols_lst):
                data_dict['ProtocolId']=protocol
                data_dict["SerialNumber"]=f"{str(row['EX'])}_{i}"
                yield data_dict
        else:
            yield data_dict
        



    