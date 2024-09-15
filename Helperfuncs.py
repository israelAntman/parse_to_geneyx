import requests
import pandas as pd
import logging
logger = logging.getLogger(__name__)
#import lift_over_dx
#from .get_full_sample_name_v2 import get_full_sample_name
#from django.conf.urls import url
api_geneyx = {"apiUserKey": "PfiiZhaPOjywmersy9uI6Zma/XmSZGGYQwRP8+LhJrmm52aecPdhzA==", 
              "apiUserId": "HdpXjtnwUF0xG+EcNXuGberA752xUDKVLkER39QKcm8=" }
url_samples ="https://analysis.geneyx.com/api/samples" 
api_find_samples=api_geneyx | {"pageSize": 10,"page": 0}
file_log = open("log.txt", "a")

class Helperfuncs:
    
    def check_phenotype(hpo):
        
        if issubclass(type(hpo), float):
            hpo = ""
        return hpo

    
    def get_full_sampleSn(dna_num,batch_patient_dict, build):
        
        if type(dna_num) !=str:
                print(dna_num)
                
                if type(dna_num) == pd.Series:
                    if len(dna_num)>0: 
                        dna_num=str(dna_num[0])
                    else:
                        logger.error("empty series")
                        return
                else:
                    try:
                        dna_num=str(dna_num)
                    except:
                        logger.error(f"cant convert {dna_num} to string")
                        return
                print(dna_num)
        sampleSn = ""
        if dna_num in batch_patient_dict.keys():
            sampleSn = batch_patient_dict[dna_num][0]
            batch_patient_dict[dna_num] = [sampleSn, 'YES']

        else:
            
            build_dict={'hg38':['hg38'],'hg19':['grch37', 'hs37d5', 'hg19']}
            build_names=build_dict[build]
            if build not in ('hg38', 'grch37', 'hs37d5','hg19'):
                logger.error ('Error: build number is invalid')
                return

            # if build == 'hg19':
            #     alt_builds = ('grch37', 'hs37d5', 'hg19') #list of aliases for build number
            # else:
            #     alt_builds = ('hg38')

            api_samples = {
                "apiUserKey": "PfiiZhaPOjywmersy9uI6Zma/XmSZGGYQwRP8+LhJrmm52aecPdhzA==",
                "apiUserId": "HdpXjtnwUF0xG+EcNXuGberA752xUDKVLkER39QKcm8=",
                "query": str(dna_num),
                "pageSize": 10,
                "page": 0
            }

            
            data_response = requests.post(url_samples, json=api_samples)
            

            #print(data_response.json())
            full_sample_names = data_response.json()["Data"]
            #print(full_sample_names)
            find=False
            for sample_id in full_sample_names:
                if any(b in sample_id for b in build_names):
                    find=True
                    
                
                    if sample_id.split('.')[0] == dna_num:
                        find=True
                        sampleSn=sample_id
                        batch_patient_dict[dna_num] = [sampleSn, 'YES']
                    elif sample_id.split('_')[0] == dna_num:
                        find=True
                        sampleSn=sample_id
                        batch_patient_dict[dna_num] = [sampleSn, 'YES']
            if not find:
                batch_patient_dict[dna_num] = [sampleSn, 'NO']
                logger.error(f"{dna_num} havent sample in {build} build")
                    
        return sampleSn
            