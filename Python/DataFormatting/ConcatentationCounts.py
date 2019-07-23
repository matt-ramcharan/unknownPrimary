import numpy as np
import pandas as pd
import os
import requests
import json
import re

fields = [
    "case_id",
    "primary_site",
    "disease_type",
    "diagnoses.vital_status",
    #"files.data_category",
    #"files.analysis.workflow_type"
    "file_name",
    "cases.submitter_id",
    #"cases.samples.sample_type",
    "cases.disease_type"
    #,    "cases.project.project_id"
    ]

fields = ",".join(fields)

cases_endpt = "https://api.gdc.cancer.gov/cases"
files_endpt = "https://api.gdc.cancer.gov/files"

breast_filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.primary_site",
            "value": ["Breast"]
            }
        },
        {
        "op": "in",
            "content":{
                "field": "files.data_category",
                "value": ["Transcriptome Profiling"]
            }
        },
        {
        "op": "in",
            "content":{
                "field": "files.analysis.workflow_type",
                "value": ["HTSeq - Counts"]
            }
        }
    ]
}


colon_filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "primary_site",
            "value": ["Colon"]
            }
        },
        {
        "op": "in",
            "content":{
                "field": "files.data_category",
                "value": ["Transcriptome Profiling"]
            }
        },
        {
        "op": "in",
            "content":{
                "field": "files.analysis.workflow_type",
                "value": ["HTSeq - Counts"]
            }
        }
    ]
}


# A POST is used, so the filter parameters can be passed directly as a Dict object.
params = {
    "filters": breast_filters,
    "fields": fields,
    "format": "TSV",
    "size": "2000"
    }

# The parameters are passed to 'json' rather than 'params' in this case
response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)

print(response.content.decode("utf-8"))

#####
# Here a GET is used, so the filter parameters should be passed as a JSON string.

params = {
    "filters": json.dumps(breast_filters),
    "fields": "file_id",
    "format": "JSON",
    "size": "1000"
    }

response = requests.get(files_endpt, params = params)

file_uuid_list = []

# This step populates the download list with the file_ids from the previous query
for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    file_uuid_list.append(file_entry["file_id"])

data_endpt = "https://api.gdc.cancer.gov/data"

params = {"ids": file_uuid_list}

response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)

