import click
import requests
import pprint

IMAGE_TAG = 'Recombination'
FIELD_LIST = [
    "analysis_accession", "study_accession", "sample_accession", "analysis_title", "analysis_type", "center_name",
    "first_public", "last_updated", "study_title", "tax_id", "scientific_name", "analysis_alias", "submitted_bytes",
    "submitted_md5", "submitted_ftp", "submitted_aspera", "submitted_galaxy", "submitted_format",
    "broker_name", "pipeline_name", "pipeline_version", "assembly_type", "accession", "description", "germline"
]




@click.command()
# TODO check single or double quotes
def main():
    """
    Main function that will import analysis data from ena
    :return:
    """
    print("Start importing IMAGE genotype data")
    url = f"http://www.ebi.ac.uk/eva/webservices/rest/v1/meta/studies/all"
    data = requests.get(url).json()
    image_datasets = list()
    for record in data['response'][0]['result']:
        if record['name'].startswith(IMAGE_TAG):
            image_datasets.append(record['id'])

    print(f"There are {len(image_datasets)} IMAGE datasets found within EVA")
    field_str = ",".join(FIELD_LIST)
    analyses = dict()

    ena_params = {
        "result": "analysis",
        "format": "JSON",
        "limit": "0",
        "fields": field_str,
        "dataPortal": "ena"
    }
    base_url = "https://www.ebi.ac.uk/ena/portal/api/search/?"
    for k, v in ena_params.items():
        base_url = f"{base_url}{k}={v}&"

    for study_accession in image_datasets:
        url = f"http://www.ebi.ac.uk/eva/webservices/rest/v1/studies/{study_accession}/summary"
        # expect always to have data from EVA as the list is retrieved live
        # get EVA summary
        eva_summary = requests.get(url).json()['response'][0]['result'][0]

        # f"https://www.ebi.ac.uk/ena/portal/api/search/?result=analysis&format=JSON&limit=0&" \
        #    f"query=study_accession%3D%22{study_accession}%22&fields={field_str}"
        # extra constraint based on study accession
        optional_str = f"query=study_accession%3D%22{study_accession}%22"
        url = f"{base_url}{optional_str}"
        response = requests.get(url)
        if response.status_code == 204:  # 204 is the status code for no content => the current term does not have match
            continue
        data = response.json()
        displayed = set()
        for record in data:
            analysis_accession = record['analysis_accession']
            if analysis_accession in analyses:
                es_doc = analyses[analysis_accession]
            else:
                es_doc = convert_analysis(record)
                if not es_doc:
                    continue
                # in ENA api, it is description in ena result, different to analysis_description in faang result portal
                es_doc['description'] = record['description']
                if eva_summary['experimentType'] != '-':
                    es_doc['experimentType'] = eva_summary['experimentType'].split(', ')
                # es_doc['program'] = eva_summary['program']
                if eva_summary['platform'] != '-':
                    es_doc['platform'] = eva_summary['platform'].split(', ')
                # imputation has not been exported in the ENA warehouse
                # use PRJEB22988 (non farm animal) as example being both imputation and phasing project
                # es_doc['imputation'] = record['imputation']
            es_doc['sampleAccessions'].append(record['sample_accession'])
            analyses[analysis_accession] = es_doc
            count = len(analyses)
            if count % 50 == 0 and str(count) not in displayed:
                displayed.add(str(count))
                logger.info(f"Processed {count} analysis records")
        # end of analysis list for one study loop
    # end of all studies loop
    pprint.pprint(analyses)


def convert_analysis(record):
    file_server_types = ['ftp', 'galaxy', 'aspera']
    file_server_type = ''
    for tmp in file_server_types:
        key_to_check = f"submitted_{tmp}"
        if key_to_check in record and record[key_to_check] != '':
            file_server_type = tmp
            break
    if len(file_server_type) == 0:
        return dict()

    es_doc = dict()
    files = record[f"submitted_{file_server_type}"].split(";")
    sizes = record["submitted_bytes"].split(";")
    formats = record["submitted_format"].lower().split(";")
    # for ENA, it is fixed to MD5 as the checksum method
    checksums = record["submitted_md5"].split(";")
    if len(files) != len(checksums) or len(files) != len(sizes) or len(files) != len(formats) or len(files) == 0:
        return dict()
    for i, file in enumerate(files):
        fullname = file.split("/")[-1]
        # filename = fullname.split(".")[0]
        suffix = fullname.split(".")[-1]
        if suffix != 'md5':
            es_doc.setdefault('fileNames', list())
            es_doc.setdefault('fileTypes', list())
            es_doc.setdefault('fileSizes', list())
            es_doc.setdefault('checksumMethods', list())
            es_doc.setdefault('checksums', list())
            es_doc.setdefault('urls', list())
            es_doc['fileNames'].append(fullname)
            es_doc['fileTypes'].append(formats[i])
            es_doc['fileSizes'].append(convert_readable(sizes[i]))
            es_doc['checksumMethods'].append('md5')
            es_doc['checksums'].append(checksums[i])
            es_doc['urls'].append(file)
    es_doc['accession'] = record['analysis_accession']
    es_doc['title'] = record['analysis_title']
    es_doc['alias'] = record['analysis_alias']
    es_doc['releaseDate'] = record['first_public']
    es_doc['updateDate'] = record['last_updated']
    es_doc.setdefault('organism', dict())
    es_doc['organism']['text'] = record['scientific_name']
    es_doc['organism']['ontologyTerms'] = f"http://purl.obolibrary.org/obo/NCBITaxon_{record['tax_id']}"

    es_doc['datasetAccession'] = record['study_accession']
    es_doc.setdefault('sampleAccessions', list())

    es_doc['analysisCenter'] = record['center_name']
    es_doc['analysisType'] = record['analysis_type']
    return es_doc


def convert_readable(size_to_convert):
    """
    This function will convert size to human readable string
    :param size_to_convert: size in bytes
    :return: human-readable string with size
    """
    i = 0
    size_to_convert = int(size_to_convert)
    units = ['B', 'kB', 'MB', 'GB', 'TB', 'PB']
    for i in range(6):
        size_to_convert /= 1024
        if size_to_convert < 1:
            break
    size_to_convert *= 1024
    if i == 0:
        return f"{size_to_convert}B"
    return f"{round(size_to_convert, 2)}{units[i]}"


if __name__ == "__main__":
    main()
