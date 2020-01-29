import click
import requests
import pprint

# the pattern used in the Project Title in EVA submission to identify IMAGE related genotype data
IMAGE_TAG = 'Recombination'
SPECIMEN = 'specimen'
ANALYSIS = 'analysis'
# the column list to retrieve data from ENA
FIELD_LIST = [
    "analysis_accession", "study_accession", "sample_accession", "analysis_title", "analysis_type", "center_name",
    "first_public", "last_updated", "study_title", "tax_id", "scientific_name", "analysis_alias", "submitted_bytes",
    "submitted_md5", "submitted_ftp", "submitted_aspera", "submitted_galaxy", "submitted_format",
    "broker_name", "pipeline_name", "pipeline_version", "assembly_type", "accession", "description", "germline"
]


@click.command()
@click.option(
    '--result_type',
    default='specimen',
    help='Specify the result type which can be either specimen (default) or analysis'
)
def main(result_type):
    """
    Main function that will import eva data
    :return: None
    """
    if result_type != SPECIMEN and result_type != ANALYSIS:
        print(f"The parameter result_type can only be '{SPECIMEN}' or '{ANALYSIS}")
        exit(1)
    print("Start importing IMAGE genotype data")
    url = f"http://www.ebi.ac.uk/eva/webservices/rest/v1/meta/studies/all"
    data = requests.get(url).json()
    image_datasets = list()
    for record in data['response'][0]['result']:
        if record['name'].startswith(IMAGE_TAG):
            image_datasets.append(record['id'])

    print(f"There are {len(image_datasets)} IMAGE datasets found within EVA")
    field_str = ",".join(FIELD_LIST)
    results = dict()

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

        # query based on study accession
        optional_str = f"query=study_accession%3D%22{study_accession}%22"
        url = f"{base_url}{optional_str}"
        response = requests.get(url)
        if response.status_code == 204:  # 204 is the status code for no content => the current term does not have match
            continue
        data = response.json()
        displayed = set()
        for record in data:
            if result_type == ANALYSIS:
                results = parse_into_analysis(record, results, eva_summary)
            elif result_type == SPECIMEN:
                results = parse_into_specimen(record, results)
            count = len(results)
            if count % 50 == 0 and str(count) not in displayed:
                displayed.add(str(count))
                print(f"Processed {count} records")
        # end of analysis list for one study loop
    # end of all studies loop

    # consume the parsing results, currently just print out
    pprint.pprint(results)


def parse_into_specimen(record, results):
    """
    Parse the EVA record stored within ENA into specimen based
    :param record: ENA API record
    :param results: existing result
    :return: updated result
    """
    sample_accession = record['sample_accession']
    if sample_accession in results:  # the sample used in more than one analysis
        es_doc = results[sample_accession]
    else:
        file_server_type = determine_file_type(record)
        if len(file_server_type) == 0:
            return results
        es_doc = extract_files(record, file_server_type)
        if not es_doc:
            return results

    es_doc.setdefault('analyses', list())
    es_doc['analyses'].append(record['analysis_accession'])
    results[sample_accession] = es_doc

    return results


def extract_files(record, file_server_type):
    """
    Extract file information from the ENA API result
    :param record: ENA API record
    :param file_server_type: the system where the files are stored
    :return: the dict containing file information
    """
    es_doc = dict()
    files = record[f"submitted_{file_server_type}"].split(";")
    sizes = record["submitted_bytes"].split(";")
    formats = record["submitted_format"].lower().split(";")
    # for ENA, it is fixed to MD5 as the checksum method
    checksums = record["submitted_md5"].split(";")
    if len(files) != len(checksums) or len(files) != len(sizes) or len(files) != len(formats) or len(files) == 0:
        return es_doc
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
    return es_doc


def parse_into_analysis(record, results, eva_summary):
    """
    Parse the EVA record stored within ENA and EVA into analysis based
    :param record: ENA API record
    :param results: existing result
    :param eva_summary: the EVA API result
    :return: updated result
    """
    analysis_accession = record['analysis_accession']
    if analysis_accession in results:
        es_doc = results[analysis_accession]
    else:
        es_doc = convert_analysis(record)
        if not es_doc:
            return results
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
    results[analysis_accession] = es_doc
    return results


def determine_file_type(record):
    """
    Determine the file type used in the ENA archive
    :param record: ENA API record
    :return: the file type
    """
    file_server_types = ['ftp', 'galaxy', 'aspera']
    for tmp in file_server_types:
        key_to_check = f"submitted_{tmp}"
        if key_to_check in record and record[key_to_check] != '':
            return tmp
    return ''


def convert_analysis(record):
    """
    Convert the ENA API result into an analysis object
    :param record: ENA API result
    :return: analysis object
    """
    file_server_type = determine_file_type(record)
    if len(file_server_type) == 0:
        return dict()

    es_doc = extract_files(record, file_server_type)
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
