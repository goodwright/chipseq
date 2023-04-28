#!/usr/bin/env python3

import os
import argparse
import logging
import json

import rich
import rich.console
import rich.logging
import rich.traceback


log = logging.getLogger()
stderr = rich.console.Console(stderr=True)
stdout = rich.console.Console()
rich.traceback.install(console=stderr, width=200, word_wrap=True, extra_lines=1)

params_ignore = [
    "monochrome_logs",
    "multiqc_config",
    "tracedir",
    "validate_params",
    "show_hidden_params",
    "enable_conda",
    "max_cpus",
    "max_memory",
    "max_time",
    "genome",
    "outdir"
]

params_gw_ignore = [
    "takes_genome",
    "description",
    "advanced",
    "name"
]

params_input = [
    "input"
]

params_genome = [
    "fasta",
    "gtf",
    "gff",
    "gene_bed",
    "bwa_index",
    "bowtie_index",
    "star_index",
    "chromap_index"
]

def print_header():
    stderr.print("\n\n", highlight=False)
    stderr.print("██████████████████████████████████████████████████████████████████████████████████████████", highlight=False)
    stderr.print("░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░", highlight=False)
    stderr.print("[white]░░░░░██████╗░░██████╗░░██████╗░██████╗░██╗░░░░██╗██████╗░██╗░██████╗░██╗░░██╗████████╗░░░░[white]", highlight=False)
    stderr.print("[white]░░░░██╔════╝░██╔═══██╗██╔═══██╗██╔══██╗██║░░░░██║██╔══██╗██║██╔════╝░██║░░██║╚══██╔══╝░░░░[white]", highlight=False)
    stderr.print("[white]░░░░██║░░███╗██║░░░██║██║░░░██║██║░░██║██║░█╗░██║██████╔╝██║██║░░███╗███████║░░░██║░░░░░░░[white]", highlight=False)
    stderr.print("[white]░░░░██║░░░██║██║░░░██║██║░░░██║██║░░██║██║███╗██║██╔══██╗██║██║░░░██║██╔══██║░░░██║░░░░░░░[white]", highlight=False)
    stderr.print("[white]░░░░╚██████╔╝╚██████╔╝╚██████╔╝██████╔╝╚███╔███╔╝██║░░██║██║╚██████╔╝██║░░██║░░░██║░░░░░░░[white]", highlight=False)
    stderr.print("[white]░░░░░╚═════╝░░╚═════╝░░╚═════╝░╚═════╝░░╚══╝╚══╝░╚═╝░░╚═╝╚═╝░╚═════╝░╚═╝░░╚═╝░░░╚═╝░░░░░░░[white]", highlight=False)
    stderr.print("░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░", highlight=False)
    stderr.print("██████████████████████████████████████████████████████████████████████████████████████████", highlight=False)
    stderr.print("\n", highlight=False)
    stderr.print(
        f"[grey25] goodwright/pipeline-tools version 0.1dev - [link=https://github.com/goodwright/pipeline-tools]https://github.com/goodwright/pipeline-tools[/]",
        highlight=False,
    )
    stderr.print("\n", highlight=False)
    stderr.print("██████████████████████████████████████████████████████████████████████████████████████████", highlight=False)
    stderr.print("\n\n", highlight=False)

def create_entry(key, data):
    new_data = {}
    new_data["name"] = key
    new_data["required"] = False

    if "description" in data: new_data["description"] = data["description"]
    if "default" in data: new_data["default"] = data["default"]
    if "pattern" in data: new_data["pattern"] = data["pattern"]
    if "type" in data: 
        new_data["type"] = data["type"]
        if new_data["type"] == "integer": new_data["type"] = "number"

    return new_data

def create(args):
    # Init
    log.info("Create Schema")
    abs_path = os.path.abspath(args.source_dir)
    log.debug(abs_path)

    # Setup paths
    nf_schema_path = os.path.join(abs_path, "nextflow_schema.json")
    gw_schema_path = os.path.join(abs_path, "flow", "schema", args.schema_name + ".json")

    # Load nextflow schema
    with open(nf_schema_path, 'r') as f:
        nf_schema = json.load(f)
    nf_schema = nf_schema[ "definitions"]

    # Init load up existing schema if it exists
    gw_schema = {}
    gw_schema["inputs"] = {}
    if os.path.exists(gw_schema_path):
        with open(gw_schema_path, 'r') as f:
            gw_schema = json.load(f)

    if "sample_options" not in gw_schema["inputs"]: 
        gw_schema["inputs"]["sample_options"] = {}
        gw_schema["inputs"]["sample_options"]["name"] = "Sample options"
        gw_schema["inputs"]["sample_options"]["description"] = "Parameters relating to the sample being analysed."
        gw_schema["inputs"]["sample_options"]["advanced"] = False
        gw_schema["inputs"]["sample_options"]["properties"] = {}

    if "genome_options" not in gw_schema["inputs"]: 
        gw_schema["inputs"]["genome_options"] = {}
        gw_schema["inputs"]["genome_options"]["name"] = "Genome options"
        gw_schema["inputs"]["genome_options"]["description"] = "The genome being aligned to."
        gw_schema["inputs"]["genome_options"]["advanced"] = False
        gw_schema["inputs"]["genome_options"]["takes_genome"] = True
        gw_schema["inputs"]["genome_options"]["properties"] = {}

    # Get a list of all existing params to make sure we dont include params that have been moved etc.
    existing_gw_params = []
    for key, value in gw_schema["inputs"].items():
        if "properties" in value:
            section_params = value["properties"]
            existing_gw_params.extend(section_params)
    existing_gw_params = [item for item in existing_gw_params if item not in params_gw_ignore]
    log.debug(existing_gw_params)

    # Get a list of all valid params in the nf-core schema
    existing_nf_params = []
    for key, value in nf_schema.items():
        section_params = value["properties"]
        for key2, value2 in section_params.items():
            if key2 not in params_ignore:
                value_keys = value2.keys()
                if "hidden" in value_keys and value2["hidden"] == "false":
                    existing_nf_params.append(key2)
                elif "hidden" not in value_keys:
                    existing_nf_params.append(key2)

    # Check for GW params that dont exist in the nf-core schema
    not_in_nf = list(set(existing_gw_params) - set(existing_nf_params))
    for item in not_in_nf:
        log.warning(f"Parameter does not exist in nf-core schema: {item}")

    # Cycle through nf params and include only if it doesnt already exist
    for key, value in nf_schema.items():
        section_params = value["properties"]
        for key2, value2 in section_params.items():
            process = False
            if key2 not in params_ignore:
                value_keys = value2.keys()
                if "hidden" in value_keys and value2["hidden"] == "false":
                    process = True
                elif "hidden" not in value_keys:
                    process = True

            if key2 in existing_gw_params:
                process = False

            if process == True:
                new_entry = create_entry(key2, value2)
                if "pattern" in new_entry: new_entry["type"] = "file"
                if key2 in params_input:
                    gw_schema["inputs"]["sample_options"]["properties"][key2] = new_entry
                    log.info(f"Adding new parameter: sample_options/{key2}")
                elif key2 in params_genome:
                    new_entry["genome_file"] = key2
                    gw_schema["inputs"]["genome_options"]["properties"][key2] = new_entry
                    log.info(f"Adding new parameter: genome_options/{key2}")
                else:
                    if key not in gw_schema["inputs"]:
                        gw_schema["inputs"][key] = {}
                        gw_schema["inputs"][key]["properties"] = {}
                    gw_schema["inputs"][key]["properties"][key2] = new_entry
                    log.info(f"Adding new parameter: {key}/{key2}")

    # Write data
    with open(gw_schema_path, 'w') as f:
        json.dump(gw_schema, f, indent=4)


def scan(args):
    pass

if __name__ == "__main__":
    # Create command args
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    subparsers.required = True

    # Create schema
    convert_args = subparsers.add_parser('create')
    convert_args.set_defaults(func=create)
    convert_args.add_argument('--source-dir', required=True)
    convert_args.add_argument('--schema_name', required=True)

    # Scan schema
    # convert_args = subparsers.add_parser('scan')
    # convert_args.set_defaults(func=scan)
    # convert_args.add_argument('--target-dir', required=True)

    # Setup logging
    log.setLevel(logging.DEBUG)
    log.addHandler(
        rich.logging.RichHandler(
            level=logging.DEBUG,
            console=rich.console.Console(stderr=True),
            show_time=False,
            show_path=True,
            markup=True,
        )
    )

    # Init logging
    print_header()

    # Parse args
    parsed_args = parser.parse_args()

    # Call functions
    parsed_args.func(parsed_args)
