# Pelops

Dedicated caller for *DUX4* rearrangements from whole genome sequencing data.

## Citing
Pelops is based on a method first described in:

Ryan, S.L., Peden, J.F., Kingsbury, Z. et al.
Whole genome sequencing provides comprehensive genetic testing in childhood B-cell acute lymphoblastic leukaemia.
Leukemia 37, 518–528 (2023). https://doi.org/10.1038/s41375-022-01806-8

Pelops itself is described and validated in:

Grobecker, P., Mijuskovic, M., et al.
Pelops: A dedicated caller for *DUX4* rearrangements from short-read whole genome sequencing data.
In preparation (2024)


## Prerequisites

- Python 3.7 and above.

## Installation

You can install the latest stable released version of Pelops using pip

```shell
pip install ilmn-pelops --upgrade
```

Note: pip/pypi will prefer stable versions (i.e. 0.7.0) over subsequent beta releases (i.e.
0.8.0b1). If you need to install a beta version for testing, please first
uninstall your current version of pelops.

```shell
pip uninstall ilmn-pelops
```

## Usage

Pelops is a tool with a command line interface (cli). Discover its usage with
```shell
pelops --help
```

### Calling _DUX4_-rearrangements
To call DUX4-rearrangements from a BAM/CRAM file, use the `dux4r` subcommand. To
see all available options run
```shell
pelops dux4r --help
```

## Inputs

The input to Pelops is a short-read whole-genome sequencing BAM or CRAM file from a tumour sample,
aligned to the GRCh38 reference genome. The BAM/CRAM file needs to be indexed.
Pelops was tested on alignments by DRAGEN (version 4.0.3), bwa (version 0.7.17), and Isaac (version SAAC01325.18.01.29).

### Systematic noise BED file

To increase specificity when calling non-_IGH_ _DUX4_-rearrangements, we recommend using a systematic noise BED file.
This file contains genomic regions that will be ignored by Pelops when identifying candidate regions involved in a
_DUX4_-rearrangement. Since such regions can be specific to the read alignment tool, reference genome,
sequencing protocol, and cancer type analysed, we recommend creating a separate systematic noise BED file
for each project. One way to obtain these genomic regions would be to run Pelops on a panel of normal samples, which are
guaranteed to have no _DUX4_-rearrangements, and generate a list of false-positive calls.

## Outputs

Pelops outputs results in a JSON file, and optionally exports supporting reads in SAM files.

### Description of JSON output
The **top level** of the JSON contains information about pelops (assumed genome reference, version, name, and CLI command).
It also contains information about the input file (number of unique and mapped reads - which can be a user input).
Finally, it contains a list of rearrangements investigated by pelops.

```json
{
    "reference": "GRCh38",
    "unique_mapped_reads": 1000000000,
    "rearrangements": [...],
    "program_name": "pelops",
    "version": "0.5.0",
    "cli_command": "pelops dux4r --total-number-reads 1000000000 --export . test.bam"
}
```

The **rearrangements** consist of a unique ID, genomic region sets "A" and "B", and the evidence for the rearrangement
between these two regions.
For the command `pelops dux4r`, ID `01` always corresponds to rearrangements between the core _DUX4_ regions and _IGH_,
while ID `02` corresponds to rearrangements of the extended _DUX4_ regions with _IGH_.
IDs `03` and beyond are potential rearrangements of the core _DUX4_ region with other genomic regions (marked as `UNNAMED`);
there can be a variable number of them.

```json
{
  "rearrangements": [
    {
      "id": "01",
      "A": {"name": "CoreDUX4"...},
      "B": {"name": "IGH"...},
      "evidence": {...}
    },
    {
      "id": "02",
      "A": {"name": "ExtendedDUX4"...},
      "B": {"name": "IGH"...},
      "evidence": {...}
    },
    {
      "id": "03",
      "A": {"name": "CoreDUX4"...},
      "B": {"name": "UNNAMED"...},
      "evidence": {...}
    }
  ]
}
```

**A** and **B** document the exact set of genomic regions used for each rearrangement.
For example, the core _DUX4_ region is shown below.
While `IGH`, `CoreDUX4` and `ExtendedDUX4` are pre-defined, each `UNNAMED` region will be different.

```json
{
  "name": "CoreDUX4",
  "regions": [
    {
      "chrom": "chr4",
      "start": 190020407,
      "end": 190023665
    },
    {
      "chrom": "chr4",
      "start": 190066935,
      "end": 190093279
    },
    {
      "chrom": "chr4",
      "start": 190172774,
      "end": 190176845
    },
    {
      "chrom": "chr10",
      "start": 133663429,
      "end": 133685936
    },
    {
      "chrom": "chr10",
      "start": 133739606,
      "end": 133762125
    }
  ]
}
```

The **evidence** for each rearrangement consists of the number of split and paired reads between region sets A and B,
and the spanning read pairs per billion (SRPB). It is calculated as
$$\text{SRPB} = 10^9 \frac{\text{paired reads} + \text{split reads}}{\text{total unique and mapped reads}}.$$

```json
{
  "paired_reads": 15,
  "split_reads": 4,
  "SRPB": 19.0
}
```

### SAM files
Optionally, for each rearrangement a SAM file can be exported which contains all paired and split reads with their mates.
The naming convention is `<id>_<name_A>-<name_B>.sam`, where `<id>`, `<name_A>`, `<name_B>` correspond to the ID and
names of genomic region sets A and B, respectively, as documented in the JSON.

## Contributing
We are not accepting pull requests into this repository at this time, as the
licence currently does not allow modifications by third parties.
For any bug report / recommendation / feature request, please open an issue.

## Credits

See [Authors](AUTHORS.md).
