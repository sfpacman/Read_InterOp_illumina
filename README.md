<div id="top"></div>




<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#configuration">Configuration</a></li>
  </ol>
</details>
<!-- ABOUT THE PROJECT -->

## Contributors
[@sfpacman](https://github.com/sfpacman)

## About The Project

This repository contains a standalone Python script designed to extract information from Illumina binary QC files and convert to a YAML file. The script serves as a refactored replacement for [illuminate](https://bitbucket.org/invitae/illuminate) module and uses [InterOp](https://illumina.github.io/interop/index.html) module, which was subsequently integrated into the whole-exome sequencing (WES) pipeline during my tenure at UCSF.

## Prerequisites

conda can be used to install ```interop``` and ```pandas``` . 

## Installation
1. Clone the repo
   ```sh
   git clone https://github.com/lawrenson-lab/Read_InterOp_illumina/
   ```
2. Install packages via conda
    ```sh
    conda install bioconda::illumina-interop
    conda install pandas
    ```
You are now ready to run the script!

<!-- USAGE EXAMPLES -->
## Usage
Execute the Python script in the terminal:
  ```sh
   python run_qc_yaml_interop_production.py <target_dir> <out_dir>
   ```
### Input
Provide a directory containing `RunInfo.xml` and an `InterOp` subdirectory containing Illumina binary files
### Output
A yaml file contains the following QC metrics: 
* lane_level_metrics
* xread_level_metrics
* read_level_metrics
* read_yield_metrics
* sample_level_metrics
* run_level_metrics

### Configuration
No additional arguments are included for modifying Illumina QC column names and metric conversion for the final report. However, you can simply change the implenetation for the following functions. 

* ```get_columns_name()```
* ```get_metrics()```

Consider implementing a YAML configuration parsing function in the future for enhanced flexibility.
