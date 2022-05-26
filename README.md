# User Guide to Oncomine Automation
## Workflow from run initiation to review of annotated variants
* Annotate variants from oncomine annotator
* Focus on variants of interest
* Visualize and interpret variants
* Link variants to relevant clinical evidence

## Software required:
## 1. [python3](https://www.python.org/download/releases/3.0/)
Python is a high-level, interpreted, general-purpose programming language. It is used to execute automate tasks and conduct data analysis.
## 2. [watchdog](https://pypi.org/project/watchdog/)
Python API library and shell utilities to monitor file system events.
## 3.  [ion_automation](https://github.com/khzhu/ion_automation)
The automated workflow to generate quality control plots/reports and select variants from Oncomine Focus and Myeloid Assay NGS sequencing data.

## Workflow Diagram
![This is a flowchart](https://github.com/khzhu/ion_automation/blob/main/docs/oncomine-workflow.png)

## Directory Structure

```
|— amplicon.dropout.dropoff 
|— downloads
|- dropoff
|- reports
```