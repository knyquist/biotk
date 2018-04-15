# Biochemistry toolkit

This repository contains a collection of tools generally useful for performing analysis of single-molecule sequencing data, especially from the perspective of an enzymologist. The tools are of varying scope and the folder structure is intended to reflect that. Tools are located within three subfolders: tools, scripts, and libs.

1. tools are self-contained analysis workflows available from the command line. 
2. scripts are short pieces of code meant to accomplish relatively simple tasks.
3. libs are loosely defined tools or data structures used externally in other workflows or analyses.

Two tools are currently available, bammend and taulysis. Taulysis is a command-line tool for performing survival analysis of polymerase activity. Bammend is a command-line tool for editing basecalls in PacBio BAM files.

Many scripts are available. Examples include a script for rapidly extracting important metrics contained in .bam.pbi files, a script for truncating bam files in time, and a script for slicing and dicing pacbio bams into smaller files.

The libs folder contains a lot of dependencies used in other tools, outside of the biochemistry-toolkit. For instance, QuickKinetics rapidly summarizes PW and IPD information from a subreadset. It's used to quickly rescale kinetics data. BurstMetrics produces metrics related to a particular insertion phenomenon. HtmlScreenShot uses a headless browser to produce static screenshots of Plotly images. The poa lib contains classes, datastructures, and routines for summarizing and visualizing consensus sequences produced from raw data.

The tools are coarsely unit tested using nose (i.e. tested at larger unit scales than they probably should be). Tests should probably be consolidated to a single top-level location. They're currently sprinkled around, near their tested pieces of code.