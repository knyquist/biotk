# `bammend`
`bammend` is a simple, generic pulse-rejection tool for testing pulsecall filters, enabling rapid hypothesis-testing without having to modify the primary analysis stack.  

## How to get started with `bammend`
`bammend` is available as a command-line executable. It should be available within the `primary-toolkit` module via  
```module load primary-toolkit```.  

`bammend` is also available for local installation. Navigate to the root directory and  
```pip install -r requirements.txt```  
```python setup.py install```.

## How to use `bammend`  
`bammend` is simple to use. It has two inputs, a pacbio bam file and a CSV file that marks the indices of basecalls to reject by ZMW. The CSV file has three columns  
| holeNumber | annotationStartIndex | annotationEndIndex |  
|------------|----------------------|--------------------|  

`annotationStartIndex` and `annotationEndIndex` are indexed from the beginning of the ZMW read (i.e. query indices, see http://pacbiofileformats.readthedocs.io/en/5.0/BAM.html). 

Disclaimer

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.