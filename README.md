### Online paper: ___`Gu C, Wang X, Wang K, et alCryoablation triggers type I interferon-dependent antitumor immunity and potentiates immunotherapy efficacy in lung cancerJournal for ImmunoTherapy of Cancer 2024;12:e008386. [doi: 10.1136/jitc-2023-008386]`___ (https://jitc.bmj.com/content/12/1/e008386)  

### `script`:  
`script` includes all codes of analysis, 
   - .sh: to preprocess fastq data  
   - .r/.py: to process GEX data by mainly Seurat/Scanpy  
   - .ipynb: to use Scanpy&CoNGA to analyze scRNA&scTCR data with output picture but large space size

Under this directory,  
  - `creatObject` includes all scripts from the early stages of the project, which are only used to show and document the generation process of objects;  
   - `finalAnalysis` (mainly under paper/) includes all scripts to visualize all results after creating the objects, which the input data can be generated directly by files from `data` directory.  
    
### `data`:  
`data` includes the data that can be directly used to analyze and re-produce the output of script from `finalAnalysis`.  
Note that some data needs to be reformatted appropriately or merged together manually by user.  
  
### `raw_gex`:  
`raw_gex` includes the unprocessed Gene expression matrix data of scRNA-seq. Use it to process by user's pipeline, or use data from `data` if you wanna re-produce the results in the paper.  
