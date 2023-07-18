# proteomics_analysis

Back up your proteomics data every time you run an experiment!

## Data management
### 1. Create a project folder for each experiment you run
The folder should be named in three parts:        
1. A project name - Name your project in a meaningful way
2. Put your initials
3. Put the month and year you run this experiment

For example, in May 2022, Alice Smith wanted to run a proteomics experiment in antibiotic resistant of KT2440.      
First, she should create a personal folder in 'C:/Xcalibur/data' as 'AliceSmith'.      
Second, in her personal folder, she should create a project folder for this experiment called 'ABRKT2440_May2022'. 
### 2. Create a metadata table
This should be the most important document that make everyone including yourself to understand the experimental design. Each project should have its own metadata table. It should be informative, and precisely labeled. 
It can be saved into any of the file formats of your choice. 

For example, Alice had an experiment with three replicates of each control and treatment group, she should have a metadata table like this:
|Date|Treatment | Replicate | Organizm| Data file name|Note|    
|--- | --- | --- | ---|--- | ---|        
|20220501| control | 1 | psudomonas| KT2440_ctrol1_0501|---|      
|20220503| control | 2 | psudomonas| KT2440_ctrol2_0503|---|      
|20220504| control | 3 | psudomonas| KT2440_ctrol3_0504|---|      
|20220505| kanamycin | 1 | psudomonas| KT2440_kan1_0505|---|      
|20220507| kanamycin | 2 | psudomonas| KT2440_kan2_0507|---|      
|20220509| kanamycin | 3 | psudomonas| KT2440_kan3_0509|---|      

If you need to repeat the run of some samples, please make notes, like:   
|Date|Treatment | Replicate | Organizm| Data file name|Note|    
|--- | --- | --- | ---|--- | ---|        
|202205010| control | 4 | psudomonas| KT2440_ctrol4_05010|repeat control1|      
|202205013| control | 5 | psudomonas| KT2440_ctrol5_05013|repeat control2|      

This should be created and stored in the project folder before you start running any of your samples. 
### 3. When running samples, adjust your metadata
Keep track of every run you did and adjust your metadata table as needed. 

### 4. When a run is finished
Check your MS peaks to see if they look good and decide if you need to repeat the sample or not. Check the file size of each run to make sure they include sufficient information. If more samples are needed, adjust the metadata table. 

### 5. Transfer data for analysis
We now transfer the data to the HPRC for database search, QC, and analysis. You can also mark in the metadata, whether or not the data has been processed.

### 6. When finishing an experiment
After database search and basic QC on the results, you should decide whether more samples needed. If yes, add samples in the metadata table and run repeats. If no, follow the below process to back up data on the lab hard drive. **Double-check your metadata table matches all the data files in an experiment**

## Data storage
This step should only happen when you are satisfied with all your data in an experiment. Transfer the data from local computer to the hard drive. 
* If this is the first time you back up your data on the hard drive, transfer the entire personal folder to the hard drive.
* If you have backed up your data before, only transfer the new project folder to the hard drive.

**Important note! To make sure the drive is used efficiently**
  * Only back up raw data and metadata table for each experiment. Delete all the intermediate files generated from the analysis. Delete proteome databases. You can keep one method file for each experiment.
  * Once the data is backed up, delete the project folder from the local computer. We don't keep the same data in both local and hard drives, to avoid confusion. 
