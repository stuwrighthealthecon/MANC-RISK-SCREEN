# MANC-RISK-SCREEN

This model was developed by Ewan Gray, Anna Donten, and Katherine Payne and validated by Stuart Wright, Gabriel Rogers, Katherine Payne, and Rob Hainsworth at the Manchester Centre for Health Economics. The model is currently maintained by Stuart Wright (stuart.j.wright@manchester.ac.uk)

This discrete event simulation model estimates the QALYs, costs, and clinical outcomes of women attending breast cancer screening delivered using one of 6 strategies or for#no screening. The available screening programmes are 
   1) risk stratified screening using the Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs defined in the PROCAS project.
   2) risk stratified screening using the Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs dividing the population into tertiles
   3) 3 yearly breast cancer screening
   4) 2 yearly breast cancer screening
   5) 5 yearly breast cancer screening
   6) 10 yearly breast cancer screening (at age 50 and 60)
   7) No breast cancer screening
   8) 5 yearly breast cancer screening for women at low risk predict using Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs defined in the PROCAS project
   9) 6 yearly breast cancer screening for women at low risk predict using Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs defined in the PROCAS project
   10) A fully stratified breast cancer screening programme with women at high or moderate risk receiving annual and bi-annual screening respectively and lower risk women receiving 5 yearly screening.

The model was originally published as an early economic evaluation in the journal Value in Health:
Gray E, Donten A, Karssemeijer N, et al. Evaluation of a Stratified National Breast Screening Program in the United Kingdom: An Early Model-Based Cost-Effectiveness Analysis. Value Heal 2017;20:1100–9.

The documentation folder contains information about the model. An overall summary is provided in the Model Summary word document. The accompanying file titled "Model Text Algorithm" was created by Stuart Wright and checked by Ewan Gray in the model validation process. This file contains a step by step text based explanation of the processes taking place in the model. The research team recommend reading this file and the published economic evaluation to gain a full understanding of the model. The documentation folder also contains diagrams outlining the structure of the model, tables of model parameters and their sources, and details about the process of updating the model parameters between version 0 and version 1.

Below are some key notes about the model for new users:
1) The model evalutes one screening strategy at a time (as defined by the user in the screen_strategy parameter). The results are therefore not incremental and need to be compared to those from another strategy which must also be run seperately
2) Like most discrete event simulation models this model is computationally expensive. We have taken steps to decrease run-time and will be looking into this further. 
3) DES models create natural variation in the outcomes and so require a large population to be sampled in the model to achieve stable results. Currently we recommend that users simulate results for at least 7 million women, particularly for the more variable risk-stratified strategies. We are investigating variance reduction techniques to reduce this number.
4) Given the large number of required model runs, the model splits the total sample into 10 sub-samples. If the model is interrupted mid-run then the results from some of these sub-samples may be retrieved to save running the whole sample again
5) Key data in this model include a sample of ~15,000 women's Volpara breast density estimates, estimated 10 year risk of breast cancer, and estimated lifetime risk of breast cancer. This data was based on real data from the PROCAS 2 study. In order to create a model version which could widely be shared while being sensitive to data protection concerns, the researchers created a synthetic dataset of these variables basedon the original data. This used the synthpop package in R which maintains the original structure of the data. 
