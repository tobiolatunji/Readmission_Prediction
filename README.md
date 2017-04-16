# Readmission_Prediction
Prediction algorithm for 10yr, 100,000 patients diabetes dataset

The dataset represents 10 years (1999-2008) of clinical care at 130 US hospitals
and integrated delivery networks. It includes over 50 features representing
patient and hospital outcomes.

Beata Strack, Jonathan P. DeShazo, Chris Gennings, Juan L. Olmo,
Sebastian Ventura, Krzysztof J. Cios, and John N. Clore,
???Impact of HbA1c Measurement on Hospital Readmission Rates: Analysis of 70,000
Clinical Database Patient Records,??? BioMed Research International,
vol. 2014, Article ID 781670, 11 pages, 2014.

<https://archive.ics.uci.edu/ml/datasets/Diabetes+130-US+hospitals+for+years+1999-2008>

The dataset for the readmission prediction was quite huge (100,000  rows) needing a lot of processing power to churn through the machine learning algorithms on my Mac. I had to let some models run overnight just to get some results! I had to stop here and put this up here given the time. The analysis and prediction models could take days to weeks but I think I got some pretty decent results- 

## Preliminary results

up to 93% prediction accuracy for NO readmission in some models. High specificity for <30-day and >30-day readmissions up to 91% in some models.
