## Random then greedy approach for biomarker selection

The script `boosting_markers.R` works through implementation of the Random-then-greedy approach applied to the FS-epsilon boosting algorithm developed by Prof. Rahul Mazumder's group, and assess it's performance across different working data covariance matrices (Inditity, Exchengable, and AR(1)). 

A common problem in biomarker selection for prediction tasks is that the models tend to 
1. be highly correlated 
2. be highly dimensional (e.g. 20k coding genes in the human genome)

Based on the properties of boosting algorithms, we can use the RtG approach to accomodate high correlation across the biomarkers as can be done in greedy approaches to boosting in the FS-epsilon algorithm, but with the speed of of randomized algorithms. 