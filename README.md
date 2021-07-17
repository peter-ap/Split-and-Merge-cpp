# SaM RANSAC-Lines Corners

author@ Aerts Peter
e-mail@ peter.aerts@kuleuven.be

2D LiDAR data is clustered, lines and corners are calculated

->2D LiDAR converted from polar to carthesian
->Clustering based on data jumps (rough clustering)     
->Split-and-Merge to further split the clusters to lines
->Lines calulated based on RANSAC and regression line. 
->Corners between consecutive lines are calculated
->Virtual corners are discarded

-> Visualization using Python within C++ (matplotlibcpp.h file)


One downside -> the JumpSize for rough clustering is a fixed number (0.35m)
             -> the beakValue for SaM is also a fixed number ()