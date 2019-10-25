rm(list=ls())

library(leaps)
library(dplyr)
library(caret)
library(data.table)
library(plotly)
library(mclust)
library(dbscan)
library(Metrics)
library(randomForest)

path<-"data"

# Load functions
source("code/utility_functions.R")

# Get train and test data
dt<-fread(paste0(path,"/train.csv"))
test<-fread(paste0(path, "/test.csv"))

# Load the structures files
structures<-fread(paste0(path, "/structures.csv"))

# Load molecular geometry features
dt_feat<-fread(paste0(path, "/simple-molecular-geometry-features/train_geom.csv"))
test_feat<-fread(paste0(path, "/simple-molecular-geometry-features/test_geom.csv"))

# Add atomic radii and electronegativity to structures
atomic_radii<-data.table(atom=c('H','C','N','O','F'),
                         radius=c(0.38, 0.77, 0.75, 0.73, 0.71))
en<-data.table(atom=c('H','C','N','O','F'),
               en=c(2.2, 2.55, 3.04, 3.44, 3.98))
fudge_factor<-0.05
atomic_radii[,":="(radius=radius+fudge_factor)]

structures<-Reduce(merge, list(structures, atomic_radii, en))
structures<-structures[order(molecule_name, atom_index)]

# Get number of bonds, bond length means and standard deviations for each each in each molecule
strucs<-structuresTransform(dt1 = structures)

# Merge all train and test both
dt_str<-mergeStructures(dt1=dt, structures = strucs, mol_geom=dt_feat)
test_str<-mergeStructures(dt1=test, structures = strucs, mol_geom=test_feat)

# Create closest
dt_str1<-createClosest(dt1=dt_str)
test_str1<-createClosest(dt1=test_str)

# Merge all
merge_names<-names(dt_str1)[names(dt_str1) %in% names(dt_str)]
dt_str<-merge(dt_str,dt_str1, by=merge_names)
test_str<-merge(test_str, test_str1, by=merge_names)

# Add cos angle features
dt_str<-add_cos_features(dt1=dt_str)
test_str<-add_cos_features(dt1=test_str)

# Create features of min, max, sd. etc
dt_str<-createFeatures(dt=dt_str)
test_str<-createFeatures(dt1=test_str)
dt_str[is.na(dt_str)]<-0
test_str[is.na(test_str)]<-0

# Write to csv
fwrite(dt_str, "data/dt_str.csv")
fwrite(test_str, 'data/test_str.csv')