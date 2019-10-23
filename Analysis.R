library(dplyr)
library(Matrix)
library(MLmetrics)
library(lightgbm)
library(data.table)

Sys.setenv(KMP_DUPLICATE_LIB_OK=TRUE)
Sys.setenv('KMP_DUPLICATE_LIB_OK'=TRUE)

path<-"~/Data/"

# Load functions
source("~/utility_functions.R")

# Load data
dt_str<-fread("dt_str.csv")
test_str<-fread('test_str.csv')

# Add some more variables
transform2<-function(dt){
  
  dt[, molecule_atom_index_0_meaninv1_mean := mean(meaninv1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_1_meaninv1_mean := mean(meaninv1), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_0_meaninv1_sd := sd(meaninv1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_1_meaninv1_sd := sd(meaninv1), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_num_bonds2_mean := mean(num_bonds2), by=list(molecule_name)]
  
  dt[, molecule_atom_index_0_cos_0_mean := mean(cos_0), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_0_sd := sd(cos_0), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_0_min := min(cos_0), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_0_max := max(cos_0), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_cos_0_mean_diff := cos_0 - molecule_atom_index_0_cos_0_mean]
  dt[, molecule_atom_index_0_cos_0_sd_diff := cos_0 - molecule_atom_index_0_cos_0_sd]
  dt[, molecule_atom_index_0_cos_0_min_diff := cos_0 - molecule_atom_index_0_cos_0_min]
  dt[, molecule_atom_index_0_cos_0_max_diff := cos_0 - molecule_atom_index_0_cos_0_max]
  
  dt[, molecule_atom_index_0_cos_0_mean_div := cos_0/molecule_atom_index_0_cos_0_mean]
  dt[, molecule_atom_index_0_cos_0_sd_div := cos_0/molecule_atom_index_0_cos_0_sd]
  dt[, molecule_atom_index_0_cos_0_min_div := cos_0/molecule_atom_index_0_cos_0_min]
  dt[, molecule_atom_index_0_cos_0_max_div := cos_0/molecule_atom_index_0_cos_0_max]
  
  dt[, molecule_atom_index_0_cos_1_mean := mean(cos_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_1_sd := sd(cos_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_1_min := min(cos_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_1_max := max(cos_1), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_cos_1_mean_diff := cos_1 - molecule_atom_index_0_cos_1_mean]
  dt[, molecule_atom_index_0_cos_1_sd_diff := cos_1 - molecule_atom_index_0_cos_1_sd]
  dt[, molecule_atom_index_0_cos_1_min_diff := cos_1 - molecule_atom_index_0_cos_1_min]
  dt[, molecule_atom_index_0_cos_1_max_diff := cos_1 - molecule_atom_index_0_cos_1_max]
  
  dt[, molecule_atom_index_0_cos_1_mean_div := cos_1/molecule_atom_index_0_cos_1_mean]
  dt[, molecule_atom_index_0_cos_1_sd_div := cos_1/molecule_atom_index_0_cos_1_sd]
  dt[, molecule_atom_index_0_cos_1_min_div := cos_1/molecule_atom_index_0_cos_1_min]
  dt[, molecule_atom_index_0_cos_1_max_div := cos_1/molecule_atom_index_0_cos_1_max]
  
  dt[, molecule_atom_index_0_cos_0_1_mean := mean(cos_0_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_0_1_sd := sd(cos_0_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_0_1_min := min(cos_0_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_cos_0_1_max := max(cos_0_1), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_cos_0_1_mean_diff := cos_0_1 - molecule_atom_index_0_cos_0_1_mean]
  dt[, molecule_atom_index_0_cos_0_1_sd_diff := cos_0_1 - molecule_atom_index_0_cos_0_1_sd]
  dt[, molecule_atom_index_0_cos_0_1_min_diff := cos_0_1 - molecule_atom_index_0_cos_0_1_min]
  dt[, molecule_atom_index_0_cos_0_1_max_diff := cos_0_1 - molecule_atom_index_0_cos_0_1_max]
  
  dt[, molecule_atom_index_0_cos_0_1_mean_div := cos_0_1/molecule_atom_index_0_cos_0_1_mean]
  dt[, molecule_atom_index_0_cos_0_1_sd_div := cos_0_1/molecule_atom_index_0_cos_0_1_sd]
  dt[, molecule_atom_index_0_cos_0_1_min_div := cos_0_1/molecule_atom_index_0_cos_0_1_min]
  dt[, molecule_atom_index_0_cos_0_1_max_div := cos_0_1/molecule_atom_index_0_cos_0_1_max]
  
  dt[, molecule_atom_index_1_cos_0_mean := mean(cos_0), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_0_sd := sd(cos_0), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_0_min := min(cos_0), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_0_max := max(cos_0), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_cos_0_mean_diff := cos_0 - molecule_atom_index_1_cos_0_mean]
  dt[, molecule_atom_index_1_cos_0_sd_diff := cos_0 - molecule_atom_index_1_cos_0_sd]
  dt[, molecule_atom_index_1_cos_0_min_diff := cos_0 - molecule_atom_index_1_cos_0_min]
  dt[, molecule_atom_index_1_cos_0_max_diff := cos_0 - molecule_atom_index_1_cos_0_max]
  
  dt[, molecule_atom_index_1_cos_0_mean_div := cos_0/molecule_atom_index_1_cos_0_mean]
  dt[, molecule_atom_index_1_cos_0_sd_div := cos_0/molecule_atom_index_1_cos_0_sd]
  dt[, molecule_atom_index_1_cos_0_min_div := cos_0/molecule_atom_index_1_cos_0_min]
  dt[, molecule_atom_index_1_cos_0_max_div := cos_0/molecule_atom_index_1_cos_0_max]
  
  dt[, molecule_atom_index_1_cos_1_mean := mean(cos_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_1_sd := sd(cos_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_1_min := min(cos_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_1_max := max(cos_1), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_cos_1_mean_diff := cos_1 - molecule_atom_index_1_cos_1_mean]
  dt[, molecule_atom_index_1_cos_1_sd_diff := cos_1 - molecule_atom_index_1_cos_1_sd]
  dt[, molecule_atom_index_1_cos_1_min_diff := cos_1 - molecule_atom_index_1_cos_1_min]
  dt[, molecule_atom_index_1_cos_1_max_diff := cos_1 - molecule_atom_index_1_cos_1_max]
  
  dt[, molecule_atom_index_1_cos_1_mean_div := cos_1/molecule_atom_index_1_cos_1_mean]
  dt[, molecule_atom_index_1_cos_1_sd_div := cos_1/molecule_atom_index_1_cos_1_sd]
  dt[, molecule_atom_index_1_cos_1_min_div := cos_1/molecule_atom_index_1_cos_1_min]
  dt[, molecule_atom_index_1_cos_1_max_div := cos_1/molecule_atom_index_1_cos_1_max]
  
  dt[, molecule_atom_index_1_cos_0_1_mean := mean(cos_0_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_0_1_sd := sd(cos_0_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_0_1_min := min(cos_0_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_cos_0_1_max := max(cos_0_1), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_cos_0_1_mean_diff := cos_0_1 - molecule_atom_index_1_cos_0_1_mean]
  dt[, molecule_atom_index_1_cos_0_1_sd_diff := cos_0_1 - molecule_atom_index_1_cos_0_1_sd]
  dt[, molecule_atom_index_1_cos_0_1_min_diff := cos_0_1 - molecule_atom_index_1_cos_0_1_min]
  dt[, molecule_atom_index_1_cos_0_1_max_diff := cos_0_1 - molecule_atom_index_1_cos_0_1_max]
  
  dt[, molecule_atom_index_1_cos_0_1_mean_div := cos_0_1/molecule_atom_index_1_cos_0_1_mean]
  dt[, molecule_atom_index_1_cos_0_1_sd_div := cos_0_1/molecule_atom_index_1_cos_0_1_sd]
  dt[, molecule_atom_index_1_cos_0_1_min_div := cos_0_1/molecule_atom_index_1_cos_0_1_min]
  dt[, molecule_atom_index_1_cos_0_1_max_div := cos_0_1/molecule_atom_index_1_cos_0_1_max]
  
  dt[, molecule_atom_index_0_abs_cos_0_mean := mean(abs_cos_0), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_0_sd := sd(abs_cos_0), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_0_min := min(abs_cos_0), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_0_max := max(abs_cos_0), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_abs_cos_0_mean_diff := abs_cos_0 - molecule_atom_index_0_abs_cos_0_mean]
  dt[, molecule_atom_index_0_abs_cos_0_sd_diff := abs_cos_0 - molecule_atom_index_0_abs_cos_0_sd]
  dt[, molecule_atom_index_0_abs_cos_0_min_diff := abs_cos_0 - molecule_atom_index_0_abs_cos_0_min]
  dt[, molecule_atom_index_0_abs_cos_0_max_diff := abs_cos_0 - molecule_atom_index_0_abs_cos_0_max]
  
  dt[, molecule_atom_index_0_abs_cos_0_mean_div := abs_cos_0/molecule_atom_index_0_abs_cos_0_mean]
  dt[, molecule_atom_index_0_abs_cos_0_sd_div := abs_cos_0/molecule_atom_index_0_abs_cos_0_sd]
  dt[, molecule_atom_index_0_abs_cos_0_min_div := abs_cos_0/molecule_atom_index_0_abs_cos_0_min]
  dt[, molecule_atom_index_0_abs_cos_0_max_div := abs_cos_0/molecule_atom_index_0_abs_cos_0_max]
  
  dt[, molecule_atom_index_0_abs_cos_1_mean := mean(abs_cos_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_1_sd := sd(abs_cos_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_1_min := min(abs_cos_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_1_max := max(abs_cos_1), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_abs_cos_1_mean_diff := abs_cos_1 - molecule_atom_index_0_abs_cos_1_mean]
  dt[, molecule_atom_index_0_abs_cos_1_sd_diff := abs_cos_1 - molecule_atom_index_0_abs_cos_1_sd]
  dt[, molecule_atom_index_0_abs_cos_1_min_diff := abs_cos_1 - molecule_atom_index_0_abs_cos_1_min]
  dt[, molecule_atom_index_0_abs_cos_1_max_diff := abs_cos_1 - molecule_atom_index_0_abs_cos_1_max]
  
  dt[, molecule_atom_index_0_abs_cos_1_mean_div := abs_cos_1/molecule_atom_index_0_abs_cos_1_mean]
  dt[, molecule_atom_index_0_abs_cos_1_sd_div := abs_cos_1/molecule_atom_index_0_abs_cos_1_sd]
  dt[, molecule_atom_index_0_abs_cos_1_min_div := abs_cos_1/molecule_atom_index_0_abs_cos_1_min]
  dt[, molecule_atom_index_0_abs_cos_1_max_div := abs_cos_1/molecule_atom_index_0_abs_cos_1_max]
  
  dt[, molecule_atom_index_0_abs_cos_0_1_mean := mean(abs_cos_0_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_0_1_sd := sd(abs_cos_0_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_0_1_min := min(abs_cos_0_1), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_abs_cos_0_1_max := max(abs_cos_0_1), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_abs_cos_0_1_mean_diff := abs_cos_0_1 - molecule_atom_index_0_abs_cos_0_1_mean]
  dt[, molecule_atom_index_0_abs_cos_0_1_sd_diff := abs_cos_0_1 - molecule_atom_index_0_abs_cos_0_1_sd]
  dt[, molecule_atom_index_0_abs_cos_0_1_min_diff := abs_cos_0_1 - molecule_atom_index_0_abs_cos_0_1_min]
  dt[, molecule_atom_index_0_abs_cos_0_1_max_diff := abs_cos_0_1 - molecule_atom_index_0_abs_cos_0_1_max]
  
  dt[, molecule_atom_index_0_abs_cos_0_1_mean_div := abs_cos_0_1/molecule_atom_index_0_abs_cos_0_1_mean]
  dt[, molecule_atom_index_0_abs_cos_0_1_sd_div := abs_cos_0_1/molecule_atom_index_0_abs_cos_0_1_sd]
  dt[, molecule_atom_index_0_abs_cos_0_1_min_div := abs_cos_0_1/molecule_atom_index_0_abs_cos_0_1_min]
  dt[, molecule_atom_index_0_abs_cos_0_1_max_div := abs_cos_0_1/molecule_atom_index_0_abs_cos_0_1_max]
  
  dt[, molecule_atom_index_1_abs_cos_0_mean := mean(abs_cos_0), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_0_sd := sd(abs_cos_0), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_0_min := min(abs_cos_0), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_0_max := max(abs_cos_0), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_abs_cos_0_mean_diff := abs_cos_0 - molecule_atom_index_1_abs_cos_0_mean]
  dt[, molecule_atom_index_1_abs_cos_0_sd_diff := abs_cos_0 - molecule_atom_index_1_abs_cos_0_sd]
  dt[, molecule_atom_index_1_abs_cos_0_min_diff := abs_cos_0 - molecule_atom_index_1_abs_cos_0_min]
  dt[, molecule_atom_index_1_abs_cos_0_max_diff := abs_cos_0 - molecule_atom_index_1_abs_cos_0_max]
  
  dt[, molecule_atom_index_1_abs_cos_0_mean_div := abs_cos_0/molecule_atom_index_1_abs_cos_0_mean]
  dt[, molecule_atom_index_1_abs_cos_0_sd_div := abs_cos_0/molecule_atom_index_1_abs_cos_0_sd]
  dt[, molecule_atom_index_1_abs_cos_0_min_div := abs_cos_0/molecule_atom_index_1_abs_cos_0_min]
  dt[, molecule_atom_index_1_abs_cos_0_max_div := abs_cos_0/molecule_atom_index_1_abs_cos_0_max]
  
  dt[, molecule_atom_index_1_abs_cos_1_mean := mean(abs_cos_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_1_sd := sd(abs_cos_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_1_min := min(abs_cos_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_1_max := max(abs_cos_1), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_abs_cos_1_mean_diff := abs_cos_1 - molecule_atom_index_1_abs_cos_1_mean]
  dt[, molecule_atom_index_1_abs_cos_1_sd_diff := abs_cos_1 - molecule_atom_index_1_abs_cos_1_sd]
  dt[, molecule_atom_index_1_abs_cos_1_min_diff := abs_cos_1 - molecule_atom_index_1_abs_cos_1_min]
  dt[, molecule_atom_index_1_abs_cos_1_max_diff := abs_cos_1 - molecule_atom_index_1_abs_cos_1_max]
  
  dt[, molecule_atom_index_1_abs_cos_1_mean_div := abs_cos_1/molecule_atom_index_1_abs_cos_1_mean]
  dt[, molecule_atom_index_1_abs_cos_1_sd_div := abs_cos_1/molecule_atom_index_1_abs_cos_1_sd]
  dt[, molecule_atom_index_1_abs_cos_1_min_div := abs_cos_1/molecule_atom_index_1_abs_cos_1_min]
  dt[, molecule_atom_index_1_abs_cos_1_max_div := abs_cos_1/molecule_atom_index_1_abs_cos_1_max]
  
  dt[, molecule_atom_index_1_abs_cos_0_1_mean := mean(abs_cos_0_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_0_1_sd := sd(abs_cos_0_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_0_1_min := min(abs_cos_0_1), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_abs_cos_0_1_max := max(abs_cos_0_1), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_abs_cos_0_1_mean_diff := abs_cos_0_1 - molecule_atom_index_1_abs_cos_0_1_mean]
  dt[, molecule_atom_index_1_abs_cos_0_1_sd_diff := abs_cos_0_1 - molecule_atom_index_1_abs_cos_0_1_sd]
  dt[, molecule_atom_index_1_abs_cos_0_1_min_diff := abs_cos_0_1 - molecule_atom_index_1_abs_cos_0_1_min]
  dt[, molecule_atom_index_1_abs_cos_0_1_max_diff := abs_cos_0_1 - molecule_atom_index_1_abs_cos_0_1_max]
  
  dt[, molecule_atom_index_1_abs_cos_0_1_mean_div := abs_cos_0_1/molecule_atom_index_1_abs_cos_0_1_mean]
  dt[, molecule_atom_index_1_abs_cos_0_1_sd_div := abs_cos_0_1/molecule_atom_index_1_abs_cos_0_1_sd]
  dt[, molecule_atom_index_1_abs_cos_0_1_min_div := abs_cos_0_1/molecule_atom_index_1_abs_cos_0_1_min]
  dt[, molecule_atom_index_1_abs_cos_0_1_max_div := abs_cos_0_1/molecule_atom_index_1_abs_cos_0_1_max]
  
  return(dt)
}


dt_str[,":="(cos_0_sq=cos_0^2,
             cos_0_inv=1/(2+cos_0),
             cos_1_sq=cos_1^2,
             cos_1_inv=1/(2+cos_1),
             cos_0_1_sq=cos_0_1^2,
             cos_0_1_inv=1/(2+cos_0_1),
             cos_0_flag=ifelse(cos_0 == 1 | cos_0 == -1, 2, 0),
             cos_1_flag=ifelse(cos_1 == 1 | cos_1 == -1, 2, 0),
             cos_0_1_flag = ifelse(cos_0_1 == 1 | cos_0_1 == -1, 2, 0),
             dist_flag=ifelse(abs_distance < 1.067, 2, 0),
             abs_cos_0=abs(cos_0),
             abs_cos_1=abs(cos_1),
             abs_cos_0_1=abs(cos_0_1),
             dist_inv_sq=dist_inv^2,
             dist_inv_cube=dist_inv^3,
             num_atoms_dist_in=num_atoms*dist_inv,
             num_bonds_dist=num_bonds2*abs_distance,
             num_bonds_distinv=num_bonds2*dist_inv,
             meaninv1=1/bond_lengths_mean1,
             meaninv2=1/bond_lengths_mean2)]

dt_str[, ":="(dist_cos_0_flag = cos_0_flag/abs_distance,
              dist_cos_1_flag = cos_1_flag/abs_distance,
              dist_cos_0_1_flag = cos_0_1_flag/abs_distance,
              
              num_bonds_cos_0_flag=num_bonds2 * cos_0_flag,
              num_bonds_cos_1_flag=num_bonds2 * cos_1_flag,
              num_bonds_cos_0_1_flag=num_bonds2 * cos_0_1_flag,
              
              meaninv2_cos_0_flag = cos_0_flag/meaninv2,
              meaninv2_cos_1_flag = cos_1_flag/meaninv2,
              meaninv2_cos_0_1_flag = cos_0_1_flag/meaninv2,
              
              meaninv1_cos_0_flag = cos_0_flag/meaninv1,
              meaninv1_cos_1_flag = cos_1_flag/meaninv1,
              meaninv1_cos_0_1_flag = cos_0_1_flag/meaninv1,
              
              abs_cos_0_sq=abs_cos_0^2,
              abs_cos_1_sq=abs_cos_1^2,
              abs_cos_0_1_sq=abs_cos_0_1^2,
              
              dist_flag_cos_0_flag = dist_flag * cos_0_flag,
              dist_flag_cos_1_flag = dist_flag * cos_1_flag,
              dist_flag_cos_0_1_flag = dist_flag * cos_0_1_flag)]

test_str[,":="(cos_0_sq=cos_0^2,
               cos_0_inv=1/(2+cos_0),
               cos_1_sq=cos_1^2,
               cos_1_inv=1/(2+cos_1),
               cos_0_1_sq=cos_0_1^2,
               cos_0_1_inv=1/(2+cos_0_1),
               cos_0_flag=ifelse(cos_0 == 1 | cos_0 == -1, 2, 0),
               cos_1_flag=ifelse(cos_1 == 1 | cos_1 == -1, 2, 0),
               cos_0_1_flag = ifelse(cos_0_1 == 1 | cos_0_1 == -1, 2, 0),
               dist_flag=ifelse(abs_distance < 1.067, 2, 0),
               abs_cos_0=abs(cos_0),
               abs_cos_1=abs(cos_1),
               abs_cos_0_1=abs(cos_0_1),
               dist_inv_sq=dist_inv^2,
               dist_inv_cube=dist_inv^3,
               num_atoms_dist_in=num_atoms*dist_inv,
               num_bonds_dist=num_bonds2*abs_distance,
               num_bonds_distinv=num_bonds2*dist_inv,
               meaninv1=1/bond_lengths_mean1,
               meaninv2=1/bond_lengths_mean2)]

test_str[, ":="(dist_cos_0_flag = cos_0_flag/abs_distance,
                dist_cos_1_flag = cos_1_flag/abs_distance,
                dist_cos_0_1_flag = cos_0_1_flag/abs_distance,
                
                num_bonds_cos_0_flag=num_bonds2 * cos_0_flag,
                num_bonds_cos_1_flag=num_bonds2 * cos_1_flag,
                num_bonds_cos_0_1_flag=num_bonds2 * cos_0_1_flag,
                
                meaninv2_cos_0_flag = cos_0_flag/meaninv2,
                meaninv2_cos_1_flag = cos_1_flag/meaninv2,
                meaninv2_cos_0_1_flag = cos_0_1_flag/meaninv2,
                
                meaninv1_cos_0_flag = cos_0_flag/meaninv1,
                meaninv1_cos_1_flag = cos_1_flag/meaninv1,
                meaninv1_cos_0_1_flag = cos_0_1_flag/meaninv1,
                
                abs_cos_0_sq=abs_cos_0^2,
                abs_cos_1_sq=abs_cos_1^2,
                abs_cos_0_1_sq=abs_cos_0_1^2,
                
                dist_flag_cos_0_flag = dist_flag * cos_0_flag,
                dist_flag_cos_1_flag = dist_flag * cos_1_flag,
                dist_flag_cos_0_1_flag = dist_flag * cos_0_1_flag)]

dt_str<-transform2(dt_str)
test_str<-transform2(test_str)

if(sum(is.na(dt_str)) > 0) dt_str[is.na(dt_str)]<-0
if(sum(is.na(test_str)) > 0) test_str[is.na(test_str)]<-0

dt_str[mapply(is.infinite, dt_str)]<-0
test_str[mapply(is.infinite, test_str)]<-0

# Load other data
dipole<-fread(paste0(path, "/dipole_moments.csv"))
magnet<-fread(paste0(path, "/magnetic_shielding_tensors.csv"))
mullik<-fread(paste0(path, "/mulliken_charges.csv"))
potent<-fread(paste0(path, "/potential_energy.csv"))
contrs<-fread(paste0(path, "/scalar_coupling_contributions.csv"))

# Create train and valid data
dtsamp<-data.table(dt_str%>%group_by(type)%>%sample_n(40000))

for(col in c('id','molecule_name','atom_index_0','atom_index_1',
             'x_closest_1','y_closest_1','z_closest_1',
             'x_closest_0','y_closest_0','z_closest_0',
             'distance_1','distance_2','atom_index_closest_0',
             'atom_index_closest_1','atom1','atom2','type0')) dtsamp[, (col) := NULL]

samprows<-sample(nrow(dtsamp), 0.6*nrow(dtsamp), replace=F)
holdoutsamp<-dtsamp[-samprows,]
dtsamp<-dtsamp[samprows,]

# lightgbm params
params <- list(boosting_type = 'gbdt',
               objective = "regression",
               metric = 'mae')
for(t in c('2JHH')){
  dtsamp_temp<-dtsamp[type==t]
  valid_temp<-holdoutsamp[type==t]
  
  dtsamp_temp$type<-NULL
  valid_temp$type<-NULL
  
  y_train<-dtsamp_temp$scalar_coupling_constant
  x_train<-dtsamp_temp[, -c('scalar_coupling_constant'), with=F]
  y_valid<-valid_temp$scalar_coupling_constant
  x_valid<-valid_temp[, -c('scalar_coupling_constant'), with=F]
  
  x_train[] <- lapply(x_train, as.numeric)
  x_valid[]<- lapply(x_valid, as.numeric)
  
  dtrain <- lgb.Dataset(as.matrix(x_train), label = y_train)
  dvalid <- lgb.Dataset(as.matrix(x_valid), label = y_valid)
  
  valids <- list(test = dvalid)
  
  clf <- lightgbm(params = params, dtrain, nrounds = 500,
                  valids, early_stopping_rounds=40)
  print(clf$best_score)
  p_test<-predict(clf, as.matrix(x_valid))
  print(MAE(p_test, y_valid))
}


# Create train and valid data
for(col in c('id','molecule_name','atom_index_0','atom_index_1',
             'x_closest_1','y_closest_1','z_closest_1',
             'x_closest_0','y_closest_0','z_closest_0',
             'distance_1','distance_2','atom_index_closest_0',
             'atom_index_closest_1','atom1','atom2','type0')) dt_str[, (col) := NULL]
for(col in c('molecule_name','atom_index_0','atom_index_1',
             'x_closest_1','y_closest_1','z_closest_1',
             'x_closest_0','y_closest_0','z_closest_0',
             'distance_1','distance_2','atom_index_closest_0',
             'atom_index_closest_1','atom1','atom2','type0')) test_str[, (col) := NULL]

# lightgbm params
params <- list(boosting_type = 'gbdt',
               objective = "regression",
               metric='mae')
testpreds<-data.table()

for(t in c('1JHN','2JHC','2JHN','2JHH','3JHC','3JHH','3JHN')){
  
  dtsamp_temp<-dt_str[type==t]
  dtsamp_temp$type<-NULL
  
  y_train<-dtsamp_temp$scalar_coupling_constant
  x_train<-dtsamp_temp[, -c('scalar_coupling_constant'), with=F]
  x_valid<-test_str[type==t][, -c('id'), with=F]
  x_valid$type<- NULL
  
  null_cols<-c()
  for(col in names(x_train)){
    if(sd(x_train[[col]]) == 0) {
      null_cols<-c(null_cols, col)
    }
  }
  
  for(col in null_cols) {
    x_train[, (col) := NULL]
    x_valid[, (col) := NULL]
  }
  
  x_train[] <- lapply(x_train, as.numeric)
  x_valid[] <- lapply(x_valid, as.numeric)
  
  dtrain <- lgb.Dataset(as.matrix(x_train), label = y_train)
  print(t)
  clf <- lightgbm(params = params, dtrain, nrounds = 4000, verbose=-1)
  
  print(MAE(predict(clf, as.matrix(x_train)), y_train))
  
  p_test<-predict(clf, as.matrix(x_valid))
  test_preds<-data.table(id=test_str[type==t]$id,
                         scalar_coupling_constant=p_test)
  
  testpreds<-rbindlist(list(testpreds, test_preds))
}

testpreds[, ":="(id=as.integer(id))]
head(testpreds)

fwrite(testpreds, "submission9.csv")

# Try keras regression
# model <- ...  # create the original model
# 
# layer_name <- 'my_layer'
# intermediate_layer_model <- keras_model(inputs = model$input,
#                                         outputs = get_layer(model, layer_name)$output)
# intermediate_output <- predict(intermediate_layer_model, data)

dtsamp<-data.table(dt_str%>%filter(type=='1JHC'))

for(col in c('id','molecule_name','atom_index_0','atom_index_1',
             'x_closest_1','y_closest_1','z_closest_1',
             'x_closest_0','y_closest_0','z_closest_0',
             'distance_1','distance_2','atom_index_closest_0',
             'atom_index_closest_1','atom1','atom2','type0','type')) dtsamp[, (col) := NULL]

samprows<-sample(nrow(dtsamp), 0.75*nrow(dtsamp), replace=F)
holdoutsamp<-dtsamp[-samprows,]
dtsamp<-dtsamp[samprows,]

train_labels<-dtsamp$scalar_coupling_constant
test_labels<-holdoutsamp$scalar_coupling_constant

dtsamp<-dtsamp[, -c('scalar_coupling_constant'), with=F]
holdoutsamp<-holdoutsamp[, -c('scalar_coupling_constant'), with=F]

null_cols<-c()
for(col in names(dtsamp)){
  if(sd(dtsamp[[col]]) == 0) {
    null_cols<-c(null_cols, col)
  }
}

for(col in null_cols) {
  dtsamp[, (col) := NULL]
  holdoutsamp[, (col) := NULL]
}

# Normalize training data
train_data <- scale(dtsamp) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(train_data, "scaled:center") 
col_stddevs_train <- attr(train_data, "scaled:scale")
test_data <- scale(holdoutsamp, 
                   center = col_means_train, scale = col_stddevs_train)

build_model <- function() {
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 90, activation = "relu",
                input_shape = dim(train_data)[2], name='dense_1',
                kernel_regularizer = regularizer_l2(l = 0.001)) %>%
    layer_dense(units = 256, activation = "relu", name='dense_2',
                kernel_regularizer = regularizer_l2(l = 0.005)) %>%
    layer_dense(units = 512, activation = "relu", name='dense_3',
                kernel_regularizer = regularizer_l2(l = 0.005)) %>%
    layer_dense(units = 128, activation = "relu", name='dense_4',
                kernel_regularizer = regularizer_l2(l = 0.005)) %>%
    layer_dense(units = 256, activation = "relu", name='dense_5',
                kernel_regularizer = regularizer_l2(l = 0.005)) %>%
    layer_dense(units = 1)
  
  model %>% compile(
    loss = "mae",
    optimizer = optimizer_adam(),
    metrics = list("mean_absolute_error")
  )
  
  model
}

model <- build_model()
model %>% summary()

# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)    
# The patience parameter is the amount of epochs to check for improvement.
early_stop <- callback_early_stopping(monitor = "val_loss", patience = 50)

epochs <- 500

# Fit the model and store training stats
history <- model %>% fit(
  train_data,
  train_labels,
  epochs = epochs,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(early_stop, print_dot_callback)
)

plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
  coord_cartesian(xlim=c(0,250), ylim=c(0,5))

tail(history$metrics$val_mean_absolute_error)
tail(history$metrics$mean_absolute_error)

holdoutpreds_nn<-model %>% predict(test_data)
print(MAE(holdoutpreds_nn[,1], test_labels))

# Predict Intermediate layers for LightGBM
for(nround in c(3000,3500,4000,5000,6000)){
  layer_name<-paste0('dense_',1)
  intermediate_layer_model1 <- keras_model(inputs = model$input,
                                           outputs = get_layer(model, 'dense_1')$output)
  intermediate_output1 <- predict(intermediate_layer_model1, train_data)
  intermediate_test1 <- predict(intermediate_layer_model1, test_data)
  
  # Train LGBM
  params <- list(boosting_type = 'gbdt',
                 objective = "regression",
                 metric='mae')
  
  dtrain <- lgb.Dataset(intermediate_output1, label = train_labels)
  
  clf <- lightgbm(params = params, dtrain, nrounds = nround, verbose=-1)
  test_preds<-predict(clf, intermediate_test1)
  print(MAE(test_preds, test_labels))
  print(MAE(predict(clf, intermediate_output1), train_labels))
}


dtrain <- lgb.Dataset(train_data, label = train_labels)

clf <- lightgbm(params = params, dtrain, nrounds = 3000, verbose=-1)
test_preds<-predict(clf, test_data)
print(MAE(test_preds, test_labels))


## 
dtsamp<-data.table(dt_str%>%filter(type=='1JHC'))

for(col in c('id','molecule_name','atom_index_0','atom_index_1',
             'x_closest_1','y_closest_1','z_closest_1',
             'x_closest_0','y_closest_0','z_closest_0',
             'distance_1','distance_2','atom_index_closest_0',
             'atom_index_closest_1','atom1','atom2','type0','type')) dtsamp[, (col) := NULL]

samprows<-sample(nrow(dtsamp), 0.75*nrow(dtsamp), replace=F)
holdoutsamp<-dtsamp[-samprows,]
dtsamp<-dtsamp[samprows,]

train_labels<-dtsamp$scalar_coupling_constant
test_labels<-holdoutsamp$scalar_coupling_constant

dtsamp<-dtsamp[, -c('scalar_coupling_constant'), with=F]
holdoutsamp<-holdoutsamp[, -c('scalar_coupling_constant'), with=F]

null_cols<-c()
for(col in names(dtsamp)){
  if(sd(dtsamp[[col]]) == 0) {
    null_cols<-c(null_cols, col)
  }
}

for(col in null_cols) {
  dtsamp[, (col) := NULL]
  holdoutsamp[, (col) := NULL]
}

# Normalize training data
train_data <- scale(dtsamp) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(train_data, "scaled:center") 
col_stddevs_train <- attr(train_data, "scaled:scale")
test_data <- scale(holdoutsamp, 
                   center = col_means_train, scale = col_stddevs_train)

build_model <- function() {
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 90, activation = "relu",
                input_shape = dim(train_data)[2], name='dense_1',
                kernel_regularizer = regularizer_l2(l = 0.001)) %>%
    layer_dense(units = 256, activation = "relu", name='dense_2',
                kernel_regularizer = regularizer_l2(l = 0.005)) %>%
    layer_dense(units = 1)
  
  model %>% compile(
    loss = "mae",
    optimizer = optimizer_adam(),
    metrics = list("mean_absolute_error")
  )
  
  model
}

model <- build_model()
model %>% summary()

# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)    
# The patience parameter is the amount of epochs to check for improvement.
#early_stop <- callback_early_stopping(monitor = "val_loss", patience = 50)

epochs <- 2000

# Fit the model and store training stats
history <- model %>% fit(
  train_data,
  train_labels,
  epochs = epochs,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(print_dot_callback)
)

plot(history, metrics = "mean_absolute_error", smooth = FALSE) +
  coord_cartesian(xlim=c(0,2000), ylim=c(0,5))

tail(history$metrics$val_mean_absolute_error)
tail(history$metrics$mean_absolute_error)

holdoutpreds_nn<-model %>% predict(test_data)
print(MAE(holdoutpreds_nn[,1], test_labels))

# Try different models of LGBM
type_counts<-dt_str[, .(Count=.N), by=list(type)]

ids<-c()
cum_sum<-0
for(i in type_counts$type){
  num_rows<-type_counts[type==i]$Count
  start=cum_sum
  end=cum_sum+num_rows
  ids<-c(ids, sample(start:end, 0.75*num_rows, replace=F))
  cum_sum<-cum_sum+num_rows
}

dtsamp<-dt_str[id %in% ids]
hosamp<-dt_str[!id %in% ids]

for(col in c('id','molecule_name','atom_index_0','atom_index_1',
             'x_closest_1','y_closest_1','z_closest_1',
             'x_closest_0','y_closest_0','z_closest_0',
             'distance_1','distance_2','atom_index_closest_0',
             'atom_index_closest_1','atom1','atom2','type0')) {
  dtsamp[, (col) := NULL]
  hosamp[, (col) := NULL]
}

params <- list(boosting_type = 'gbdt',
               objective = "regression",
               metric='mae')

for(rounds in c(3000, 4000, 5000, 6000)){
  train_maes<-c()
  test_maes<-c()
  
  importances<-list()
  clf<-list()
  
  for(t in c('1JHC','1JHN','2JHC','2JHN','2JHH','3JHC','3JHH','3JHN')){
    dtsampt<-dtsamp[type==t]
    hosampt<-hosamp[type==t]
    
    dtsampt$type<-NULL
    hosampt$type<-NULL
    
    train_labels<-dtsampt$scalar_coupling_constant
    test_labels<-hosampt$scalar_coupling_constant
    
    dtsampt<-dtsampt[, -c('scalar_coupling_constant'), with=F]
    hosampt<-hosampt[, -c('scalar_coupling_constant'), with=F]
    
    null_cols<-c()
    for(col in names(dtsampt)){
      if(sd(dtsampt[[col]]) == 0) {
        null_cols<-c(null_cols, col)
      }
    }
    
    for(col in null_cols) {
      dtsampt[, (col) := NULL]
      hosampt[, (col) := NULL]
    }
    
    # Normalize training data
    train_data <- scale(dtsampt) 
    
    # Use means and standard deviations from training set to normalize test set
    col_means_train <- attr(train_data, "scaled:center") 
    col_stddevs_train <- attr(train_data, "scaled:scale")
    test_data <- scale(hosampt, 
                       center = col_means_train, scale = col_stddevs_train)
    
    dtrain <- lgb.Dataset(train_data, label = train_labels)
    
    clf[[t]] <- lightgbm(params = params, dtrain, nrounds = rounds, verbose=-1)
    
    importances[[t]]<-lgb.importance(clf[[t]], percentage = T)
    
    temp1<-predict(clf[[t]], train_data)
    temp2<-predict(clf[[t]], test_data)
    
    train_maes<-c(train_maes, MAE(temp1, train_labels))
    test_maes<-c(test_maes, MAE(temp2, test_labels))
    
    print(paste0(t, "-Train: ", MAE(temp1, train_labels)))
    print(paste0(t, "-Test: ", MAE(temp2, test_labels)))
  }
  
  print(0.125*sum(log(train_maes)))
  print(0.125*sum(log(test_maes)))
}

imps1<-importances
clfs1<-clf
trmaes1<-train_maes
tsmaes1<-test_maes


for(rounds in c(7000, 8000, 10000)){
  train_maes<-c()
  test_maes<-c()
  
  importances<-list()
  clf<-list()
  
  for(t in c('1JHC','2JHC','3JHC')){
    dtsampt<-dtsamp[type==t]
    hosampt<-hosamp[type==t]
    
    dtsampt$type<-NULL
    hosampt$type<-NULL
    
    train_labels<-dtsampt$scalar_coupling_constant
    test_labels<-hosampt$scalar_coupling_constant
    
    dtsampt<-dtsampt[, -c('scalar_coupling_constant'), with=F]
    hosampt<-hosampt[, -c('scalar_coupling_constant'), with=F]
    
    null_cols<-c()
    for(col in names(dtsampt)){
      if(sd(dtsampt[[col]]) == 0) {
        null_cols<-c(null_cols, col)
      }
    }
    
    for(col in null_cols) {
      dtsampt[, (col) := NULL]
      hosampt[, (col) := NULL]
    }
    
    # Normalize training data
    train_data <- scale(dtsampt) 
    
    # Use means and standard deviations from training set to normalize test set
    col_means_train <- attr(train_data, "scaled:center") 
    col_stddevs_train <- attr(train_data, "scaled:scale")
    test_data <- scale(hosampt, 
                       center = col_means_train, scale = col_stddevs_train)
    
    dtrain <- lgb.Dataset(train_data, label = train_labels)
    
    clf[[t]] <- lightgbm(params = params, dtrain, nrounds = rounds, verbose=-1)
    
    importances[[t]]<-lgb.importance(clf[[t]], percentage = T)
    
    temp1<-predict(clf[[t]], train_data)
    temp2<-predict(clf[[t]], test_data)
    
    train_maes<-c(train_maes, MAE(temp1, train_labels))
    test_maes<-c(test_maes, MAE(temp2, test_labels))
    
    print(paste0(t, "-Train: ", MAE(temp1, train_labels)))
    print(paste0(t, "-Test: ", MAE(temp2, test_labels)))
  }
  
  print(0.125*sum(log(train_maes)))
  print(0.125*sum(log(test_maes)))
}

# 1JHN
dtsampt<-dtsamp[type=='1JHN']
hosampt<-hosamp[type=='1JHN']

dtsampt$type<-NULL
hosampt$type<-NULL

train_labels<-dtsampt$scalar_coupling_constant
test_labels<-hosampt$scalar_coupling_constant

dtsampt<-dtsampt[, -c('scalar_coupling_constant'), with=F]
hosampt<-hosampt[, -c('scalar_coupling_constant'), with=F]

null_cols<-c()
for(col in names(dtsampt)){
  if(sd(dtsampt[[col]]) == 0) {
    null_cols<-c(null_cols, col)
  }
}

for(col in null_cols) {
  dtsampt[, (col) := NULL]
  hosampt[, (col) := NULL]
}

# Normalize training data
train_data <- scale(dtsampt) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(train_data, "scaled:center") 
col_stddevs_train <- attr(train_data, "scaled:scale")
test_data <- scale(hosampt, 
                   center = col_means_train, scale = col_stddevs_train)

dtrain <- lgb.Dataset(train_data, label = train_labels)


maesdt<-data.table()
for(rounds in c(100, 500, 1000, 3000)){
  for(lr in c(0.0001, 0.001, 0.005, 0.01, 0.05)){
    for(max_depth in c(6,10,15,30)){
      for(max_bin in c(10,50,100,150,250,350)){
        for(num_leaves in c(50, 70, 100, 150, 200)){
          for(min_data_in_leaf in c(10,50,100,120)){
            params <- list(boosting_type = 'gbdt',
                           objective = "regression_l2",
                           metric='mae',
                           # feature_fraction=0.8,
                           # reg_alpha=0.1, 
                           # reg_lambda=0.3,
                           max_bin=max_bin,
                           min_data_in_leaf=min_data_in_leaf,
                           learning_rate=lr,
                           max_depth=max_depth,
                           num_leaves=num_leaves)
            clf1 <- lightgbm(params = params, dtrain, nrounds = rounds, verbose=-1)
            
            temp1<-predict(clf1, train_data)
            temp2<-predict(clf1, test_data)
            
            print(paste0("Train: ", round(MAE(temp1, train_labels), 4)))
            print(paste0("Test: ", round(MAE(temp2, test_labels), 4)))
            
            tempdt<-data.table(rounds=rounds, learning_rate=lr, max_depth=max_depth,
                               max_bin=max_bin, num_leaves=num_leaves, 
                               min_data_in_leaf=min_data_in_leaf, mae=MAE(temp1, train_labels))
            maesdt<-rbindlist(list(maesdt, tempdt))
          }
        }
      }
    }
  }
}
