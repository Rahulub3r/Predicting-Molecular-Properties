# Create features
createFeatures<-function(dt){
  dt[, type0 := substr(type, 1, 1)]
  dt[, molecular_couples := length(id), by='molecule_name']
  
  dt[, molecule_dist_mean := mean(abs_distance), by='molecule_name']
  dt[, molecule_dist_min := min(abs_distance), by='molecule_name']
  dt[, molecule_dist_man := max(abs_distance), by='molecule_name']
  
  dt[, atom_0_couples_count := length(id), by=list(molecule_name, atom_index_0)]
  dt[, atom_1_couples_count := length(id), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_0_x_2_std := sd(x2), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_y_2_mean := mean(y2), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_y_2_mean_diff := molecule_atom_index_0_y_2_mean-y2]
  dt[, molecule_atom_index_0_y_2_mean_div := molecule_atom_index_0_y_2_mean/y2]
  dt[, molecule_atom_index_0_y_2_max := max(y2), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_y_2_max_diff := molecule_atom_index_0_y_2_max-y2]
  
  dt[, molecule_atom_index_0_y_2_sd := sd(y2), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_z_2_sd := sd(z2), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_dist_mean := mean(abs_distance), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_mean_diff := molecule_atom_index_0_dist_mean-abs_distance]
  dt[, molecule_atom_index_0_dist_mean_div := molecule_atom_index_0_dist_mean/abs_distance]
  
  dt[, molecule_atom_index_0_dist_max := max(abs_distance), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_max_diff := molecule_atom_index_0_dist_max-abs_distance]
  dt[, molecule_atom_index_0_dist_max_div := molecule_atom_index_0_dist_max/abs_distance]
  
  dt[, molecule_atom_index_0_dist_min := min(abs_distance), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_min_diff := molecule_atom_index_0_dist_min-abs_distance]
  dt[, molecule_atom_index_0_dist_min_div := molecule_atom_index_0_dist_min/abs_distance]
  
  dt[, molecule_atom_index_0_dist_sd := sd(abs_distance), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_sd_diff := molecule_atom_index_0_dist_sd-abs_distance]
  dt[, molecule_atom_index_0_dist_sd_div := molecule_atom_index_0_dist_sd/abs_distance]
  
  dt[, molecule_atom_index_1_dist_mean := mean(abs_distance), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_mean_diff := molecule_atom_index_1_dist_mean-abs_distance]
  dt[, molecule_atom_index_1_dist_mean_div := molecule_atom_index_1_dist_mean/abs_distance]
  
  dt[, molecule_atom_index_1_dist_max := max(abs_distance), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_max_diff := molecule_atom_index_1_dist_max-abs_distance]
  dt[, molecule_atom_index_1_dist_max_div := molecule_atom_index_1_dist_max/abs_distance]
  
  dt[, molecule_atom_index_1_dist_min := min(abs_distance), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_min_diff := molecule_atom_index_1_dist_min-abs_distance]
  dt[, molecule_atom_index_1_dist_min_div := molecule_atom_index_1_dist_min/abs_distance]
  
  dt[, molecule_atom_index_1_dist_sd := sd(abs_distance), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_sd_diff := molecule_atom_index_1_dist_sd-abs_distance]
  dt[, molecule_atom_index_1_dist_sd_div := molecule_atom_index_1_dist_sd/abs_distance]
  
  dt[, molecule_atom_2_dist_mean := mean(abs_distance), by=list(molecule_name, atom2)]
  
  dt[, molecule_atom_2_dist_min := min(abs_distance), by=list(molecule_name, atom2)]
  dt[, molecule_atom_2_dist_min_diff := molecule_atom_2_dist_min-abs_distance]
  dt[, molecule_atom_2_dist_min_div := molecule_atom_2_dist_min/abs_distance]
  
  dt[, molecule_atom_2_dist_sd := sd(abs_distance), by=list(molecule_name, atom2)]
  dt[, molecule_atom_2_dist_sd_diff := molecule_atom_2_dist_sd-abs_distance]
  
  dt[, molecule_type_0_dist_sd := sd(abs_distance), by=list(molecule_name, type0)]
  dt[, molecule_type_0_dist_sd_diff := molecule_type_0_dist_sd - abs_distance]
  
  dt[, molecule_type_dist_mean := mean(abs_distance), by=list(molecule_name, type)]
  dt[, molecule_type_dist_min := min(abs_distance), by=list(molecule_name, type)]
  dt[, molecule_type_dist_max := max(abs_distance), by=list(molecule_name, type)]
  dt[, molecule_type_dist_sd := sd(abs_distance), by=list(molecule_name, type)]
  
  dt[, molecule_type_dist_mean_diff := molecule_type_dist_mean - abs_distance]
  dt[, molecule_type_dist_mean_div := molecule_type_dist_mean/abs_distance]
  dt[, molecule_type_dist_min_diff := molecule_type_dist_min - abs_distance]
  dt[, molecule_type_dist_min_div := molecule_type_dist_min/abs_distance]
  dt[, molecule_type_dist_sd_diff := molecule_type_dist_sd - abs_distance]
  dt[, molecule_type_dist_sd_div := molecule_type_dist_sd/abs_distance]
  dt[, molecule_type_dist_max_diff := molecule_type_dist_max - abs_distance]
  dt[, molecule_type_dist_max_div := molecule_type_dist_max/abs_distance]
  
  dt[, molecule_type_dist_inv_mean := mean(dist_inv), by=list(molecule_name, type)]
  dt[, molecule_type_dist_inv_sd := sd(dist_inv), by=list(molecule_name, type)]
  dt[, molecule_type_dist_inv_min := min(dist_inv), by=list(molecule_name, type)]
  dt[, molecule_type_dist_inv_max := max(dist_inv), by=list(molecule_name, type)]
  
  dt[, molecule_type_dist_inv_mean_diff := dist_inv - molecule_type_dist_inv_mean]
  dt[, molecule_type_dist_inv_sd_diff := dist_inv - molecule_type_dist_inv_sd]
  dt[, molecule_type_dist_inv_min_diff := dist_inv - molecule_type_dist_inv_min]
  dt[, molecule_type_dist_inv_max_diff := dist_inv - molecule_type_dist_inv_max]
  
  dt[, molecule_type_dist_inv_mean_div := dist_inv/molecule_type_dist_inv_mean]
  dt[, molecule_type_dist_inv_sd_div := dist_inv/molecule_type_dist_inv_sd]
  dt[, molecule_type_dist_inv_min_div := dist_inv/molecule_type_dist_inv_min]
  dt[, molecule_type_dist_inv_max_div := dist_inv/molecule_type_dist_inv_max]
  
  dt[, molecule_atom_index_0_dist_inv_mean := mean(dist_inv), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_inv_sd := sd(dist_inv), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_inv_min := min(dist_inv), by=list(molecule_name, atom_index_0)]
  dt[, molecule_atom_index_0_dist_inv_max := max(dist_inv), by=list(molecule_name, atom_index_0)]
  
  dt[, molecule_atom_index_0_dist_inv_mean_diff := dist_inv - molecule_atom_index_0_dist_inv_mean]
  dt[, molecule_atom_index_0_dist_inv_sd_diff := dist_inv - molecule_atom_index_0_dist_inv_sd]
  dt[, molecule_atom_index_0_dist_inv_min_diff := dist_inv - molecule_atom_index_0_dist_inv_min]
  dt[, molecule_atom_index_0_dist_inv_max_diff := dist_inv - molecule_atom_index_0_dist_inv_max]
  
  dt[, molecule_atom_index_0_dist_inv_mean_div := dist_inv/molecule_atom_index_0_dist_inv_mean]
  dt[, molecule_atom_index_0_dist_inv_sd_div := dist_inv/molecule_atom_index_0_dist_inv_sd]
  dt[, molecule_atom_index_0_dist_inv_min_div := dist_inv/molecule_atom_index_0_dist_inv_min]
  dt[, molecule_atom_index_0_dist_inv_max_div := dist_inv/molecule_atom_index_0_dist_inv_max]
  
  dt[, molecule_atom_index_1_dist_inv_mean := mean(dist_inv), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_inv_sd := sd(dist_inv), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_inv_min := min(dist_inv), by=list(molecule_name, atom_index_1)]
  dt[, molecule_atom_index_1_dist_inv_max := max(dist_inv), by=list(molecule_name, atom_index_1)]
  
  dt[, molecule_atom_index_1_dist_inv_mean_diff := dist_inv - molecule_atom_index_1_dist_inv_mean]
  dt[, molecule_atom_index_1_dist_inv_sd_diff := dist_inv - molecule_atom_index_1_dist_inv_sd]
  dt[, molecule_atom_index_1_dist_inv_min_diff := dist_inv - molecule_atom_index_1_dist_inv_min]
  dt[, molecule_atom_index_1_dist_inv_max_diff := dist_inv - molecule_atom_index_1_dist_inv_max]
  
  dt[, molecule_atom_index_1_dist_inv_mean_div := dist_inv/molecule_atom_index_1_dist_inv_mean]
  dt[, molecule_atom_index_1_dist_inv_sd_div := dist_inv/molecule_atom_index_1_dist_inv_sd]
  dt[, molecule_atom_index_1_dist_inv_min_div := dist_inv/molecule_atom_index_1_dist_inv_min]
  dt[, molecule_atom_index_1_dist_inv_max_div := dist_inv/molecule_atom_index_1_dist_inv_max]
  
  dt[,":="(cos_0_sq=cos_0^2,
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
  
  dt[, ":="(dist_cos_0_flag = cos_0_flag/abs_distance,
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
  
  dt[, molecule_type_cos_0_mean := mean(cos_0), by=list(molecule_name, type)]
  dt[, molecule_type_cos_0_sd := sd(cos_0), by=list(molecule_name, type)]
  dt[, molecule_type_cos_0_min := min(cos_0), by=list(molecule_name, type)]
  dt[, molecule_type_cos_0_max := max(cos_0), by=list(molecule_name, type)]
  
  dt[, molecule_type_cos_0_mean_diff := cos_0 - molecule_type_cos_0_mean]
  dt[, molecule_type_cos_0_sd_diff := cos_0 - molecule_type_cos_0_sd]
  dt[, molecule_type_cos_0_min_diff := cos_0 - molecule_type_cos_0_min]
  dt[, molecule_type_cos_0_max_diff := cos_0 - molecule_type_cos_0_max]
  
  dt[, molecule_type_cos_0_mean_div := cos_0/molecule_type_cos_0_mean]
  dt[, molecule_type_cos_0_sd_div := cos_0/molecule_type_cos_0_sd]
  dt[, molecule_type_cos_0_min_div := cos_0/molecule_type_cos_0_min]
  dt[, molecule_type_cos_0_max_div := cos_0/molecule_type_cos_0_max]
  
  dt[, molecule_type_cos_1_mean := mean(cos_1), by=list(molecule_name, type)]
  dt[, molecule_type_cos_1_sd := sd(cos_1), by=list(molecule_name, type)]
  dt[, molecule_type_cos_1_min := min(cos_1), by=list(molecule_name, type)]
  dt[, molecule_type_cos_1_max := max(cos_1), by=list(molecule_name, type)]
  
  dt[, molecule_type_cos_1_mean_diff := cos_1 - molecule_type_cos_1_mean]
  dt[, molecule_type_cos_1_sd_diff := cos_1 - molecule_type_cos_1_sd]
  dt[, molecule_type_cos_1_min_diff := cos_1 - molecule_type_cos_1_min]
  dt[, molecule_type_cos_1_max_diff := cos_1 - molecule_type_cos_1_max]
  
  dt[, molecule_type_cos_1_mean_div := cos_1/molecule_type_cos_1_mean]
  dt[, molecule_type_cos_1_sd_div := cos_1/molecule_type_cos_1_sd]
  dt[, molecule_type_cos_1_min_div := cos_1/molecule_type_cos_1_min]
  dt[, molecule_type_cos_1_max_div := cos_1/molecule_type_cos_1_max]
  
  dt[, molecule_type_cos_0_1_mean := mean(cos_0_1), by=list(molecule_name, type)]
  dt[, molecule_type_cos_0_1_sd := sd(cos_0_1), by=list(molecule_name, type)]
  dt[, molecule_type_cos_0_1_min := min(cos_0_1), by=list(molecule_name, type)]
  dt[, molecule_type_cos_0_1_max := max(cos_0_1), by=list(molecule_name, type)]
  
  dt[, molecule_type_cos_0_1_mean_diff := cos_0_1 - molecule_type_cos_0_1_mean]
  dt[, molecule_type_cos_0_1_sd_diff := cos_0_1 - molecule_type_cos_0_1_sd]
  dt[, molecule_type_cos_0_1_min_diff := cos_0_1 - molecule_type_cos_0_1_min]
  dt[, molecule_type_cos_0_1_max_diff := cos_0_1 - molecule_type_cos_0_1_max]
  
  dt[, molecule_type_cos_0_1_mean_div := cos_0_1/molecule_type_cos_0_1_mean]
  dt[, molecule_type_cos_0_1_sd_div := cos_0_1/molecule_type_cos_0_1_sd]
  dt[, molecule_type_cos_0_1_min_div := cos_0_1/molecule_type_cos_0_1_min]
  dt[, molecule_type_cos_0_1_max_div := cos_0_1/molecule_type_cos_0_1_max]
  
  dt[, molecule_type_abs_cos_0_mean := mean(abs_cos_0), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_0_sd := sd(abs_cos_0), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_0_min := min(abs_cos_0), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_0_max := max(abs_cos_0), by=list(molecule_name, type)]
  
  dt[, molecule_type_abs_cos_0_mean_diff := abs_cos_0 - molecule_type_abs_cos_0_mean]
  dt[, molecule_type_abs_cos_0_sd_diff := abs_cos_0 - molecule_type_abs_cos_0_sd]
  dt[, molecule_type_abs_cos_0_min_diff := abs_cos_0 - molecule_type_abs_cos_0_min]
  dt[, molecule_type_abs_cos_0_max_diff := abs_cos_0 - molecule_type_abs_cos_0_max]
  
  dt[, molecule_type_abs_cos_0_mean_div := abs_cos_0/molecule_type_abs_cos_0_mean]
  dt[, molecule_type_abs_cos_0_sd_div := abs_cos_0/molecule_type_abs_cos_0_sd]
  dt[, molecule_type_abs_cos_0_min_div := abs_cos_0/molecule_type_abs_cos_0_min]
  dt[, molecule_type_abs_cos_0_max_div := abs_cos_0/molecule_type_abs_cos_0_max]
  
  dt[, molecule_type_abs_cos_1_mean := mean(abs_cos_1), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_1_sd := sd(abs_cos_1), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_1_min := min(abs_cos_1), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_1_max := max(abs_cos_1), by=list(molecule_name, type)]
  
  dt[, molecule_type_abs_cos_1_mean_diff := abs_cos_1 - molecule_type_abs_cos_1_mean]
  dt[, molecule_type_abs_cos_1_sd_diff := abs_cos_1 - molecule_type_abs_cos_1_sd]
  dt[, molecule_type_abs_cos_1_min_diff := abs_cos_1 - molecule_type_abs_cos_1_min]
  dt[, molecule_type_abs_cos_1_max_diff := abs_cos_1 - molecule_type_abs_cos_1_max]
  
  dt[, molecule_type_abs_cos_1_mean_div := abs_cos_1/molecule_type_abs_cos_1_mean]
  dt[, molecule_type_abs_cos_1_sd_div := abs_cos_1/molecule_type_abs_cos_1_sd]
  dt[, molecule_type_abs_cos_1_min_div := abs_cos_1/molecule_type_abs_cos_1_min]
  dt[, molecule_type_abs_cos_1_max_div := abs_cos_1/molecule_type_abs_cos_1_max]
  
  dt[, molecule_type_abs_cos_0_1_mean := mean(abs_cos_0_1), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_0_1_sd := sd(abs_cos_0_1), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_0_1_min := min(abs_cos_0_1), by=list(molecule_name, type)]
  dt[, molecule_type_abs_cos_0_1_max := max(abs_cos_0_1), by=list(molecule_name, type)]
  
  dt[, molecule_type_abs_cos_0_1_mean_diff := abs_cos_0_1 - molecule_type_abs_cos_0_1_mean]
  dt[, molecule_type_abs_cos_0_1_sd_diff := abs_cos_0_1 - molecule_type_abs_cos_0_1_sd]
  dt[, molecule_type_abs_cos_0_1_min_diff := abs_cos_0_1 - molecule_type_abs_cos_0_1_min]
  dt[, molecule_type_abs_cos_0_1_max_diff := abs_cos_0_1 - molecule_type_abs_cos_0_1_max]
  
  dt[, molecule_type_abs_cos_0_1_mean_div := abs_cos_0_1/molecule_type_abs_cos_0_1_mean]
  dt[, molecule_type_abs_cos_0_1_sd_div := abs_cos_0_1/molecule_type_abs_cos_0_1_sd]
  dt[, molecule_type_abs_cos_0_1_min_div := abs_cos_0_1/molecule_type_abs_cos_0_1_min]
  dt[, molecule_type_abs_cos_0_1_max_div := abs_cos_0_1/molecule_type_abs_cos_0_1_max]
  
  return(dt)
}

createClosest<-function(dt1){
  dt<-data.table(dt1)
  temp1<-dt[,.(molecule_name, atom_index_0,
               atom_index_1, abs_distance, x1, y1, z1,
               x2, y2, z2)]
  temp2<-temp1
  setnames(temp2, c('atom_index_0','atom_index_1', 'x1', 'y1', 'z1', 'x2', 'y2','z2'),
           c('atom_index_1','atom_index_0','x2','y2','z2','x1','y1','z1'))
  temp1<-dt[,.(molecule_name, atom_index_0,
               atom_index_1, abs_distance, x1, y1, z1,
               x2, y2, z2)]
  
  temp1<-rbindlist(list(temp1, temp2), use.names = T)
  
  temp1[, min_distance := min(abs_distance), by=list(molecule_name, atom_index_0)]
  
  temp1<-temp1[min_distance==abs_distance]
  
  for(col in c('x1','y1','z1','min_distance')) temp1[, (col) := NULL]
  setnames(temp1, c('atom_index_0', 'atom_index_1', 'abs_distance', 'x2','y2','z2'),
           c('atom_index','atom_index_closest',
             'distance_closest','x_closest','y_closest','z_closest'))
  
  temp_train <- data.table(dt[,.(molecule_name, atom_index_0,
                                atom_index_1, abs_distance, x1, y1, z1,
                                x2, y2, z2)])
  
  for (atom_idx in c(0,1)){
    temp_train<-merge(temp_train, temp1, by.x=c('molecule_name',
                                                paste0('atom_index_',atom_idx)),
                      by.y=c('molecule_name','atom_index'))
    cols_to_rename<-c('atom_index_closest','distance_closest',
                      'x_closest','y_closest','z_closest')
    setnames(temp_train, cols_to_rename, 
             paste0(cols_to_rename,"_",atom_idx))
  }
  return(temp_train)
}

add_cos_features<-function(dt1){
  dt<-data.table(dt1)
  dt[,distance_1 := sqrt((x1-x_closest_0)^2+(y1-y_closest_0)^2+(z1-z_closest_0)^2)]
  dt[,distance_2 := sqrt((x2-x_closest_1)^2+(y2-y_closest_1)^2+(z2-z_closest_1)^2)]
  dt[,vec_0_x := (x1-x_closest_0)/distance_1]
  dt[,vec_0_y := (y1-y_closest_0)/distance_1]
  dt[,vec_0_z := (z1-z_closest_0)/distance_1]
  dt[,vec_1_x := (x2-x_closest_1)/distance_2]
  dt[,vec_1_y := (y2-y_closest_1)/distance_2]
  dt[,vec_1_z := (z2-z_closest_1)/distance_2]
  dt[,vec_x := (x2-x1)/abs_distance]
  dt[,vec_y := (y2-y1)/abs_distance]
  dt[,vec_z := (z2-z1)/abs_distance]
  dt[,cos_0_1 := vec_0_x*vec_1_x + vec_0_y*vec_1_y + vec_0_z*vec_1_z]
  dt[,cos_0 := vec_0_x*vec_x + vec_0_y*vec_y + vec_0_z*vec_z]
  dt[,cos_1 := vec_1_x*vec_x + vec_1_y*vec_y + vec_1_z*vec_z]
  for(col in c('vec_0_x','vec_0_y','vec_0_z',
               'vec_1_x','vec_1_y','vec_1_z',
               'vec_x','vec_y','vec_z')) dt[, (col) := NULL]
  return(dt)
}

mergeStructures<-function(dt1, structures, mol_geom){
  dt<-data.table(dt1)
  dt<-merge(dt, structures, by.x=c('molecule_name','atom_index_0'),
            by.y=c('molecule_name', 'atom_index'))
  setnames(dt, c('atom','x','y','z','radius','en',
                 'num_bonds','bond_lengths_mean','bond_lengths_sd'),
           c('atom1', 'x1', 'y1', 'z1','radius1','en1','num_bonds1',
             'bond_lengths_mean1','bond_lengths_sd1'))
  
  dt<-merge(dt, structures, by.x=c('molecule_name','atom_index_1'),
            by.y=c('molecule_name', 'atom_index'))
  setnames(dt, c('atom','x','y','z','radius','en',
                 'num_bonds','bond_lengths_mean','bond_lengths_sd'),
           c('atom2','x2','y2','z2','radius2','en2', 'num_bonds2',
             'bond_lengths_mean2','bond_lengths_sd2'))
  
  dt<-merge(dt, mol_geom[,.(molecule_name, atom_index_0,
                            atom_index_1, num_atoms,
                            flatness_metric, bond_angle_plane,
                            bond_angle_axis)],
            by=c('molecule_name', 'atom_index_0',
                 'atom_index_1'))
  dt[,":="(abs_distance=sqrt(((x1-x2)^2) + ((y1-y2)^2) + ((z1-z2)^2)))]
  dt[,":="(dist_inv=1/abs_distance)]
  
  dt[,":="(xdist=(x1-x2)^2,
           ydist=(y1-y2)^2,
           zdist=(z1-z2)^2)]
  
  dt[is.na(dt)]<-0
  
  for(col in c('num_bonds1', 'en1', 'radius1','bond_lengths_sd1')) dt[, (col) := NULL]
  return(dt)
}

structuresTransform<-function(dt1){
  dt<-data.table(dt1)
  dt$num_bonds<-0
  sdcols=c('x','y','z','radius')
  
  for(col in paste0('BL', 1:30)) dt[, (col) := NA]
  
  for(i in 1:30){
    b<-dt[, shift(.SD, i, type='lead'), .SDcols=sdcols, by='molecule_name']
    names(b)<-c('molecule_name', sdcols)
    
    c<-b[,.(x,y,z)] - dt[,.(x,y,z)]
    c<-c[,.(sqrt(x^2+y^2+z^2))]
    
    dt$dist<-c$V1
    dt$radius_sum<-b$radius+dt$radius
    
    dt[, ":="(num_bonds = ifelse(dist>0.0001 & !is.na(dist) & !is.na(radius_sum) &
                                   dist < radius_sum,
                                 num_bonds+1, num_bonds))]
    
    dt[, paste0('BL',i) := ifelse(dist>0.0001 & !is.na(dist) & !is.na(radius_sum) &
                                    dist < radius_sum,
                                  dist, get(paste0('BL', i)))]
    
  }
  for(i in 1:30){
    b<-dt[, shift(.SD, i, type='lag'), .SDcols=sdcols, by='molecule_name']
    names(b)<-c('molecule_name', sdcols)
    
    c<-b[,.(x,y,z)] - dt[,.(x,y,z)]
    c<-c[,.(sqrt(x^2+y^2+z^2))]
    
    dt$dist<-c$V1
    dt$radius_sum<-b$radius+dt$radius
    
    dt[, ":="(num_bonds = ifelse(dist>0.0001 & !is.na(dist) & !is.na(radius_sum) &
                                   dist < radius_sum,
                                 num_bonds+1, num_bonds))]
    dt[, paste0('BL',i) := ifelse(dist>0.0001 & !is.na(dist) & !is.na(radius_sum) &
                                    dist < radius_sum, dist, get(paste0('BL', i)))]
  }
  
  dt[, bond_lengths_mean := apply(.SD, 1, mean, na.rm=T), .SDcols=paste0('BL',1:30)]
  dt[, bond_lengths_sd := apply(.SD, 1, sd, na.rm=T), .SDcols=paste0('BL', 1:30)]
  
  for(col in paste0('BL',1:30)) dt[, (col) := NULL]
  dt$dist<-NULL
  dt$radius_sum<-NULL
  dt[is.na(bond_lengths_sd)]$bond_lengths_sd<-0
  
  return(dt)
}

# RFs lists
getrfs<-function(dtmain, dtadd, ntree=200, mtry=30){
  rfslist<-list()
  
  # For Mulliken Charges and Magnetic Tensors
  if('atom_index' %in% names(dtadd)){
    reqnames<-names(dtadd)
    reqnames<-reqnames[!reqnames %in% c('molecule_name','atom_index')]
    
    dttemp<-merge(dtmain, dtadd, by.x=c('molecule_name','atom_index_0'),
                  by.y=c('molecule_name','atom_index'))
    setnames(dttemp, reqnames, paste0(reqnames, 1))
    dttemp<-merge(dttemp, dtadd, by.x=c('molecule_name','atom_index_1'),
                  by.y=c('molecule_name','atom_index'))
    setnames(dttemp, reqnames, paste0(reqnames, 2))
    
    targets<-c(paste0(reqnames, 1), paste0(reqnames, 2))
  }
  # For Contributions
  else if ('atom_index_0' %in% names(dtadd)) {
    reqnames<-names(dtadd)
    reqnames<-reqnames[!reqnames %in% c('molecule_name','atom_index_0','atom_index_1',
                                        'type')]
    
    dttemp<-merge(dtmain, dtadd, by=c('molecule_name','atom_index_0','atom_index_1','type'))
    
    targets<-reqnames
  }
  
  # For Dipole Moments and Potential Energy
  else {
    reqnames<-names(dtadd)
    reqnames<-reqnames[!reqnames %in% c('molecule_name','atom_index')]
    
    dttemp<-merge(dtmain, dtadd, by=c('molecule_name'))
    
    targets<-reqnames
  }
  for(col in c('id','molecule_name','atom_index_0','atom_index_1',
               'x_closest_1','y_closest_1','z_closest_1',
               'x_closest_0','y_closest_0','z_closest_0',
               'distance_1','distance_2','atom_index_closest_0',
               'atom_index_closest_1','atom1','atom2','type0')) dttemp[, (col) := NULL]
  
  dttemp[,":="(distance_closest_0_inv=1/(distance_closest_0+2),
               distance_closest_1_inv=1/(distance_closest_1+2),
               molecule_atom_index_1_dist_min_inv=1/(molecule_atom_index_1_dist_min+2),
               molecule_atom_index_0_dist_min_inv=1/(molecule_atom_index_0_dist_min+2),
               bond_lengths_mean1_inv=1/(bond_lengths_mean1+2),
               molecule_atom_index_0_dist_max_diff_inv=1/(molecule_atom_index_0_dist_max_diff+2),
               num_bonds2_inv=1/num_bonds2,
               num_bonds2_sq=num_bonds2^2)]
  
  dttemp$type<-NULL
  
  for(target in targets){
    vars1<-names(dttemp)
    vars1<-vars1[!vars1 %in% c('scalar_coupling_constant',targets)]
    vars<-paste(vars1, collapse = "+")
    form <- as.formula(paste(target, '~', vars))
    #print(form)
    rfmod<-randomForest(form, method='anova',
                        data=dttemp,
                        ntree=ntree, mtry=mtry, nodesize=10)
    
    rfslist[[target]]<-rfmod
  }
  return(rfslist)
}