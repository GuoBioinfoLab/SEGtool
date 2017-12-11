#



SEGtool <- function (X, exp_cutoff = 3,multi_cpu = 4,detect_mod=2,result_outdir='SEGtool_result',draw_heatmap=TRUE,draw_pca=FALSE,draw_plot=FALSE,html_report=FALSE){
	.makedir(result_outdir)
	filter_bool_list <- apply(X, 1, .filter,exp_cutoff=exp_cutoff,detect_mod=detect_mod)
	Candidate_TS_G <- X[filter_bool_list,]
	Candidate_TS_names <- rownames(Candidate_TS_G)
	Uniform_list <- X[!filter_bool_list,]
	tmp_list <- split(Candidate_TS_G,rownames(Candidate_TS_G))
	if (grepl('linux',R.version$os)){
		library(parallel)
		clusterd_id_list <- mclapply(tmp_list,.FCM_cluster_id,mc.cores=multi_cpu,detect_mod=detect_mod)
		all_matrix <- do.call(rbind,mclapply(Candidate_TS_names,function(x){.judge_TS(Candidate_TS_G[x,],clusterd_id_list[x][[1]][[1]],clusterd_id_list[x][[1]][[2]],detect_mod,clusterd_id_list[x][[1]][[3]])},mc.cores=multi_cpu))
	}else if (grepl('mingw',R.version$os)){
		clusterd_id_list <- lapply(tmp_list,.FCM_cluster_id,detect_mod=detect_mod)
		all_matrix <- do.call(rbind,lapply(Candidate_TS_names,function(x){.judge_TS(Candidate_TS_G[x,],clusterd_id_list[x][[1]][[1]],clusterd_id_list[x][[1]][[2]],detect_mod,clusterd_id_list[x][[1]][[3]])}))
	}else{
		stop('This package could only run at linux or windows os!')
	}
#	matrix_cts_names <- matrix(Candidate_TS_names)
#	candidate_TS_matrix <- do.call(rbind,mclapply(Candidate_TS_names,function(x){.judge_TS(Candidate_TS_G[x,],clusterd_id_list[x][[1]][[1]],clusterd_id_list[x][[1]][[2]],detect_mod,clusterd_id_list[x][[1]][[3]])},mc.cores=multi_cpu))
#	candidate_TS_matrix <- apply(matrix_cts_names,1,function(x){.judge_TS(Candidate_TS_G[x,],clusterd_id_list[x][[1]][[1]],clusterd_id_list[x][[1]][[2]],detect_mod,clusterd_id_list[x][[1]][[3]])})
#	summary_TS <- .get_TS_result(candidate_TS_matrix,Candidate_TS_G,result_outdir,draw_plot=draw_plot,draw_heatmap=draw_heatmap,draw_pca=draw_pca)
#        all_matrix<-t(candidate_TS_matrix)
	candidate_TS_matrix <- t(all_matrix)
        all_df <- as.data.frame(all_matrix,row.names=rownames(Candidate_TS_G))
        colnames(all_df)<-colnames(Candidate_TS_G)
        tissue_df <- as.data.frame(candidate_TS_matrix,row.names=colnames(Candidate_TS_G))
        colnames(tissue_df) <- rownames(Candidate_TS_G)
        candicate_tissue_df <- tissue_df[apply(tissue_df,1,function(z){return(any(z!=0))}),]
        TS_df <- all_df[apply(all_df,1,function(z){return(any(z!=0))}),]
        TS_df_names <- names(TS_df)
        TS_gene_names <- rownames(TS_df)
	TS_exp <- Candidate_TS_G[TS_gene_names,]
        if (grepl('linux',R.version$os)){
                TS_pvalue <- mclapply(TS_gene_names,function(x){return(.get_p_value(expression=TS_exp[x,],TS_vector=TS_df[x,]))},mc.cores=multi_cpu)
	}else if (grepl('mingw',R.version$os)){
		TS_pvalue <- lapply(TS_gene_names,function(x){return(.get_p_value(expression=TS_exp[x,],TS_vector=TS_df[x,]))})
	}	
	TS_P_V <- data.frame(do.call(rbind,TS_pvalue),stringsAsFactors=F)
	TS_P_V$max_exp<-as.numeric(TS_P_V$max_exp)
	TS_P_V$max_TS_p_value<-as.numeric(TS_P_V$max_TS_p_value)
	TS_P_V$p_len<-as.numeric(TS_P_V$p_len)
	rownames(TS_P_V)<-TS_gene_names
	order_TS_P_V <- TS_P_V[order(-TS_P_V$max_TS_p_value,TS_P_V$p_len,-TS_P_V$max_exp),]
	p_value_file <- paste(result_outdir,"SEGs.pvalue",sep="/")
	write.table(order_TS_P_V,file=p_value_file,col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)
	summary_all_samples<- dim(X)
	order_TS_df <- TS_df[rownames(order_TS_P_V),]
	summary_TS <- .get_TS_result(candidate_TS_matrix,Candidate_TS_G,result_outdir,draw_plot=draw_plot,draw_heatmap=draw_heatmap,draw_pca=draw_pca,order_name=rownames(order_TS_P_V))
	summary_result <- c(summary_all_samples,summary_TS)
	summary_result <-data.frame(matrix(summary_result,1))
	names(summary_result) <- c("total_genes","total_samples","total_SEGs","samples_have_SEG","high_SEG","low_SEG","overlap_high_low_SEG")
	summary_result_file <- paste(result_outdir,"summary_result.txt",sep="/")
	write.table(summary_result,file=summary_result_file,col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)
	if (html_report){
		.get_visualisation_html(indir=result_outdir,page.name="SEGtool_result", page.title="Specifically expressed gene analysis results", draw_heatmap=draw_heatmap, draw_pca=draw_pca)
	}
	return(list(Allsummary=summary_result,SEGinSample=order_TS_df,p_value=order_TS_P_V))
}

.makedir <- function(target_dir = getwd()){
	if(file.exists(target_dir)){
		if(file.info(target_dir)$isdir){
	                print(paste(target_dir," is not empty!!!"))
		}else{
			dir.create(target_dir, recursive=TRUE)
		}
        }else{
                dir.create(target_dir, recursive=TRUE)
                message(sprintf("'%s' has been created.", target_dir))
        }
}

.tukeybw <- function (x, c = 5, epsilon = 1e-04,mod=1){
	x <- as.numeric(x)
        m <- median(x)
        s <- median(abs(x - m))
        u <- (x - m)/(c * s + epsilon)
        w <- rep(0, length(x))
        i <- abs(u) <= 1
        w[i] <- ((1 - u^2)^2)[i]
        t.bi <- sum(w * x)/sum(w)
	if (mod==1){
		t.bi <- ifelse (t.bi>=0.5,t.bi,0.5)
	}
        return(t.bi)
}

########filter the low exp vector which don't step into the follow analysis###########
.filter<-function(x,exp_cutoff,detect_mod){
	na_sum <- sum(is.na(x))
	if (na_sum / length(x) > 0.2){
		return(FALSE)
	}
	if (any(x>exp_cutoff)){
		max_x <- max(x)
		min_x <- min(x)
		flat_cutoff <- .get_fold_constant(max_x,detect_mod)
		tbw_value <- .tukeybw(x)
		if (min_x==0){
			min_x <- 1e-04
		}
		flat_flag <- ifelse(max_x/tbw_value>=flat_cutoff | tbw_value/min_x>=flat_cutoff,TRUE,FALSE)
	}else{
		flat_flag <- FALSE
	}
	return(flat_flag)
}

.percent_judge_4_HTS <- function(expression,high_group,detect_mod){
	exp_length <- length(expression)
	high_len <- length(high_group)
	high_index_ratio <- high_len / exp_length
	if ( high_index_ratio <= 0.25 ){
		high_exp_ratio_cutoff <- (0.75 - (0.25-high_index_ratio))*((detect_mod/3)*(1+high_index_ratio))
		non_high_group<-setdiff(expression,high_group)
		lower_high_exp <- max(non_high_group)
		m_h_g <- min(high_group)
		max_non_high_group_cutoff <- ifelse(lower_high_exp<=5,5,5*lower_high_exp)
		if (all(non_high_group<=1) & m_h_g>=max_non_high_group_cutoff){return(high_group)}
		lower_high_exp <- max(non_high_group)
		expression[expression<=1] <- 2^expression[expression<=1]
		total_exp <- sum(expression)
		high_exp_ratio <- sum(high_group)/total_exp
#		print(c(high_group,high_exp_ratio,high_exp_ratio_cutoff))
		if (high_exp_ratio >=high_exp_ratio_cutoff){
			m_h_g_cutoff <- m_h_g/(.get_fold_constant(m_h_g,detect_mod)*(detect_mod/3))*(1+detect_mod/15)
#			print(c(lower_high_exp,m_h_g_cutoff,m_h_g))
#			print(high_group)
			if (lower_high_exp<=m_h_g_cutoff){
#				print('ok')
				return(high_group)
			}else{
				if (length(high_group)>1){
					return(.percent_judge_4_HTS(expression,high_group[-which.min(high_group)],detect_mod))
				}else{
					return(FALSE)
				}
			}
		}else{
			return(FALSE)
		}
	}else{
		return(FALSE)
	}
}

.tbw_std <- function(x){
	x.std <- sqrt(sum((x - .tukeybw(x))^2))
	return(x.std)
}

.get_p_value <- function(expression,TS_vector){
	max_exp <- max(expression)
	tbw_value <- .tukeybw(expression,mod=2)
	sum_exp <- as.integer(sum(expression))
	new_TS_vector <- as.vector(TS_vector)
	new_TS_vector[new_TS_vector>0] <- "greater"
	new_TS_vector[new_TS_vector<0] <- "less"
	TS_index <- which(new_TS_vector!=0)
	p_len <- length(TS_index)
	TS_p_value <- unlist(lapply(TS_index,function(x){single_p_value <- binom.test(as.integer(expression[x]),sum_exp,tbw_value/sum_exp,alternative=as.character(new_TS_vector[x]))$p.value;single_p_value <- -log10(single_p_value);return(single_p_value)}))
#	TS_p_value <- paste(TS_p_value,collapse=";")
	max_TS_p_value <- max(TS_p_value)
	TS_p_value <- paste(TS_p_value,collapse=";")
	return(cbind(max_exp,TS_p_value,max_TS_p_value,p_len))
}

.percent_judge_4_LTS <- function(expression,exp_vector,tbw_value,detect_mod,low_flag,old_l_group){
	TS_l_flag <- FALSE
	.judge_low_group<-function(e_vector){
		more_than_1_e_vector <- e_vector[which(e_vector>1)]
		less_than_1_e_vector <- e_vector[which(e_vector<=1)]
		if (length(more_than_1_e_vector)>=1 & length(less_than_1_e_vector)>=1){
			min_m_1_element <- min(more_than_1_e_vector)
			max_l_1_element <- max(less_than_1_e_vector)
			back_judge_threshold <- min_m_1_element/.get_fold_constant(min_m_1_element,detect_mod)
			if (max_l_1_element<=back_judge_threshold){
				return(less_than_1_e_vector)
			}else{
				return(FALSE)
			}
		}else{
			return(FALSE)
		}
	}
#	cand_low <- exp_vector[,-which.max(exp_vector)]
	low_ratio_cutff <- 0.2 - 0.05*detect_mod
	low_exp_cutoff <- 0.015 - 0.0045*detect_mod	
        tot_length <- length(expression)
        low_length <- length(old_l_group)
        low_ratio <- round(low_length/tot_length,detect_mod)
	min_oldlg<-min(old_l_group)
	min_oldlg <- ifelse(min_oldlg<=1,2^min_oldlg,min_oldlg)
	tot_exp <- sum(expression)
#	print(c(low_ratio_cutff,low_ratio,low_exp_cutoff,tot_exp*low_exp_cutoff))
#	if (min_oldlg>tot_exp*low_exp_cutoff | low_ratio > low_ratio_cutff ){
#		return(FALSE)
#	}
	low_cutoff <- tbw_value/.get_fold_constant(tbw_value,detect_mod)
	l_group <- old_l_group
	exp_l <- exp_vector
#	print(c(low_cutoff,))
#	print(l_group)
#	print(exp_l)
	if (any(old_l_group<1) | any(exp_vector <1)){
		l_group <- old_l_group
		l_group[,l_group<1] <- 2^l_group[,l_group<1]
		exp_l[,exp_l<1] <- 2^exp_l[,exp_l<1]
	}
#	print(c(low_cutoff,l_group[l_group>= low_cutoff]))
	if (any(l_group>= low_cutoff)){
		l_tbw_value <- .tukeybw(exp_l)
		l_tbw_std <- .tbw_std(exp_l)
#		print(l_tbw_value)
#		print(max(l_group)/.get_fold_constant(max(l_group),3))
		low_index_logic<- l_group <= max(l_group)/.get_fold_constant(max(l_group),detect_mod)
		if (length(which(low_index_logic))==0){
			low_result <- c(0)
		}else{
			low_result <- l_group[,which(low_index_logic)]
		}
		if (length(l_group)==1){
			min_more_than_low<-l_group
		}else{
			min_more_than_low <- min(l_group[,l_group > max(l_group)/.get_fold_constant(max(l_group),detect_mod)])
		}
		if (max(low_result)<=min_more_than_low/.get_fold_constant(min_more_than_low,detect_mod)){
			TS_l_flag <- low_result
		}
	}else{
		if (sum(l_group)<=tot_exp*low_exp_cutoff){
			TS_l_flag <- old_l_group
		}
	}
	if (any(TS_l_flag==FALSE)){
		TS_l_flag <- .judge_low_group(exp_vector)
	}	
	return(TS_l_flag)
}

##################################################
.FCM_cluster_id<- function(data_matrix, group_number=3, beta_index = 2, iter_max = 100, cutoff = 0.7, epsilon = 1e-06,detect_mod=2) {
	.Generate_centers <- function(data_matrix,group_number,tbw_value,iter_n_time){
		sd_tbw <- sqrt(sum((data_matrix-tbw_value)^2)/(length(data_matrix)-1))
		center_var <- iter_n_time*sd_tbw*0.5
		center_up_value <- tbw_value+center_var
		center_down_value <- tbw_value-center_var
		center_down_value <- ifelse ( center_down_value <0.001,0,center_down_value)
		out_matrix <- as.matrix(c(center_down_value,tbw_value,center_up_value))
		return(out_matrix)
	}
	.FCM_1 <- function(data_matrix, init_center, beta_index = 2, iter_max = 100, epsilon = 1e-06) {
		.judge_cluster_number <- function(Cluster_id,data_matrix){
			new_id <- rep(0,length(Cluster_id))
			max_value_id <- Cluster_id[which.max(data_matrix)]
			new_id[Cluster_id==max_value_id] <- 3
			min_value_id <- Cluster_id[which.min(data_matrix)]
			new_id[Cluster_id==min_value_id] <- 1
			new_id[new_id==0]<- 2
			return(new_id)
		}
		center_number <- dim(init_center)[1]
		histJ <- c()
		iter_flag <- iter_count <- 1
		dist_PM_1 <- Inf
		dist_matrix <- matrix(0, n_row, center_number)
		for (i in 1:center_number) {
			dist_matrix[, i] = rowSums(sweep(data_matrix, 2, init_center[i, ], "-")^2)
		}
		while(iter_flag){
			para_x <- (1/(dist_matrix + epsilon))^(1/(beta_index - 1))
			para_u <- para_x/(para_x %*% matrix(1, center_number, center_number))
			t_u1 <- t(para_u^beta_index) %*% data_matrix
			t_u2 <- t(para_u^beta_index) %*% matrix(1, n_row, n_col)
			set_coefficients <- t_u1 / t_u2
			new_dist_matrix <- matrix(0, n_row,center_number)
			for (i in 1:center_number) {
				new_dist_matrix[, i] = rowSums(sweep(data_matrix, 2, init_center[i, ], "-")^2)
			}
			dist_PM <- sum(para_u^beta_index * new_dist_matrix)
			iter_flag <- abs(dist_PM - dist_PM_1) > 1e-04 && (iter_count < iter_max)
			dist_PM_1 <- dist_PM
			histJ <- c(histJ, iter_flag)
			iter_count <- iter_count + 1
		}
		Cluster_id <- apply(para_u, 1, which.max)
		Cluster_id <- .judge_cluster_number(Cluster_id,data_matrix)
		Cluster_result <- list(para_u, dist_PM, histJ, set_coefficients, Cluster_id)
		names(Cluster_result) <- c("para_u", "dist_PM_1", "histJ", "init_center", "Cluster_id")
		return(Cluster_result)
	}
	.judge_cluster_id <- function(gene_exp_vector){
		group_FCM_time <- length(gene_exp_vector)
		flag_freq <- table(gene_exp_vector)/ group_FCM_time
		order_index <- order(flag_freq,decreasing=T)
		High_Freq_index <- as.integer(names(flag_freq[order_index[1]]))
		High_Freq <- as.numeric(flag_freq[order_index[1]])
		sample_flag <- ifelse(cutoff <= High_Freq,High_Freq_index,0)
		if (!sample_flag ){
			tmp_flag <- sort(c(names(flag_freq[order_index[1]]),names(flag_freq[order_index[2]])),decreasing=T)
			sample_flag <- as.integer(paste(tmp_flag[1],tmp_flag[2],sep=""))
		}
		return(sample_flag)
	}
	.Group_analysis <- function(group_matrix,cutoff = cutoff){
		new_matrix <- apply(group_matrix,1,.judge_cluster_id)
		return(new_matrix)
	}
	data_matrix <- as.matrix(as.numeric(data_matrix))
        detect_direction <- 2
        tbw_value <- .tukeybw(data_matrix)
	min_x <- min(data_matrix)
	detect_tbw_cutoff<-tbw_value/(.get_fold_constant(tbw_value,detect_mod)*(0.6+detect_mod/10))
        if (min_x>detect_tbw_cutoff){
                detect_direction <- 1
        }
	n_row <- nrow(data_matrix)
	n_col <- ncol(data_matrix)
	potent_iter_time <- round(n_row*0.4)
	n_begin <- ifelse(potent_iter_time>=10,10,potent_iter_time)
	data_matrix_center <- .Generate_centers(data_matrix,group_number,tbw_value,1)
	Find_Cluster_id <- .FCM_1(data_matrix,data_matrix_center,beta_index = beta_index,iter_max = iter_max,epsilon = epsilon)
	if (n_begin > 1){
		n_time = 2
		group_matrix <- matrix(0,n_row,n_begin)
		group_matrix[,1] <- Find_Cluster_id$Cluster_id
		while(n_time <= n_begin){
			new_data_matrix_center <- .Generate_centers(data_matrix,group_number,tbw_value,n_time)
			new_result <- .FCM_1(data_matrix,new_data_matrix_center,beta_index = beta_index,iter_max = iter_max,epsilon = epsilon)
			group_matrix[,n_time] <- new_result$Cluster_id
			if (new_result$dist_PM_1 <= Find_Cluster_id$dist_PM_1 ) {
				Find_Cluster_id <- new_result
			}
			n_time <- n_time + 1
		}
		Find_Cluster_id$Cluster_id <- .Group_analysis(group_matrix,cutoff = cutoff)
	}
	return(list(Find_Cluster_id$Cluster_id,detect_direction,tbw_value))
}

.get_pow <- function(y,z){
        pow_result <- 10^(log(y,10)/z)
        return(pow_result)
}

.get_fold_constant <- function(x,detect_mod){
	ten_index <- log(x,10)
	fold_constant <- 3/ten_index*(1.8*(1+(ten_index-1)/10))
	if (detect_mod==1){
		fold_constant <- exp(1)/ten_index*((exp(1)-1)*(1+(ten_index-1)/10))
	}else if(detect_mod==3){
		fold_constant <- exp(1)/ten_index*(exp(1)*(1+(ten_index-1)/10))
	}
	return(fold_constant)
}

######find the potential TS group judge #####
.get_Th_id <- function(order_index,exp_vector,index_matrix,out_flag,tbw_value,detect_mod){
#####judge if there is a TS group#####
####NO group####
#	print(order_index)
#	print(exp_vector)
#	print(index_matrix)
#	print('in')
	if( length(index_matrix) == 2){
		return(order_index[1:index_matrix[,2][1]])
#####more than 2 meas it maybe exist a group######
	}else{
########The first one is HTS then find the buddies#####
		if( index_matrix[,2][1]==1 ){
			potent_index <- order_index[2:index_matrix[,2][2]]
			potent_group <- exp_vector[potent_index]
			if ((max(potent_group) / min(potent_group) <1.8) | .get_pow(max(potent_group) / min(potent_group),length(potent_group)) <=1.15){
				return(c(order_index[1:index_matrix[,2][1]],potent_index))
			}else if (.tukeybw(potent_group)>=tbw_value*exp(.get_fold_constant(max(potent_group),detect_mod))){
				return(c(order_index[1],potent_index))
			}else{
				return(c(order_index[1:index_matrix[,2][1]]))
			}
		}else{
			return(c(order_index[1:index_matrix[,2][1]]))
		}
	}
}

.SANN_32 <- function(max_2_value,group_vector,detect_mod){
        s_v <- max(group_vector)
        n_s <- length(group_vector)
	fold_constant <- .get_fold_constant(s_v,detect_mod)
	stop_value <- fold_constant* as.numeric(max_2_value)
	order_index <- order(group_vector,decreasing=T)
	if (s_v <= stop_value ){
		out_index <- which(group_vector>=stop_value*0.5)
		return(out_index)
	}else{
		if (s_v>stop_value & length(group_vector)==1){
			return(1)
		}
	}
        e_v <- stop_value
        c_v <- s_v
        n_k <- length( group_vector[group_vector>=e_v] )
        o_k <- length( group_vector[group_vector<e_v] )
        potent_item <- NULL
        iter_n <- 1
        if (n_k == 1){
            potent_item <- iter_n
        }else if(n_k == 0 ){
		return(FALSE)
	}
        prob_e <- .get_pow(s_v/e_v,n_k)
        while(c_v >= e_v & iter_n < n_k+1 & iter_n<=n_s){
		iter_n <- iter_n + 1
                iter_prob <- sqrt(1-(( iter_n-1 )/n_k ))
                c_i <- order_index[iter_n]
		if (iter_n >n_s) {
			break
		}
                c_v <- group_vector[c_i]
                last_c_v <- group_vector[ order_index[iter_n-1] ]
		if (last_c_v == c_v ){
			potent_item <- iter_n
			next
		}
                if (last_c_v/c_v >=fold_constant ){
			potent_item <- iter_n
		}
	}
        if (length(potent_item)){
                return(order_index[1:potent_item])
        }else{
                return(FALSE)
        }
}

.SANN_21 <- function(start_value,group_vector,detect_mod){
	n_s <- length(group_vector)
	lower_cutoff <- sum(group_vector<=1)/n_s
	if (lower_cutoff >=0.1){
		group_vector[group_vector<=1] <- 2^group_vector[group_vector<=1]
	}
	min_v <- min(group_vector)
	fold_constant <- .get_fold_constant(start_value,detect_mod)
	stop_value <- start_value / fold_constant
	if ( stop_value < min_v ){
		return(FALSE)
	}
	order_index <- order(group_vector,decreasing=T)
	c_v <- start_value
	potent_item <- NULL
	iter_n <- 0
	while(c_v > stop_value & iter_n < n_s){
		iter_n <- iter_n + 1
		c_i <- order_index[iter_n]
		c_v <- group_vector[c_i]
		if (iter_n == 1 ){
			last_c_i <- 0
			last_c_v <- start_value
		}else{
			last_c_i <- order_index[iter_n-1]
			last_c_v <- group_vector[last_c_i]
		}
		if (last_c_v /c_v >=fold_constant){
			potent_item <- iter_n-1
		}
	}
	if (length(potent_item)){
		if (potent_item == 0 ){
			return('min_group2')
		}else{
			return(order_index[potent_item:n_s])
		}
	}else{
		return(FALSE)
	}
}
		
.get_T_h <- function( exp_vector ,detect_mod,tbw_value){
	max_h <- max( exp_vector )
	threshold_value <- .get_fold_constant( max_h ,detect_mod)
	e_v <- threshold_value * min(exp_vector)
	n_e <- length( exp_vector )
	n_h <- length( exp_vector[ exp_vector >= e_v ] )
	alive_matrix <- NULL
	if( n_h == 0 ){
		return( FALSE )
#	}else if( n_h == 1 ){
#		alive_matrix <- rbind(alive_matrix,c(1,0))
#		return(list(1,0))
	}else{
		prob_h <- .get_pow( max_h/e_v , n_h )
	}
	o_e <- exp_vector[ exp_vector < e_v ]
	back_flag=0
	if( e_v >= threshold_value*max(o_e) ){
		back_flag=1
	}
	o_k <- length( exp_vector[ exp_vector < e_v ] )
	if( o_k ){
		prob_o <- .get_pow(e_v/min(exp_vector),o_k )
	}
	potential_TS <- NULL
	iter_n <- 0
	try_n <- 1
	record_n <- 0
	order_index <- order( exp_vector , decreasing=T )
    	for( c_i in order_index[-1] ){
		iter_n <- iter_n + 1
		c_v <- exp_vector[c_i]
		iter_prob  <- ifelse(c_v >= e_v,sqrt(1-((iter_n)/n_h )),sqrt(1-((iter_n-n_h)/o_k )))
		last_c_i <- order_index[ iter_n ]
		last_c_v <- exp_vector[last_c_i ]
		if( try_n | c_v >= e_v ){
			if ( c_v < e_v ){ try_n <- try_n - 1 }
 #               	if ( c_v <= max_h/(2*prob_h^(iter_n) ) | last_c_v/c_v >2 ){
#				print(c(c_v,last_c_v,prob_h,max_h/(1.2*prob_h^(iter_n-1)), max_h/(2*prob_h^(iter_n-1) )))
			if( last_c_v / c_v >= 0.75*threshold_value ){
				if( last_c_v / c_v >= threshold_value ){
					alive_matrix <- rbind(alive_matrix,c(last_c_i,iter_n))
				}else{
					potential_TS <- rbind(potential_TS,c(last_c_i,iter_n))
				}
			}
		}
#		if( length(alive_matrix) >3 | length(potential_TS) > 3){break}
	}
	if(length(alive_matrix) ==0 ){
		if ( detect_mod <3 ){
			if (length(potential_TS) == 0){
				return(FALSE)
			}else{
				return(.get_Th_id(order_index,exp_vector,potential_TS,2,tbw_value,detect_mod))
			}
		}else{
			return(FALSE)
		}	
	}else{
		return(.get_Th_id(order_index,exp_vector,alive_matrix,1,tbw_value,detect_mod))
	}
}	

.get_T_l<-function(expression,exp_vector){
	m_l <- min(exp_vector)
	m_o <- setdiff(exp_vector,m_l)



}

.judge_TS <- function(expression,cluster_id,detect_direction,detect_mod,tbw_value){
	names_sample <- names(expression)
	sorted_expression <- sort(expression)
	flag <- 0
	index_3 <- which(cluster_id==3)
	index_2 <- which(cluster_id==2)
	index_1 <- which(cluster_id==1)
	if (length(index_3)>0){
		high_group <- expression[index_3]
	}else{
		high_group <- max(expression)
	}
	if (length(index_1)>0){
		low_group <- expression[index_1]
	}else{
		low_group <- min(expression)
	}
	if (length(index_2)>0){
		max_2 <- max( expression[index_2] )
		min_2 <- min( expression[index_2] )
	}else{
		expression_len <- length(expression)
		max_2_index <- round(expression_len/2)+1
		min_2_index <- round(expression_len/2)
		max_2 <- as.numeric(sorted_expression[max_2_index])
		min_2 <- as.numeric(sorted_expression[min_2_index])
	}
	max_3 <- max( high_group )
	max_1 <- max( low_group )
	old_h_grp <- high_group
	old_l_grp <- low_group
	if (any(cluster_id>4)){
		if ( 32 %in% cluster_id ){
			index_32 <- which( cluster_id==32 )
			cnm_32 <- colnames(expression)[index_32]
			expr_32 <- expression[index_32]
			order_index_32 <- order(expr_32,decreasing=T)
			SANN_32_index <- .SANN_32(max_2,expr_32,detect_mod)
			if (length(SANN_32_index)){
				cnm <- colnames(high_group)
				rnm <- rownames(high_group)
				high_group <- merge(high_group,expr_32[SANN_32_index])
				if (length(SANN_32_index)==1){
					name_SANN_32_index <-cnm_32[SANN_32_index]
					cnm <- append(cnm,name_SANN_32_index)
					colnames(high_group)<-cnm
				}
				rownames(high_group) <- rnm
			}
		}
		if ( 21 %in% cluster_id ){
			index_21 <- which( cluster_id==21 )
			expr_21 <- expression[,index_21]
			order_index_21 <- order(expr_21,decreasing=T)
			SANN_21_index <- .SANN_21(max_2,expr_21,detect_mod)
			if (length(SANN_21_index)>0 & SANN_21_index != FALSE){ 
				rnm <- rownames(low_group)
				if (SANN_21_index == "min_group2"){
					low_group <- merge(expression[index_2][which.min(expression[index_2])],expression[c(index_21,index_1)])
				}else {
					low_group <- merge(low_group,expr_21[SANN_21_index])
				}
				rownames(low_group)<-rnm
			}
		}
	}else{
		high_cutoff <- max_2 - 1e-10
		if (high_cutoff <=0.2){
			high_cutoff <- tbw_value
		}
		high_group <- expression[,expression >=high_cutoff]
		if (length(low_group)<2){
			low_cutoff <- min_2 + 1e-10
			low_group <- expression[,expression <= low_cutoff]
			low_flag <- 1
		}else{
			 low_flag <-0  
		}
	}
	T_h_names <- names(high_group)
	T_l_names <- names(low_group)
	TS_H <-NULL
	TS_L <-NULL
#	print('h')
#	print(high_group)
#	print('o_h')
#	print(old_h_grp)
	TS_H_id <- .get_T_h(high_group,detect_mod,tbw_value)
#	print(TS_H_id)
	.judge_better_percent <- function(expression,x,y,detect_mod){
                TS_H_o <- .percent_judge_4_HTS(expression,x,detect_mod)
                TS_H_h <- .percent_judge_4_HTS(expression,y,detect_mod)
                if (any(TS_H_o==FALSE) & any(TS_H_h==FALSE)){
                        TS_h_names<-FALSE
                }else{
			if (any(TS_H_o==FALSE)){
				return(TS_H_h)
			}else if (any(TS_H_h==FALSE)){
				return(TS_H_o)
			}else{
				if (length(TS_H_o) > length(TS_H_h)){
					return(TS_H_o)
				}else{
					return(TS_H_h)
				}
			}
		}
	}
	if (any(TS_H_id==FALSE)){
		TS_H <- .judge_better_percent(expression,old_h_grp,high_group,detect_mod)
		if ( all(TS_H!=FALSE) & length(TS_H)>0 ){
				TS_h_names <- names(TS_H)
		}else{
			TS_h_names<-FALSE
		}
	}else{
		tmp_TS_h_names<-T_h_names[TS_H_id]
#		tmp_TS_h_names <- T_h_names[order(high_group,decreasing=T)][TS_H_id]
#		print(tmp_TS_h_names)
		tmp_TS_h_exp <- as.vector(high_group[,order(high_group,decreasing=T)[TS_H_id]])
		tmp_TS_h_baseline <- 0.45+0.03*detect_mod
		tmp_TS_h_threshold <- (1/length(expression))*detect_mod
		TS_h_names <- NULL
		sum_try_time <- 4 - detect_mod
		get_sum_index <- 1
		for(i in seq(0,length(tmp_TS_h_exp)-1)){
			tmp_ts_h_index <- i + 1
#			print(tmp_TS_h_exp[tmp_ts_h_index])
#			print(tmp_TS_h_baseline)
#			print(sum(tmp_TS_h_exp))
#			print(sum(tmp_TS_h_exp)*(tmp_TS_h_baseline+tmp_TS_h_threshold*i))
			if (sum(tmp_TS_h_exp[1:tmp_ts_h_index]) >= sum(tmp_TS_h_exp)*(tmp_TS_h_baseline+tmp_TS_h_threshold*i)){
				TS_h_names <- tmp_TS_h_names[1:tmp_ts_h_index]
			}else{
				if (tmp_ts_h_index <= sum_try_time){
					next
				}else{
					break
				}
			}
		}
#		print(TS_h_names)
		if (length(TS_h_names)==0){
			TS_h_names <- FALSE
		}
		if (length(TS_h_names)==1 & length(old_h_grp)>1){
#			print('in length ')
	                tmp_TS_H <- .judge_better_percent(expression,old_h_grp,high_group,detect_mod)
#               }else{
#                       TS_H <- .percent_judge_4_HTS(expression,old_h_grp,detect_mod)
#               }
	                if ( all(tmp_TS_H!=FALSE) & length(tmp_TS_H)>0){
                	                tmp_TS_h_names <- names(tmp_TS_H)
                	    
	                }else{
        	                tmp_TS_h_names<-TS_h_names
                	}
			TS_h_names<-tmp_TS_h_names
		}	
	}
	TS_l_names<-FALSE
	if (detect_direction>1){
		TS_L <- .percent_judge_4_LTS(expression,low_group,tbw_value,detect_mod,low_flag,old_l_grp)
		if (length(TS_L)>0 & class(TS_L)!="logical"){
			TS_l_names <- names(TS_L)
		}
	}else{
		TS_l_names <-FALSE
	}
	putative_result <- rep(0,length(names_sample))
	if (FALSE %in% TS_h_names & FALSE %in% TS_l_names){
		return(putative_result)
	}else{
		if (! FALSE %in% TS_h_names){
			h_index <- match(TS_h_names,names_sample)
			putative_result[h_index]<- 1
		}
		if (! FALSE %in% TS_l_names){
			l_index <- match(TS_l_names,names_sample)
			putative_result[l_index]<- -1
		}
	}
#	if(TS_H==FALSE){
#		potent_HTS<-.percent_judge_4_HTS(expression,high_group,detect_mod)
#
	return(putative_result)
}

.get_TS_result <- function(candidate_TS_matrix,Candidate_TS_G,outdir,draw_plot,draw_heatmap,draw_pca,order_name=rownames(order_TS_P_V)){
	library(ggplot2)
	all_matrix<-t(candidate_TS_matrix)
	all_df <- as.data.frame(all_matrix,row.names=rownames(Candidate_TS_G))
	colnames(all_df)<-colnames(Candidate_TS_G)
	tissue_df <- as.data.frame(candidate_TS_matrix,row.names=colnames(Candidate_TS_G))
	colnames(tissue_df) <- rownames(Candidate_TS_G)
	candicate_tissue_df <- tissue_df[apply(tissue_df,1,function(z){return(any(z!=0))}),]
	TS_df <- all_df[apply(all_df,1,function(z){return(any(z!=0))}),]
	TS_df <- TS_df[order_name,]
	tot_TS_number <- dim(TS_df)
	TS_df_names <- names(TS_df)
	TS_gene_names <- rownames(TS_df)
#	TS_df<-TS_df[order_name,]
	.get_gene_2_ts_function <- function(z){
		z_h <- "-"
		z_h_n <- 0
		z_l <- "-"
		z_l_n <- 0
		if (1 %in% z){
			z_h <- TS_df_names[which(z>0)]
			z_h_n <- length(z_h)
			if(z_h_n>1){
				z_h <- paste(z_h,sep=",",collapse=",")
			}
		}
		if ( -1 %in% z){
			z_l <- TS_df_names[which(z<0)]
			z_l_n <- length(z_l)
			if(z_l_n>1){
				z_l <- paste(z_l,sep=",",collapse=",")
			}
		}
		total_spe_n <- z_h_n + z_l_n
		return(c(z_h,z_h_n,z_l,z_l_n,total_spe_n))
	}
	gene_2_TS_matrix <- t(apply(TS_df,1,.get_gene_2_ts_function))
	.my_table <- function(z){
		k_up<-sum(z>0)
		k_down<-sum(z<0)
		return(c(k_up,k_down))
	}
	sample_TS_freq<-apply(TS_df,2,.my_table)
	rownames(sample_TS_freq)<-c('SEG_H','SEG_L')
	statics_sample_TS <- paste(outdir,"sample_TS_statistics.xls",sep="/")
	cat('Sample\t',paste(colnames(sample_TS_freq),collapse="\t"),"\n",file=statics_sample_TS)
	write.table(sample_TS_freq,file=statics_sample_TS,append=T,quote=F,col.names=F,sep="\t")
	gene_2_TS_df <-data.frame(cbind(TS_gene_names,gene_2_TS_matrix),stringsAsFactors=FALSE)
	ALL_TS_gene <- Candidate_TS_G[gene_2_TS_df$TS_gene_names,]
	colnames(gene_2_TS_df) <- c("SEG_name","H_SEG_sample","H_SEG_sample_number","L_SEG_sample","L_SEG_sample_number","Total_SEG_sample_number")
	gene_2_TS_df<-gene_2_TS_df[order_name,]
	TS_gene_sample_detail <- paste(outdir,"SEG_sample_detail.xls",sep="/")
	write.table(gene_2_TS_df,file=TS_gene_sample_detail,quote=FALSE,sep="\t",row.names=FALSE)
	statics_gene_TS <- paste(outdir,"SEG_in_various_statistics.xls",sep="/")
	out_gene_TS_table <- table(gene_2_TS_df[,6])
	out_gene_TS_names <- names(out_gene_TS_table)
	out_gene_TS_order <- order(as.integer(out_gene_TS_names))
	out_gene_TS_names <- out_gene_TS_names[out_gene_TS_order]
	out_gene_TS_number <- as.integer(out_gene_TS_table[out_gene_TS_order])
	out_gene_TS_matrix <- rbind(c('SEG_in_samples',out_gene_TS_names),c("Gene_number",out_gene_TS_number))
	write.table(out_gene_TS_matrix,file=statics_gene_TS,sep="\t",quote=F,row.names=F,col.names=F)
	png(paste(outdir,"SEG_in_various_statistics.png",sep="/"))
	barplot(out_gene_TS_number,names.arg=out_gene_TS_names,ylim=c(0,max(out_gene_TS_number)*1.2),col="blue",xlab="SEG_in_samples",ylab="Gene_number",main="SEG expressed in diffirent conditions/samples")
	dev.off()
	if (draw_heatmap){
		heatmap_prefix <- paste(outdir,"SEGtool_result_heatmap",sep="/")
		.get_heatmap(ALL_TS_gene,heatmap_prefix)
	}
	if (draw_pca){
		pca_prefix <- paste(outdir,"SEGtool_result_pca",sep="/")
		.get_pca(ALL_TS_gene,pca_prefix)
	}
	only_TS_hdf <- TS_df[apply(TS_df,1,function(z){return(all(z>-1))}),]
	only_TS_hdf_number <- dim(only_TS_hdf)[1]
	only_TS_ldf <- TS_df[apply(TS_df,1,function(z){return(all(z<1))}),]
	only_TS_ldf_number <- dim(only_TS_ldf)[1]
	cross_gene <- TS_gene_names[-match(c(rownames(only_TS_hdf),rownames(only_TS_ldf)),TS_gene_names)]
	cross_gene_number <- length(cross_gene)
	spe_TS_hdf <- only_TS_hdf[apply(only_TS_hdf,1,function(z){if(sum(z)==1){return(TRUE)}else{return(FALSE)}}),]
	spe_TS_ldf <- only_TS_ldf[apply(only_TS_ldf,1,function(z){if(sum(z)==-1){return(TRUE)}else{return(FALSE)}}),]
	multi_TS_hdf <- only_TS_hdf[apply(only_TS_hdf,1,function(z){return(sum(z>- 1) > 1 & all(z> -1))}),]
	multi_TS_ldf <- only_TS_ldf[apply(only_TS_ldf,1,function(z){return(sum(z< 0)>1 & all(z< 1))}),]
	multi_TS_hdf_samplename <- NULL
	multi_TS_hdf_gene_name <- NULL
	multi_TS_ldf_samplename <- NULL
	multi_TS_ldf_gene_name <- NULL
	if (length(multi_TS_hdf)>1){
		multi_TS_hdf_samplename <- colnames(multi_TS_hdf)
		multi_TS_hdf_gene_name <- rownames(multi_TS_hdf)
	}
	if (length(multi_TS_ldf)>1){
		multi_TS_ldf_samplename <- colnames(multi_TS_ldf)
		multi_TS_ldf_gene_name <- rownames(multi_TS_ldf)
	}
	single_sample_outdir<-paste(outdir,'single_sample_result',sep="/")
	.makedir(single_sample_outdir)
	for(i in TS_df_names){
		tmp_result_prefix <- paste(single_sample_outdir,i,sep="/")
		tmp_result_txt <- paste(tmp_result_prefix,".SEGresult.txt",sep="")
		cat(paste("######## ",i,' high SEG detail#####',"\n",sep=""),file=tmp_result_txt)
		tmp_ts_h_cache <- TS_gene_names[which(TS_df[,i]>0)]
		if (length(tmp_ts_h_cache)>0){
			tmp_h_cb<- data.frame(TS_high_gene=tmp_ts_h_cache,TS_high_sample=gene_2_TS_df[tmp_ts_h_cache,"H_SEG_sample"],TS_high_number=gene_2_TS_df[tmp_ts_h_cache,"H_SEG_sample_number"],stringsAsFactors=FALSE)
			tmp_h_cb[,2][tmp_h_cb[,2]==i] <- '-'
			cat(paste(colnames(tmp_h_cb),collapse="\t"),"\n",sep="",file=tmp_result_txt,append=T)
			write.table(tmp_h_cb,file=tmp_result_txt,sep="\t",quote=F,row.names=F,col.names=F,append=T)
		}
		tmp_ts_l_cache <- TS_gene_names[which(TS_df[,i]<0)]
		if (length(tmp_ts_l_cache)>0){
			tmp_l_cb<- data.frame(TS_low_gene=tmp_ts_l_cache,TS_low_sample=gene_2_TS_df[tmp_ts_l_cache,"L_SEG_sample"],TS_low_number=gene_2_TS_df[tmp_ts_l_cache,"L_SEG_sample_number"],stringsAsFactors=FALSE)
#			print(gene_2_TS_df[tmp_ts_l_cache,])

			tmp_l_cb[,2][tmp_l_cb[,2]==i] <- '-'
			cat(paste("######## ",i,' low SEG detail#####',sep=""),"\n",file=tmp_result_txt,append=T)
			cat(paste(colnames(tmp_l_cb),collapse="\t"),"\n",sep="",file=tmp_result_txt,append=T)
			write.table(tmp_l_cb,file=tmp_result_txt,sep="\t",quote=F,row.names=F,col.names=F,append=T)
		}
	}
	single_gene_outdir<-paste(outdir,'single_gene_result',sep="/")
	.makedir(single_gene_outdir)
	if (dim(Candidate_TS_G)[2]>10){
		display_flag=1
	}else{
		display_flag=0
	}
	show_colour <- c("#000000","#D55E00","#01DF3A","#E6E6E6")
	.get_gene_png <- function(gene_id,show_colour=show_colour){
		i<-gene_id
                tmp_figure_prefix <- paste(single_gene_outdir,i,sep="/")
                tmp_result_figure <- paste(tmp_figure_prefix,'.expr.png',sep="")
                tmp_gene_exp_df <- Candidate_TS_G[i,]
                Sample_Name <- names(tmp_gene_exp_df)
                Gene_Name <- rownames(tmp_gene_exp_df)
                Gene_Expression <- as.numeric(tmp_gene_exp_df)
                Sample_classification <- as.numeric(TS_df[i,])
		Sample_classification[Sample_classification<0]<-2
                order_index <- order(Gene_Expression,decreasing=T)
                display_TS_index<-which(Sample_classification!=0)
		display_TS_names<- Sample_Name[display_TS_index]
		zero_index <- which(Sample_classification==0)
		order_zero_index <- zero_index[order(-Gene_Expression[zero_index])]
		display_other_exp <- Gene_Expression[zero_index]
		display_other_index <- order(display_other_exp,decreasing=T)
		length_display_TS_index <- length(display_TS_index)
		if (length_display_TS_index>9){
			display_sample_index <- display_TS_index
		}else{
			need_zero_index <- 10 - length_display_TS_index 
			display_sample_index <- c(display_TS_index,order_zero_index[seq(1,need_zero_index)])
		}
                display_cb <- cbind(Sample_Name=Sample_Name[display_sample_index],Gene_Expression=Gene_Expression[display_sample_index],Sample_classification=Sample_classification[display_sample_index])
		if (length(display_sample_index)<length(Gene_Expression)){
			display_cb <-rbind(display_cb,c('other',-0.001,3))
			display_other_box_index <- setdiff(order_index,display_sample_index)
			display_other_box_df <- data.frame(other=Gene_Expression[display_other_box_index])
			
		}
                display_df<-data.frame(display_cb,stringsAsFactors=FALSE)
                display_df$Gene_Expression <- as.numeric(display_df$Gene_Expression)
		tmp_class_factor <- factor(display_df$Sample_classification)
		colour_levels <- levels(tmp_class_factor)
		colour_display <- show_colour[match(colour_levels,c(0,1,2,3))]
		display_df$Sample_classification <- factor(display_df$Sample_classification,levels=colour_levels,labels=colour_display)
                png(tmp_result_figure,800,600)
                tmp_title <- paste(i," expression in diffrent samples",sep="")
		gg_order_exp_index <- order(-display_df$Gene_Expression)
		gg_order_exp <- display_df$Gene_Expression[gg_order_exp_index]
		gg_order_name <- display_df$Sample_Name[gg_order_exp_index]
		TS_gg_x <- match(display_TS_names,gg_order_name)
		TS_gg_y <- gg_order_exp[TS_gg_x]*1.05
		TS_gg_name <- gg_order_name[TS_gg_x]
		p<-ggplot()+geom_point(data=display_df,aes(x=reorder(Sample_Name,-Gene_Expression),y=Gene_Expression,colour=Sample_classification),size=5)+guides(colour=FALSE)+annotate("text",x=TS_gg_x,y=TS_gg_y,label=TS_gg_name,alpha=.5,size=5,fontface="italic",colour="darkred")+ggtitle(tmp_title)+theme(axis.title.x=element_text(angle=0,face="italic",size=rel(1.2),))+scale_colour_manual(values=colour_display)+labs(x="Samples",y="Gene Expression")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2))
		if ('other' %in% display_df[,1]){
			p<-p+geom_boxplot(data=display_other_box_df,aes(x=dim(display_df)[1],y=other),width=.3,notch=F)

		}
		print(p)
                dev.off()
        }
	if (draw_plot){
		null_result <- apply(as.matrix(TS_gene_names),1,.get_gene_png,show_colour=show_colour)
	}
	TS_sample_number <- sum(apply(TS_df,2,function(x){return(any(x!=0))}))
	summary_TS <- c(tot_TS_number,only_TS_hdf_number,only_TS_ldf_number,cross_gene_number)
	return(summary_TS)
	
}		


.get_heatmap <- function(SEG_df,file_name){
        library(pheatmap)
        filename <- paste(file_name,c("pdf","png"),sep=".")
        SEG_dim <- dim(SEG_df)
        SEG_dim[SEG_dim >=50] <- 50
        select_SEG_df <- SEG_df[1:SEG_dim[1],1:SEG_dim[2]]
        for(i in filename){
                if (grepl("\\.png$",i)){
                        png(i,width=800,height=800)
                }else{
                        pdf(i)
                }
                pheatmap(select_SEG_df,scale="row", color = c("green", "white", "red"),border_color=1,legend=F,fontsize=7)
                dev.off()
        }
}
 
.get_pca <- function(SEG_df,file_name){
        filename <- paste(file_name,c("png","pdf"),sep=".")
        N_samp <- dim(SEG_df)[2] - 1
        pc.cr <- princomp(SEG_df, cor = TRUE)
        pc.var <- pc.cr$sdev ** 2
        pc.var.scaled <- 100 * pc.var/sum(pc.var)
        pc.var.scaled <- sprintf("%.2f", pc.var.scaled)
        axis.names <- paste(names(pc.var), pc.var.scaled, sep = "(")
        axis.names <- paste(axis.names, "%)", sep = "")
        load <- loadings(pc.cr)
        pc <- data.frame(load[,1:2])
	for(i in filename){
		if (grepl("\\.png$",i)){
			png(i,width=800,height=480)
		}else{
			 pdf(i)
		}
        	plot(pc,cex.lab=0.8,tck=-0.015,cex.lab=0.8,cex.axis=0.8,mgp=c(1.3,0.4,0),xlab=axis.names[1],ylab=axis.names[2],pch = "*",col = rainbow(N_samp), main = "loadings of PC1 and PC2 in samples")
        	max_xlim <- max(pc[,1])
        	max_ylim <- max(pc[,2])
        	max_xlim <- abs(max_xlim)*0.24 + max_xlim
        	max_ylim <- abs(max_ylim)*0.22 + max_ylim
        	rnam <- rownames(pc)
        	legend(max_xlim,max_ylim, pch = "*", col = rainbow(N_samp), cex = 0.6, legend = rnam,border=2,xpd=T,ncol=2,bty='n')
        	dev.off()
	}
}

.get_visualisation_html<- function(indir="SEGtool_result",page.name="SEGtool_result", page.title="specific expressed gene analysis results", draw_heatmap=FALSE, draw_pca=FALSE){
	outdir <- indir
	library(hwriterPlus)
	.makedir(outdir)
	setwd(outdir)
	p <- newPage('SEGtools_reuslt.html',link.css = 'main.css',title = "SEGtools_reuslt")
#	hwrite('<link rel="stylesheet" type="text/css" href="main.css"></link>',p)
	cat(".high-expression{color:red}.low-expression{color:green}a{text-decoration:none}body{padding-top:10px;padding-left:10%;padding-right:10%;font-family:Times New Roman;background-color:#fff}table{text-align:center; border-collapse:collapse;border-collapse:collapse;margin-left:auto;margin-right:auto;background-color:#f6f6f6;border-color:gray}table td{padding-top:5px;padding-bottom:5px;padding-left:5px;padding-right:5px;border-color:gray}.title{font-family:Times New Roman,serif;font-size:25pt;font-style:normal;font-weight:900;font-variant:small-caps;text-align:center;margin:10pt 2.5%}.subtitle{font-family:Times New Roman,serif;font-size:15pt;font-style:normal;font-weight:900;font-variant:small-caps;text-align:center;margin:10pt 2.5%}.tafig{font-family:Times New Roman,serif;font-size:12pt;font-style:normal;font-weight:500;text-align:center;margin:5pt 2.5%}",file="main.css")
	hwrite('SEGtools analysis report', p, center = TRUE, heading = 1,class="title")
	hwrite(paste("<span class ='subtitle'> report time :",Sys.time(), "</span>",sep = ""),p, center = TRUE, br = TRUE)
	hwrite("\n\n", p, br = TRUE)
	hwrite("Summary Of SEGtools Analysis Result",p,heading=2)
	summary_total_result <- read.table("summary_result.txt",sep="\t",header=T)
	table_number <- 1
	figure_number <- 1
	hwrite(paste("<p class='tafig'><strong>Table ",table_number,":</strong> The survey for the analysis result </p>",sep=""),p,br=FALSE)
	hwrite(summary_total_result,p,table.style='text-align:center; border-collapse:collapse',col.width=rep('75px',7),row.bgcolor='##ffffaa',br=FALSE)
	table_number <- table_number + 1 
	hwrite("SEGs in samples",p,heading=2)
	summary_gene_TS_result <- read.table("SEG_in_various_statistics.xls",sep="\t")
	colnames(summary_gene_TS_result) <-c()
	hwrite(paste("<p class='tafig'><strong>Table ",table_number,":</strong> Number of genes detected as SEGs in different number of samples</p>",sep=""),p,br=F)
	hwrite(summary_gene_TS_result,p,table.style='text-align:center; border-collapse:collapse', col.bgcolor='#ffaaaa',col.width=c('120px',rep('50px',ncol(summary_gene_TS_result-1))))
	table_number <- table_number + 1
	hwrite("",p,br=TRUE)
	hwriteImage("SEG_in_various_statics.png",p,br=TRUE,center=TRUE)
	#hwrite(sprintf("<center><img src='TS_in_various_statics.png' alt='TS_in_various_statics.png'></img></center>"))
	hwrite(paste("<p class='tafig'><strong> Figure ",figure_number,":</strong> Genes detected as SEG in the different number of samples</p>",sep=""),p,br=TRUE)
	figure_number <- figure_number +1
	hwrite("",p,br=TRUE)
	library(svglite)
	library(stringr)
#	par(mar=c(4,3,0.5,0.5),mgp=c(2,0.5,0))
	main_txt <- sprintf("Number of SEGs for each sample")
	svg_name <- "All_Samples_SEGs_barplot.svg"
	out_svg <- svg_name
	hwrite(sprintf("<center><img border=\"0\" src=\"%s\" alt=\"%s\" width=\"720\" height=\"576\" usemap=\"#%s\" hidefocus=\"true\"></img></center>",svg_name,svg_name,svg_name),p,br=TRUE)
	svg_matrix_file <- "sample_TS_statistics.xls"
	data_matrix <- as.matrix(read.table(svg_matrix_file,header=T,sep="\t",row.names=1))
	svglite(file = out_svg, width = 10, height = 8,bg = "white", pointsize=12)
	position<-barplot(data_matrix,ylim=c(0,max(data_matrix)*1.2),main=main_txt,col=c("red","green"),beside=T,axisnames = F,legend.text=c("SEG_high","SEG_low"),args.legend=list("topright",box.col="white",inset=0.02),ylab="SEG number",offset=0)
	x<-(position[1,]+position[2,])/2
	y<- -max(data_matrix)/20
	text(x,y,labels=substr(colnames(data_matrix),1,8),srt=90,adj=0.8,xpd=TRUE,cex=0.7)
	dev.off()
	svg_context <- readLines(svg_name)
#	tmp_pl_file <- paste(out_svg,".pl",sep="")
#	cat(sprintf("#! /usr/bin/perl\n@l=();\n$f = shift;\nopen(IN,$f);\nwhile(<IN>){\nchomp;\nif(/rect\\sx=\'([^\']+)\'\\sy=\'([^\']+)\'\\swidth=\'([^\']+)\'\\sheight=\'([^\']+)\'/){push @l,$1,$2,$3,$4}};print \"@l\\n\"; "),file=tmp_pl_file)
#	t1 <- try(system(sprintf("perl %s %s",tmp_pl_file,out_svg), intern = TRUE))
	names_data_matrix <- colnames(data_matrix)
	tmp_reg_pattern <- "'[0-9.]+'"
	reg_pattern <- '"[0-9.]+"'
	tmp_map_loci <- as.numeric(gsub("'","",unlist(str_extract_all(svg_context[grepl("<rect x='[^']{1,}' y='[^']{1,}' width='[^']{1,}' height='[^']{1,}'",svg_context)],tmp_reg_pattern))))
	map_loci <- as.numeric(gsub('"',"",unlist(str_extract_all(svg_context[grepl('<rect x="[^"]{1,}" y="[^"]{1,}" width="[^"]{1,}" height="[^"]{1,}"',svg_context)],reg_pattern))))
	hwrite(sprintf("<map name=\"%s\" id= \"%s\"> ",svg_name,svg_name),p,br=FALSE)
	loci_n <- 0
	for(j in seq(0,length(names_data_matrix)-1)){
		inner_index <- j + 1
		href_name <- paste('single_sample_result',sprintf("%s.SEGresult.txt",names_data_matrix[inner_index]),sep="/")
		begin_loci <- j*8+1
		end_loci <- inner_index*8
		hist_map_loci <- tmp_map_loci[begin_loci:end_loci]
		hist_map_up_loci <- paste(hist_map_loci[1],hist_map_loci[2],hist_map_loci[1]+hist_map_loci[3],hist_map_loci[2]+hist_map_loci[4],sep=',')
		hist_map_up_loci <- paste('"',hist_map_up_loci,'"',sep="")
		hist_map_down_loci <- paste(hist_map_loci[5],hist_map_loci[6],hist_map_loci[5]+hist_map_loci[7],hist_map_loci[6]+hist_map_loci[8],sep=',')
		hist_map_down_loci <- paste('"',hist_map_down_loci,'"',sep="")
		name_hist_map <- names_data_matrix[inner_index]
		hwrite(sprintf("<area shape=\"rect\" coords=%s href=\"%s\" target=\"_blank\"/>\n",hist_map_up_loci,href_name),p,br=FALSE)
		hwrite(sprintf("<area shape=\"rect\" coords=%s href=\"%s\" target=\"_blank\"/>\n",hist_map_down_loci,href_name),p,br=FALSE)
	}
	hwrite(sprintf("</map>",p,br=TRUE))
	hwrite("",p,br=TRUE)
        hwrite(paste("<p class='tafig'><strong>Figure ",figure_number,":</strong> Genes detected as SEG in each samples, click the bar graph to see more detail about SEGs in that sample</p>",sep=""),p,br=FALSE)
	figure_number <- figure_number + 1
	if (draw_heatmap){
		fig_heatmap_name<-"SEGtool_result_heatmap.png"
		hwrite("SEGs Heatmap profile(TOP 50)",p,heading=2)
		hwrite("Profiles; green means the gene is low SEG in this sample,red represents high SEG ",p,heading=3)
		hwriteImage(fig_heatmap_name,p,center=TRUE)
		hwrite(paste("<p class='tafig'><strong>Figure ",figure_number,":</strong> This heatmap shows the SEGs(50 numbers,x axis) which having the similar expression pattern in samples(y axis), the colors correspond to the SEG level; red menas high SEG, white means the gene has no special expression pattern in the corresponse sample and the blue represense the low SEG pattern. The color bar on the right indicate the SEG levels.</p> ", sep=""),p,br=TRUE)
		figure_number <- figure_number + 1
	}
	if (draw_pca){
		fig_pca_name<-"SEGtool_result_pca.png"
		hwrite("PCA figure of samples with SEG profiles",p,heading=2)
		hwriteImage(fig_pca_name,p,center=TRUE)
		hwrite(paste("<p class = 'tafig'><strong>Figure ",figure_number,":</strong> This PCA figure shows some kind of relationship between samples, The distance of samples shows the simility of samples.</p>",sep=""),p,br=TRUE)
		figure_number <- figure_number + 1
	}
	hwrite("SEGs table result",p,heading=2)
	SEG_result_detaili_file <- "SEG_sample_detail.xls"
	SEG_result_detail <- read.table(SEG_result_detaili_file,sep="\t",header=T)
#	tmp_png_link_name <- rep('expr figure',nrow(SEG_result_detail))
#	apply(matrix(all_pngNames,ncol=1),1,function(x){return(hwrite('Expression Figure', link = x))})
	all_pngNames <- paste(paste('single_gene_result',SEG_result_detail$SEG_name,sep="/"),".expr.png",sep="")
	tmp_png_link_name<-apply(matrix(all_pngNames,ncol=1),1,function(x){return(hwrite('Expression Figure', link = x,target='_blank'))})
	SEG_result_detail$pnglink <- tmp_png_link_name
#	h_link <- hwrite(SEG_result_detail[,"pnglink"],link=all_pngNames,target='_blank',table=FALSE)
#	h_table_link <- hwrite(c("SEG expr link",h_link),dim=c(nrow(SEG_result_detail)+1,1),row.bgcolor="#ffffaa",br=TRUE)
	p_exp_detail <- newPage('SEG_expression_detail.html',link.css = 'main.css',title = "SEG_expression_detail")
	hwrite('SEG expression detail',p_exp_detail,heading=1,class='title')
#	h_table_result<-hwrite(SEG_result_detail[,c(1:6)],br=TRUE,center=TRUE, row.names=FALSE, row.bgcolor='#ffffaa',table.style='border-collapse: collapse; text-align:center;',row.style=list('text-align:center'))
	hwrite(paste("Table 1: All detected SEGs. The number of samples in which corresponce SEG is specially expressed(up, down and all),detail information please click the hyperlinks at the end of each line" , sep=""),p_exp_detail,br=TRUE)
	hwrite(SEG_result_detail,p_exp_detail,border=1,br=TRUE,center=TRUE, row.names=FALSE, row.bgcolor='#ffffaa',table.style='border-collapse: collapse; text-align:center;',row.style=list('text-align:center'))
	hwrite(paste("Page created by the ", hwrite("SEGtool",link="http://bioinfo.life.hust.edu.cn/SEGtool/")," package using hwriterPlus under ",sessionInfo()[1]$R.version$version.string,collapse=""), p_exp_detail, table.style='font-size:8pt')
	closePage(p_exp_detail)
	show_table <- SEG_result_detail[seq(1,20),]
#	print(show_table)
#	show_png_link_name <- rep('detail',nrow(show_table))
#	show_pngNames <- paste(paste('single_gene_result',show_table$SEG_name,sep="/"),".expr.png",sep="")
#	print(show_pngNames)
#	show_png_link_name<-apply(matrix(rownames(show_pngNames),ncol=1),1,function(x){return(hwrite('Expression Figure', link = x,target="_blank"))})
#	show_table$pnglink <- show_png_link_name
#	h_show_link <- hwrite(show_table[,"pnglink"],link=show_pngNames,target='_blank',table=FALSE)
#	h_show_table_link <- hwrite(c("SEG expr link",h_show_link),dim=c(nrow(show_table)+1,1),row.bgcolor="#ffffaa",br=TRUE)
#	h_show_table_result<-hwrite(show_table[,c(1:6)],br=TRUE,center=TRUE, row.names=FALSE, row.bgcolor='#ffffaa',table.style='border-collapse: collapse; text-align:center;',row.style=list('text-align:center'),col.style=list('text-align:center'),br=TRUE)
	hwrite(show_table,p,border=1,br=TRUE,center=TRUE, row.names=FALSE, row.bgcolor='#ffffaa',table.style='border-collapse: collapse; text-align:center;',row.style=list('text-align:center'),col.style=list('text-align:center'))
	hwrite(paste("Table ",table_number, ": part of detected SEGs : The number of samples in which corresponce SEG is specially expressed(up, down and all),detail information please click the hyperlinks at the end of each line" , sep=""),p,br=TRUE)
	hwrite("",p,br=TRUE)
	hwrite("All SEGs detail link ",p,link='SEG_expression_detail.html',target='_blank',table=FALSE)
	hwrite(paste("The corresponding SEGs expression file is: ",SEG_result_detaili_file,sep=""),p,br=TRUE)
	hwrite(paste("Each SEGs expression figure is at: ",paste(indir,"single_gene_result",sep="/"),sep=""),p,br=TRUE)
	hwrite(paste("Each sample's SEGs information is at: ",paste(indir,"single_sample_result",sep="/"),sep=""),p,br=TRUE)
	hwrite(paste("Page created by the ", hwrite("SEGtool",link="http://bioinfo.life.hust.edu.cn/SEGtool/")," package using hwriterPlus under ",sessionInfo()[1]$R.version$version.string,collapse=""), p, table.style='font-size:8pt')
	closePage(p)
}

