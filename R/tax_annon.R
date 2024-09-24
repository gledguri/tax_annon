#' Print as dataframe
#'
#' This function prints a multiple vectors in the form of a dataframe. Mainly for esthetic reasons
#' @param vector_var The variable/vector that you want to print
#' @param wide The number of wide space you want from the first variable to the next 
#' @return The formated output in Console
#' @export
#' @examples
#' pr_df(tax$id,50)
#' pr_df(tax$sequences,50,new_row = T)
#' sliced_sequences <- pr_df(tax$sequence,50) #if you want to slice all the variables in that vector to 50 characters 
#' cat(pr_df(tax$genus_name,30),pr_df(tax$species_name,30),pr_df(tax$taxid,10),sep = "\n")
pr_df <- function(vector_var,wide=40,new_row=F){
	temp <- wide-nchar(vector_var)
	
	if (is.null(dim(vector_var))) {
		for (i in 1:length(temp)) {
			if (temp[i]>=0) {
				temp[i] <- paste(rep(" ",temp[i]), collapse = "")
			}else{
				vector_var[i] <- substr(vector_var[i],0,wide)
				temp[i] <- ""
			}
		}
		if (new_row==T) return(cat(paste0(vector_var,temp),sep = "\n"))
		return(paste0(vector_var,temp))
	}
}

#' Sequence developer
#'
#' This function creates a sequence of numbers from 1:k when i=0 and k+1:2k 
#' This is very handy when you want to slice data into fragments 
#' @param k the length of the sequence of numbers
#' @param i the loop factor
#' @return A vector of a sequence of numbers k-long
#' @export
#' @examples
#' seq_dev(100,0)
#' seq_dev(100,1)
#' for(i in 0:3) vector[i] <- seq_dev(100,i)
seq_dev <- function(k,i) ((k*i)+1):((k*i)+k)

#' Row Paste
#'
#' This function pastes rows one by one
#' @param inn the dataframe where you want the columns should be collapsed
#' @return A vector of pasted materials
#' @export
#' @examples
#' rowpaste(tax[,1:6])
#' metazoans <- rowpaste(tax[,tax$superkindom=="Metazoa"])
rowpaste <- function(inn){
	vvv <- vector(length=length(nrow(inn)))
	for (k in 1:nrow(inn)) {
		vvv[k] <- paste(inn[k,],collapse = " ")
	}
	return(vvv)
}

#' eDNA data to fasta
#'
#' This function creates fasta sequences from edna reads
#' This is very handy when you want to reblast some sequences
#' @param qid The vector containig the query ID
#' @param seq The vector containig the sequences
#' @param output the name of the output file
#' @return A fasta file with extenxtion .fasta
#' @export
#' @examples
#' # edna_to_fasta(edna$id,edna$sequence,"Plate_NANB.fasta")
edna_to_fasta <- function(qid,seq,output="out.fasta"){
	if (!require("seqinr")) {install.packages("seqinr",dependencies = TRUE);require("seqinr")}
	if (!require("utils")) {install.packages("utils",dependencies = TRUE);require("utils")}
	capture.output(
		for (i in 1:length(qid)) writeLines(paste0(">query",qid[i],";","\n",seq[i])),
		append = F,type = "output",file = "temp.fasta")
	st <- read.fasta("temp.fasta")
	write.fasta(st, qid,file.out = output);unlink("temp.fasta")
}

#' Combine tax list
#'
#' This function combines different tax dataframes (each as a different element on a list) into a singular dataframe
#' @param tax_list the list containing multiple tax dataframes. Each tax dataframe needs to contain 7 columns (Kingdom,Phylum,Class,Order,Family,Genus,Species)
#' @return A dataframe where all the taxonomic information is merged with " | ". Example tax_list[[1]] = "Gadus Morhua" & tax_list[[1]] = "Arctogadus glacialis" the result will be "Gadus Morhua | Arctogadus glacialis" 
#' @export
#' @examples
#' comb_tax <- combine_tax_df(tax_list)
#' 
combine_tax_df <- function(tax_list){
	y <- vector("numeric",length(tax_list))
	for (i in 1:length(tax_list)) {
		y[i] <- dim(tax_list[[i]])[1]
	}
	if(sum(!duplicated(y))==1){
		print("row numbers between tax data are the same")
	}
	if(sum(!duplicated(y))>1){
		print("row numbers between tax data are NOT the same");
		print(y)
	}
	for (i in 1:length(tax_list)) {
		y[i] <- dim(tax_list[[i]])[2]
	}
	if(sum(!duplicated(y))==1){
		print("column numbers between tax data are the same")
	}
	if(sum(!duplicated(y))>1){
		print("column numbers between tax data are NOT the same");
		print(y)
	}
	temp_tax <- as.data.frame(matrix(NA,nrow(tax_list[[1]]),ncol(tax_list[[1]])))
	colnames(temp_tax) <- colnames(tax_list[[1]])
	for (j in 1:ncol(tax_list[[1]])) {
		for (i in 1:nrow(tax_list[[1]])){
			q <- vector(length = length(tax_list))
			for (k in 1:length(tax_list)) {
				q[k] <- !is.na(tax_list[[k]][i,j])
			}
			if(sum(q)==1){
				temp_tax[i,j] <- tax_list[[which(q==T)]][i,j]
			}
			if(sum(q)>1){
				n <- which(q==T)
				w <- vector("character",length(n))
				for (m in 1:length(n)) {
					w[m] <- tax_list[[n[m]]][i,j]
				}
				if(sum(!duplicated(w))==1){
					temp_tax[i,j] <- w[1]
				}
				if(sum(!duplicated(w))>1){
					temp_tax[i,j] <- paste(w[!duplicated(w)],collapse=" | ")
				}
			}
		}
	}
	return(temp_tax)
}


#' Inconsistencies in combining tax lists
#'
#' This function creates a dataframe with only the taxa that are not agreeing between blasted and local database results. This function is embedded in tax_blast_compare() function.
#' @param comb_tax The combined tax dataframe created from combine_tax_df()
#' @return A dataframe with only the rows where the taxonomic information has multiple information. Example, "Gadus Morhua | Arctogadus glacialis"
#' @export
#' @examples
#' inco_comb_tax <- incos_comb_tax_df(comb_tax)
#' 
incos_comb_tax_df <- function(comb_tax){
	v <- vector("numeric",nrow(comb_tax))
	for (i in 1:nrow(comb_tax)) {
		v[i] <- sum(grepl("\\ \\|\\ ",comb_tax[i,]))
	}
	comb_tax[v>0,]
}


#' Selecting the right tax info between combined results
#'
#' This function is an interactive script asking for input on chosing the correct taxonomic information when the blasted and the local results are not agreeing
#' @param comb_tax The combined tax dataframe created from combine_tax_df()
#' @return A dataframe where all the taxonomic information is merged with " | ". Example tax_list[[1]] = "Gadus Morhua" & tax_list[[1]] = "Arctogadus glacialis" the result will be "Gadus Morhua | Arctogadus glacialis" 
#' @export
#' @examples
#' comb_tax_compared <- tax_blast_compare(comb_tax)
#' 
tax_blast_compare <- function(comb_tax){
	require(dplyr)
	require(stringr)
	inco_comb_tax <- incos_comb_tax_df(comb_tax)
	
	for (i in 1:nrow(inco_comb_tax)) {
		xx <- inco_comb_tax[i,] %>% select(colnames(.)[grepl("\\ \\|\\ ",.)])
		qq <- grepl("\\ \\|\\ ",inco_comb_tax[i,])
		pp <- as.data.frame(matrix(NA,length(tax_list),ncol(tax_list[[1]])+1)) %>% setNames(c(colnames(tax_list[[1]]),"db"))
		for (p in 1:length(tax_list)) {
			pp[p,] <- cbind(tax_list[[p]][as.numeric(rownames(inco_comb_tax[i,])),],names(tax_list)[p]) %>% setNames(c(colnames(.)[-ncol(.)],"db"))
		}
		print(pp);cat("\n")
		vv <- vector("character",ncol(xx))
		for (j in 1:ncol(xx)) {
			l <- str_split(xx[j],"\\ \\|\\ ") %>% unlist() %>% as.vector()
			questions(names(xx[j]));cat("\n")
			sl <- multi.menu(c(l,"Enter manually"))
			if(length(sl)==1){
				if(sl==(length(l)+1)){
					oo <- readline("Enter manually \n")
					l <- c(l,oo)
				}
				if(sl==0){
					l <- c(l,"")
					sl <- ncol(xx)+1
				}
				vv[j] <- l[sl]
			}
			if(length(sl)>1){
				vv[j] <- paste(l[sl],collapse = " | ")
			}
		}
		inco_comb_tax[i,qq] <- vv
	}
	comb_tax[as.numeric(rownames(inco_comb_tax)),] <- inco_comb_tax %>% replace_with_na_all(condition = ~.x == "")
	return(comb_tax)
}

#' Creates a single vector with all the taxonomic info
#'
#' This function collapses all the taxonomic information (all 7 columns mentioned in combine_tax_df()) into 1 string separated with "/". Can be used on dataframes as well where each row becomes one string of a vector (of length = nrow(df))
#' @param comb_tax The combined tax dataframe created from combine_tax_df()
#' @return A vector of length = nrow(input_df) where all the taxonomic information is concatenated into 1 string separated with "/"
#' @export
#' @examples
#' comb_tax[!duplicated(tax_to_vector(comb_tax[,1:7])),]
#' 
tax_to_vector <- function(comb_tax){
	vv <- vector("character",nrow(comb_tax))
	for (i in 1:nrow(comb_tax)) {
		vv[i] <- paste(comb_tax[i,],collapse = "/")
	}
	return(vv)
}

#' Converts edna dataframe into OTU dataframe
#'
#' This function selects only the columns where the samples are using "taxid" ans "seq_length"
#' @param edna_df The edna dataframe
#' @return OTU dataframe
#' @export
#' @examples
#' otu <- coll_sim2(edna_to_otu(edna),tax_to_vector(comb_tax_compared))
#' colSums(edna_to_otu(edna))
edna_to_otu <- function(edna_df){
	require(dplyr)
	otu <- edna_df %>% select(colnames(.)[(which(colnames(.)=="taxid")+1):(which(colnames(.)=="seq_length")-1)])
	return(otu)
}

#' Annotate blast results
#'
#' This function reads blasts results and automatically or manually annotates the blasted sequences
#' to the best match resutl
#' @param input the blast output dataframe from NCBI
#' @param method The annotation method. 1=Fully manual; 2=Semi-automatic; 3=Fully Automatic
#' @return prints a csv file in your working directory named "manually_annotated_function_outcome.csv"
#' @return returns a dataframe containing query id, species names of query match and % match
#' @export
#' @examples
#' manually_annotate(blast_output,method = 3)
#' @name Packages asd

manually_annotate <- function(input,method=1,skip=T,print_blast=T,query_vector=NA){ #Method 1: "Conservative; Methood 2: "Semiconservative"; Method 3: "Automatic"
	
	if (!require("devtools")) {install.packages("devtools",dependencies = TRUE);require("devtools")}
	if (!require("dplyr")) {install.packages("dplyr",dependencies = TRUE);require("dplyr")}
	if (!require("coda")) {install.packages("coda",dependencies = TRUE);require("coda")}
	if (!require("vctrs")) {install.packages("vctrs",dependencies = TRUE);require("vctrs")}
	if (!require("stringr")) {install.packages("stringr",dependencies = TRUE);require("stringr")}
	
	magic.pident <- function(input,pident){
		input[pident>=99.5] <-             paste0("\033[0;", 32, "m", input[pident>=99.5],"\033[0m")
		input[pident<99.5&pident>=98.5] <- paste0("\033[0;", 36, "m", input[pident<99.5&pident>=98.5],"\033[0m")
		input[pident<98.5&pident>=97.5] <- paste0("\033[0;", 34, "m", input[pident<98.5&pident>=97.5],"\033[0m")
		input[pident<97.5&pident>=96.5] <- paste0("\033[0;", 33, "m", input[pident<97.5&pident>=96.5],"\033[0m")
		input[pident<96.5] <-             paste0("\033[0;", 31, "m", input[pident<96.5],"\033[0m")
		return(input)}
	
	function_col <- function(x = character()) {
		vec_assert(x, character())
		new_vctr(x, class = "vctrs_function_col")}
	
	format.vctrs_function_col <- function(x,...) {
		gsub("function",crayon::red("function"),vec_data(x))}
	
	pr_df <- function(vector_var,wide=40){
		wide=wide
		temp <- wide-nchar(vector_var)
		for (i in 1:length(temp)) {
			if (temp[i]>=0) {
				temp[i] <- paste(rep(" ",temp[i]), collapse = "")
			}else{
				vector_var[i] <- substr(vector_var[i],0,wide)
				temp[i] <- ""
			}
		}
		return(paste0(vector_var,temp))
	}
	
	allduplicates <- function(vector){
		vector%in%unique(vector[duplicated(vector)])
	}
	
	unlist_uneven_list <- function(list){
		v <- vector("numeric")
		for (i in 1:length(list)) {
			v[i] <- length(list[[i]])
		}
		x <- as.data.frame(matrix(NA,length(list),max(v)))
		for (i in 1:length(list)) {
			x[i,1:length(list[[i]])] <- list[[i]]
		}
		return(x)
	}
	
	is.defined <- function(sym) {
		sym <- deparse(substitute(sym))
		env <- parent.frame()
		exists(sym, env)
	}
	
	rowpaste <- function(inn){
		vvv <- vector(length=length(nrow(inn)))
		for (k in 1:nrow(inn)) {
			vvv[k] <- paste(inn[k,],collapse = " ")
		}
		return(vvv)
	}
	
	input <- tibble::tibble(input)
	columns <- c("Ssciname","scommname","qseqid","sseqid","pident",
							 "length","mismatch","gapopen","qcovus","qstart","qend",
							 "sstart","send","evalue","staxids","qlen","qcovs")
	if (!"annotated_tax"%in%colnames(input)) input$annotated_tax <- NA
	if (!"pmatchsel"%in%colnames(input)) input$pmatchsel <- NA
	cat("The data frame should contain these columns \n",columns)
	cat("checking the names of the input column")
	col_missing <- sum(!columns%in%colnames(input))
	if(col_missing>0)cat(columns[!columns%in%colnames(input)]," column is missing")
	input$Ssciname[grep(";",input$Ssciname)] <- gsub(";", " | ", input$Ssciname[grep(";",input$Ssciname)])
	if(!sum(is.na(query_vector))){input <- input[input$qseqid%in%query_vector,]}
	
	if (method==1) {
		j=1 #build the loop
		while (j <= length(unique(input$qseqid))){ #Loop of method 1
			pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
			i <- unique(input$qseqid)[j] #select the unique query ID
			list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
			
			if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 1
			else { #Not skipping method 1
				#Esthetics of the list
				list_of_query_print <- tibble::tibble(list_of_query);
				list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
				list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
				list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
				list_of_query_print$pident <- function_col(list_of_query_print$pident)
				
				#Print
				if (print_blast==T) {
					cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
					print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
				}
				
				#Create a summary list
				summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
				summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
				
				#Menu and choices
				chsl <- c(paste0(summary_of_list$Ssciname),"Undefined","Enter manually") #Create taxa selection options
				cat("\n");questions("Select the subject that your query belongs to ::");cat("\n")
				sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
													 "Unidentified","Enter manually"))
				
				if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa ")}
				#Writing on the fasta file
				input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") #Assigns the selected taxa to the input database
				input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
				
				#Esthetics
				setTxtProgressBar(pb,j);cat("\n") #progress bar
				j <- j+1 #loop mechanism
				# write.csv(input,"manually_annotated_function_outcome.csv",row.names = FALSE)
			} #Not skipping method 1
		} #Loop of method 1
	} #Method 1
	else if (method==2) { #Method 2
		j=1 #build the loop
		while (j <= length(unique(input$qseqid))){ #Loop of method 2
			pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
			i <- unique(input$qseqid)[j] #select the unique query ID
			list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
			
			#Esthetics of the list
			if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 2
			else { #Not skipping method 2
				list_of_query_print <- tibble::tibble(list_of_query);
				list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
				list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
				list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
				list_of_query_print$pident <- function_col(list_of_query_print$pident)
				
				#Print
				if (print_blast==T) {
					cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
					print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
				}
				
				#Create a summary list
				summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
				summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
				
				chsl <- c(paste0(summary_of_list$Ssciname),"Undefined","Enter manually") #Create taxa selection options
				
				match_id <- 100
				l <- summary_of_list$pident_max==match_id
				if (sum(l)==0) {
					l <- summary_of_list$pident_max==max(summary_of_list$pident_max)
					match_id <- max(summary_of_list$pident_max)
				}
				if (sum(l)==1) {#If only 1 subject is 100%
					cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
											magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
											magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
					cat("\n");questions("The algorithm selected ::")
					cat("\n");cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname[l],30),summary_of_list$pident_max[l]),"::",
																magic.pident(pr_df(summary_of_list$pident_max[l],5),summary_of_list$pident_max[l]), 
																magic.pident(pr_df(summary_of_list$pident[l],5),summary_of_list$pident[l]))),sep="\n")
					cat("\n");questions("Do you agree ?")
					slyn <- menu(c("yes","no"))
					if (slyn==1) {#If yes agree to algorithm's suggestion on annotating to the only subject that is 100#
						input$annotated_tax[input$qseqid==i] <- paste(summary_of_list$Ssciname[l])
						input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%summary_of_list$Ssciname[l]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
						
					}#If yes agree to algorithm's suggestion on annotating to the only subject that is 100#
					else if (slyn==2) {#If no agree to algorithm's suggestion on annotating to the only subject that is 100#
						cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
						sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																		 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																		 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
															 "Unidentified","Enter manually"))
						
						if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa ")}
						
						#Writing on the fasta file
						input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") #!!!
						input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
						
					}#If no agree to algorithm's suggestion on annotating to the only subject that is 100#
				}#If only 1 subject is 100%
				
				
				else{#If multiple subjects are 100%
					multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
					colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
					multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus,multiple_sbj$species)),]
					multiple_sbj <- multiple_sbj[,1:2]
					multiple_sbj$comb_name <- rowpaste(multiple_sbj)
					multiple_sbj$pident_max <- NA
					
					#      multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
					#      colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
					if (prod(allduplicates(multiple_sbj$genus))==1) { #If multiple subjects are the same genus
						
						cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
						cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
												magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
												magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
						cat("\n");questions("The algorithm suggests to annotate to genus rank ::");cat("\n")
						cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus)," spp."),30),mean(summary_of_list$pident_max[l])),"::",
												magic.pident(pr_df(mean(summary_of_list$pident_max[l]),5),mean(summary_of_list$pident_max[l])), 
												magic.pident(pr_df(mean(summary_of_list$pident[l]),5),mean(summary_of_list$pident[l])))),sep = "\n")
						cat("\n");questions("Do you agree ?")
						slyn <- menu(c("yes","no"))
						if (slyn==1) {#If yes agree to algorithm's suggestion on annotating to the genus rank
							sl_0 <- paste0(paste(unique(multiple_sbj$genus))," spp.",collapse = " | ")
							#         sl <- paste(sl_0," || ", paste(chsl[grep(strsplit(sl_0," ")[[1]][1],chsl)],collapse = " | "))
							sl <- paste(sl_0," || ", paste(multiple_sbj$comb_name,collapse = " | "))
							input$annotated_tax[input$qseqid==i] <- sl
							
							for (ii in length(multiple_sbj$comb_name)) {multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii],summary_of_list$Ssciname))]}
							input$pmatchsel[input$qseqid==i] <- paste(round((multiple_sbj$pident_max/100),3),collapse = " | ")
							
							#          input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[grep(strsplit(sl_0," ")[[1]][1],chsl)]]/100),3),collapse = " | ")
						}#If yes agree to algorithm's suggestion on annotating to the genus rank
						
						else if (slyn==2) { #If no agree to algorithm's suggestion on annotating to the genus rank
							cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
							sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																			 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																			 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
																 "Unidentified","Enter manually"))
							if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")}
							
							#Writing on the fasta file
							input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") 
							input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
							
						}#If no agree to algorithm's suggestion on annotating to the genus rank
					} #If multiple subjects are the same genus
					else if (prod(allduplicates(multiple_sbj$genus))==0) {#If multiple subjects are not the same genus
						cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
						cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
												magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
												magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
						cat("\n");questions("The algorithm suggests to annotate to all subjects with 100% match ::");cat("\n")
						sl <- paste(unique(paste(multiple_sbj$genus,multiple_sbj$species)),collapse = " | ")
						cat(sl);cat("\n")
						cat("\n");questions("Do you agree ?")
						slyn <- menu(c("yes","no"))#If yes agree to algorithm's suggestion on annotating all subjects that are 100
						if (slyn==1) {
							input$annotated_tax[input$qseqid==i] <- sl
							input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%str_trim(unlist(strsplit(sl, "[|]")))]/100),3),collapse = " | ")
						}#If yes agree to algorithm's suggestion on annotating all subjects that are 100
						
						else if (slyn==2) {#If no agree to algorithm's suggestion on annotating all subjects that are 100
							cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
							sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																			 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																			 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
																 "Unidentified","Enter manually"))
							if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")}
							
							#Writing on the fasta file
							input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") 
							input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
						}#If no agree to algorithm's suggestion on annotating all subjects that are 100
					}#If multiple subjects are not the same genus
				}#If multiple subjects are 100%
				
				#Esthetics
				setTxtProgressBar(pb,j);cat("\n") #progress bar
				j <- j+1 #loop mechanism
				# write.csv(input,"manually_annotated_function_outcome.csv",row.names = F)
			} #Not skipping method 2
		} #Loop of method 2
	} #Method 2
	
	else if (method==3) { #Method 3
		j=1 #build the loop
		while (j <= length(unique(input$qseqid))){ #Loop of method 3
			pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
			i <- unique(input$qseqid)[j] #select the unique query ID
			list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
			
			#Esthetics of the list
			if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 2
			else { #Not skipping method 3
				list_of_query_print <- tibble::tibble(list_of_query);
				list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
				list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
				list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
				list_of_query_print$pident <- function_col(list_of_query_print$pident)
				
				#Print
				if (print_blast==T) {
					cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
					print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
				}
				
				#Create a summary list
				summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
				summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
				
				chsl <- summary_of_list$Ssciname
				
				match_id <- 100
				l <- summary_of_list$pident_max==match_id
				if (sum(l)==0) {
					l <- summary_of_list$pident_max==max(summary_of_list$pident_max)
					match_id <- max(summary_of_list$pident_max)
				}
				if (sum(l)==1) {#If only 1 subject is 100%
					cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
											magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
											magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
					cat("\n");questions("The algorithm selected ::")
					cat("\n");cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname[l],30),summary_of_list$pident_max[l]),"::",
																magic.pident(pr_df(summary_of_list$pident_max[l],5),summary_of_list$pident_max[l]), 
																magic.pident(pr_df(summary_of_list$pident[l],5),summary_of_list$pident[l]))),sep="\n")
					
					input$annotated_tax[input$qseqid==i] <- paste(summary_of_list$Ssciname[l])
					input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%summary_of_list$Ssciname[l]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
					
				}#If only 1 subject is 100%
				
				else{#If multiple subjects are 100%
					multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
					colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
					multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus,multiple_sbj$species)),]
					multiple_sbj <- multiple_sbj[,1:2]
					multiple_sbj$comb_name <- rowpaste(multiple_sbj)
					multiple_sbj$pident_max <- NA
					if (prod(allduplicates(multiple_sbj$genus))==1) { #If multiple subjects are the same genus
						chsl <- multiple_sbj$comb_name
						cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
						cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
												magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
												magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
						cat("\n");questions("The algorithm suggests to annotate to genus rank ::");cat("\n")
						cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus)," spp."),30),mean(summary_of_list$pident_max[l])),"::",
												magic.pident(pr_df(mean(summary_of_list$pident_max[l]),5),mean(summary_of_list$pident_max[l])), 
												magic.pident(pr_df(mean(summary_of_list$pident[l]),5),mean(summary_of_list$pident[l])))),sep = "\n")
						
						sl_0 <- paste0(paste(unique(multiple_sbj$genus))," spp.",collapse = " | ")
						sl <- paste(sl_0," || ", paste(multiple_sbj$comb_name,collapse = " | "))
						input$annotated_tax[input$qseqid==i] <- sl
						
						for (ii in length(multiple_sbj$comb_name)) {multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii],summary_of_list$Ssciname))]}
						input$pmatchsel[input$qseqid==i] <- paste(round((multiple_sbj$pident_max/100),3),collapse = " | ")
					} #If multiple subjects are the same genus
					
					else if (prod(allduplicates(multiple_sbj$genus))==0) {#If multiple subjects are not the same genus
						cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
						cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,40),summary_of_list$pident_max),"::",
												magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
												magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
						cat("\n");questions("The algorithm suggests to annotate to all subjects with 100% match ::");cat("\n")
						sl <- paste(unique(paste(multiple_sbj$genus,multiple_sbj$species)),collapse = " | ")
						input$annotated_tax[input$qseqid==i] <- sl
						input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%str_trim(unlist(strsplit(sl, "[|]")))]/100),3),collapse = " | ")
					}#If multiple subjects are not the same genus
				}#If multiple subjects are 100%
				
				#Esthetics
				setTxtProgressBar(pb,j);cat("\n") #progress bar
				j <- j+1 #loop mechanism
				# write.csv(input,"manually_annotated_function_outcome.csv",row.names = F)
			} #Not skipping method 3
		} #Loop of method 3
	} #Method 3
	ret_df <- setNames(data.frame(matrix(NA,length(unique(input$qseqid)),3)),c("qseqid","ssciname","pident"))
	ret_df$qseqid <- c(unique(input$qseqid))
	for (ii in unique(input$qseqid)) {
		ret_df[ret_df$qseqid==ii,2] <- unique(input$annotated_tax[input$qseqid==ii])  
		ret_df[ret_df$qseqid==ii,3] <- unique(input$pmatchsel[input$qseqid==ii])  
	}
	write.csv(input,"3_manually_annotated_function_outcome.csv",row.names = F)
	return(ret_df)
}#function


correct_manually_annotate <- function (df_3, print_no=10,df_4,query_seq)
{
	if (!require("devtools")) {
		install.packages("devtools", dependencies = TRUE)
		require("devtools")
	}
	if (!require("dplyr")) {
		install.packages("dplyr", dependencies = TRUE)
		require("dplyr")
	}
	if (!require("coda")) {
		install.packages("coda", dependencies = TRUE)
		require("coda")
	}
	if (!require("vctrs")) {
		install.packages("vctrs", dependencies = TRUE)
		require("vctrs")
	}
	if (!require("stringr")) {
		install.packages("stringr", dependencies = TRUE)
		require("stringr")
	}
	magic.pident <- function(input, pident) {
		input[pident >= 99.5] <- paste0("\033[0;", 32, "m", 
																		input[pident >= 99.5], "\033[0m")
		input[pident < 99.5 & pident >= 98.5] <- paste0("\033[0;", 
																										36, "m", input[pident < 99.5 & pident >= 98.5], 
																										"\033[0m")
		input[pident < 98.5 & pident >= 97.5] <- paste0("\033[0;", 
																										34, "m", input[pident < 98.5 & pident >= 97.5], 
																										"\033[0m")
		input[pident < 97.5 & pident >= 96.5] <- paste0("\033[0;", 
																										33, "m", input[pident < 97.5 & pident >= 96.5], 
																										"\033[0m")
		input[pident < 96.5] <- paste0("\033[0;", 31, "m", input[pident < 
																														 	96.5], "\033[0m")
		return(input)
	}
	function_col <- function(x = character()) {
		vec_assert(x, character())
		new_vctr(x, class = "vctrs_function_col")
	}
	format.vctrs_function_col <- function(x, ...) {
		gsub("function", crayon::red("function"), vec_data(x))
	}
	pr_df <- function(vector_var, wide = 40) {
		wide = wide
		temp <- wide - nchar(vector_var)
		for (i in 1:length(temp)) {
			if (temp[i] >= 0) {
				temp[i] <- paste(rep(" ", temp[i]), collapse = "")
			}
			else {
				vector_var[i] <- substr(vector_var[i], 0, wide)
				temp[i] <- ""
			}
		}
		return(paste0(vector_var, temp))
	}
	allduplicates <- function(vector) {
		vector %in% unique(vector[duplicated(vector)])
	}
	unlist_uneven_list <- function(list) {
		v <- vector("numeric")
		for (i in 1:length(list)) {
			v[i] <- length(list[[i]])
		}
		x <- as.data.frame(matrix(NA, length(list), max(v)))
		for (i in 1:length(list)) {
			x[i, 1:length(list[[i]])] <- list[[i]]
		}
		return(x)
	}
	is.defined <- function(sym) {
		sym <- deparse(substitute(sym))
		env <- parent.frame()
		exists(sym, env)
	}
	rowpaste <- function(inn) {
		vvv <- vector(length = length(nrow(inn)))
		for (k in 1:nrow(inn)) {
			vvv[k] <- paste(inn[k, ], collapse = " ")
		}
		return(vvv)
	}
	input <- df3 %>% filter(qseqid%in%query_seq)
	input <- tibble::tibble(input)
	columns <- c("Ssciname", "scommname", "qseqid", "sseqid",
							 "pident", "length", "mismatch", "gapopen", "qcovus",
							 "qstart", "qend", "sstart", "send", "evalue", "bitscore",
							 "staxids", "qlen", "qcovs","annotated_tax","pmatchsel")
	cat("The data frame should contain these columns \n", columns)
	cat("checking the names of the input column")
	col_missing <- sum(!columns %in% colnames(input))
	if (col_missing > 0) 
		cat(columns[!columns %in% colnames(input)], " column is missing")
	input$Ssciname[grep(";", input$Ssciname)] <- gsub(";", " | ", 
																										input$Ssciname[grep(";", input$Ssciname)])
	j = 1
	while (j <= length(unique(input$qseqid))) {
		if(length(unique(input$qseqid))>1){
			pb <- txtProgressBar(1, length(unique(input$qseqid)), style = 3)
		}
		i <- unique(input$qseqid)[j]
		list_of_query <- input[input$qseqid == i, ]
		
		list_of_query_print <- tibble::tibble(list_of_query)
		list_of_query_print$pident <- magic.pident(list_of_query$pident, 
																							 list_of_query$pident)
		list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch, 
																								 list_of_query$pident)
		list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch)
		list_of_query_print$pident <- function_col(list_of_query_print$pident)
		
		#Printing the blasted results for 1 sequence/hash
		{cat("\n")
			questions(paste0("The query \"", i, "\" blast resutls ::"))
			cat("\n")
			print(list_of_query_print[, c(1:2, 5, 7:13)],n=print_no) #Important
			cat("\n")}
		
		summary_of_list <- list_of_query %>% group_by(Ssciname) %>% 
			summarize(pident_max = max(pident), pident = mean(pident)) %>% 
			arrange(desc(pident_max))
		
		chsl <- c(paste0(summary_of_list$Ssciname), 
							"Undefined", "Enter manually")
		match_id <- 100
		l <- summary_of_list$pident_max == match_id
		if (sum(l) == 0) {
			l <- summary_of_list$pident_max == max(summary_of_list$pident_max)
			match_id <- max(summary_of_list$pident_max)
		}
		if (sum(l) == 1) {
			
			{#Printing the summary of the blast result
				cat("\n")
				questions("A summary of subjects that your query has been blasted towards ::")
				cat("\n")
				cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname, 
																			 30), summary_of_list$pident_max), "::", 
										magic.pident(pr_df(summary_of_list$pident_max, 
																			 5), summary_of_list$pident_max), magic.pident(pr_df(summary_of_list$pident, 
																			 																										5), summary_of_list$pident))), sep = "\n")
				cat("\n")
				cat("\n")}
			
			{#Printing the previous selection
				questions("You previously selected ::")
				cat("\n")
				cat(c(paste(magic.pident(pr_df(unique(input$annotated_tax[input$qseqid==i]),
																			 30), as.numeric(unique(input$pmatchsel[input$qseqid==i]))*100), "::", 
										magic.pident(pr_df(as.numeric(unique(input$pmatchsel[input$qseqid==i]))*100, 
																			 5), as.numeric(unique(input$pmatchsel[input$qseqid==i]))*100))), sep = "\n")
				cat("\n")}
			
			
			cat("\n")
			questions("Would you like to keep the same decision ?")
			slyn <- menu(c("yes", "no"))
			if (slyn == 2) {
				cat("\n")
				questions("Select the taxa that the query should be annotated or enter the name manually ::")
				cat("\n")
				sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname, 
																										30), summary_of_list$pident_max), "::", 
																 magic.pident(pr_df(summary_of_list$pident_max, 
																 									 5), summary_of_list$pident_max), magic.pident(pr_df(summary_of_list$pident, 
																 									 																										5), summary_of_list$pident)), "Unidentified", 
													 "Enter manually"))
				if (sum(chsl[sl] == "Enter manually") > 0) {
					chsl[sl] <- readline("Enter the name of taxa ")
				}
				input$annotated_tax[input$qseqid == i] <- paste(chsl[sl], collapse = " | ")
				input$pmatchsel[input$qseqid == i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname %in% 
																																										 	chsl[sl]]/100, 3), collapse = " | ")
			}
		}
		else {
			multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], 
																																					 "[|]"))), "[ ]+"))
			colnames(multiple_sbj) <- c("genus", "species", 
																	rep("sub.sp", ncol(multiple_sbj) - 2))
			multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus, 
																										 multiple_sbj$species)), ]
			multiple_sbj <- multiple_sbj[, 1:2]
			multiple_sbj$comb_name <- rowpaste(multiple_sbj)
			multiple_sbj$pident_max <- NA
			if (prod(allduplicates(multiple_sbj$genus)) == 1) {
				cat("\n")
				questions("A summary of subjects that your query has been blasted towards ::")
				cat("\n")
				cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname, 
																			 30), summary_of_list$pident_max), "::", 
										magic.pident(pr_df(summary_of_list$pident_max, 
																			 5), summary_of_list$pident_max), magic.pident(pr_df(summary_of_list$pident, 
																			 																										5), summary_of_list$pident))), sep = "\n")
				cat("\n")
				cat("\n")
				questions("The algorithm suggests to annotate to genus rank ::")
				cat("\n")
				cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus), 
																							" spp."), 30), mean(summary_of_list$pident_max[l])), 
										"::", magic.pident(pr_df(mean(summary_of_list$pident_max[l]), 
																						 5), mean(summary_of_list$pident_max[l])), 
										magic.pident(pr_df(mean(summary_of_list$pident[l]), 
																			 5), mean(summary_of_list$pident[l])))), 
						sep = "\n")
				cat("\n")
				questions("Do you agree ?")
				slyn <- menu(c("yes", "no"))
				if (slyn == 1) {
					sl_0 <- paste0(paste(unique(multiple_sbj$genus)), 
												 " spp.", collapse = " | ")
					sl <- paste(sl_0, " || ", paste(multiple_sbj$comb_name, 
																					collapse = " | "))
					input$annotated_tax[input$qseqid == i] <- sl
					for (ii in length(multiple_sbj$comb_name)) {
						multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii], 
																																					 summary_of_list$Ssciname))]
					}
					input$pmatchsel[input$qseqid == i] <- paste(round((multiple_sbj$pident_max/100), 
																														3), collapse = " | ")
				}
				else if (slyn == 2) {
					cat("\n")
					questions("Select the taxa that the query should be annotated or enter the name manually ::")
					cat("\n")
					sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname, 
																											30), summary_of_list$pident_max), "::", 
																	 magic.pident(pr_df(summary_of_list$pident_max, 
																	 									 5), summary_of_list$pident_max), magic.pident(pr_df(summary_of_list$pident, 
																	 									 																										5), summary_of_list$pident)), "Unidentified", 
														 "Enter manually"))
					if (sum(chsl[sl] == "Enter manually") > 
							0) {
						chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")
					}
					input$annotated_tax[input$qseqid == i] <- paste(chsl[sl], 
																													collapse = " | ")
					input$pmatchsel[input$qseqid == i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname %in% 
																																											 	chsl[sl]]/100, 3), collapse = " | ")
				}
			}
			else if (prod(allduplicates(multiple_sbj$genus)) == 0) {
				cat("\n")
				questions("A summary of subjects that your query has been blasted towards ::")
				cat("\n")
				cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname, 
																			 30), summary_of_list$pident_max), "::", 
										magic.pident(pr_df(summary_of_list$pident_max, 
																			 5), summary_of_list$pident_max), magic.pident(pr_df(summary_of_list$pident, 
																			 																										5), summary_of_list$pident))), sep = "\n")
				cat("\n")
				cat("\n")
				questions("You previously selected ::")
				cat("\n")
				sl <- unique(input$annotated_tax[input$qseqid==i])
				# sl <- paste(unique(paste(multiple_sbj$genus, 
				# 												 multiple_sbj$species)), collapse = " | ")
				cat(sl)
				cat("\n")
				cat("\n")
				questions("Would you like to keep the same decision ?")
				slyn <- menu(c("yes", "no"))
				if (slyn == 1) {
					input$annotated_tax[input$qseqid == i] <- sl
					input$pmatchsel[input$qseqid == i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname %in% 
																																													str_trim(unlist(strsplit(sl, "[|]")))]/100), 
																														3), collapse = " | ")
				}
				else if (slyn == 2) {
					cat("\n")
					questions("Select the taxa that the query should be annotated or enter the name manually ::")
					cat("\n")
					sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname, 
																											30), summary_of_list$pident_max), "::", 
																	 magic.pident(pr_df(summary_of_list$pident_max, 
																	 									 5), summary_of_list$pident_max), magic.pident(pr_df(summary_of_list$pident, 
																	 									 																										5), summary_of_list$pident)), "Unidentified", 
														 "Enter manually"))
					if (sum(chsl[sl] == "Enter manually") > 0) {
						chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")
					}
					input$annotated_tax[input$qseqid == i] <- paste(chsl[sl], 
																													collapse = " | ")
					input$pmatchsel[input$qseqid == i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname %in% 
																																											 	chsl[sl]]/100, 3), collapse = " | ")
				}
			}
		}
		if(length(unique(input$qseqid))>1){
			setTxtProgressBar(pb, j)
		}
		cat("\n")
		j <- j + 1
	}
	ret_df <- setNames(data.frame(matrix(NA, length(unique(input$qseqid)), 
																			 3)), c("qseqid", "ssciname", "pident"))
	ret_df$qseqid <- c(unique(input$qseqid))
	for (ii in unique(input$qseqid)) {
		ret_df[ret_df$qseqid == ii, 2] <- unique(input$annotated_tax[input$qseqid == 
																																 	ii])
		ret_df[ret_df$qseqid == ii, 3] <- unique(input$pmatchsel[input$qseqid == 
																														 	ii])
	}
	# write.csv(input, "manually_annotated_function_outcome.csv", 
	# 					row.names = F)
	df3[df3$qseqid%in%query_seq,] <- input
	df4[df4$qseqid%in%query_seq,] <- ret_df
	
	return(list(df_3=df3,
							df_4=df4))
}


#' Annotate blast results
#'
#' This function reads blasts results and automatically or manually annotates the blasted sequences
#' to the best match resutl
#' @param input the blast output dataframe from NCBI
#' @param method The annotation method. 1=Fully manual; 2=Semi-automatic; 3=Fully Automatic
#' @return prints a csv file in your working directory named "manually_annotated_function_outcome.csv"
#' @return returns a dataframe containing query id, species names of query match and % match
#' @export
#' @examples
#' manually_annotate(blast_output,method = 3)
#' @name Packages asd
manually_annotate_2 <- function(input,method=1,read_data,sequence_id,reads){ #Method 1: "Conservative; Methood 2: "Semiconservative"; Method 3: "Automatic"
	
	if (!require("devtools")) {install.packages("devtools",dependencies = TRUE);require("devtools")}
	if (!require("ggu.base.fun")) {devtools::install_github("gledguri/ggednasd/ggu.base.fun")}
	if (!require("dplyr")) {install.packages("dplyr",dependencies = TRUE);require("dplyr")}
	if (!require("coda")) {install.packages("coda",dependencies = TRUE);require("coda")}
	if (!require("vctrs")) {install.packages("vctrs",dependencies = TRUE);require("vctrs")}
	if (!require("stringr")) {install.packages("stringr",dependencies = TRUE);require("stringr")}
	
	magic.pident <- function(input,pident){
		input[pident>=99.5] <-             paste0("\033[0;", 32, "m", input[pident>=99.5],"\033[0m")
		input[pident<99.5&pident>=98.5] <- paste0("\033[0;", 36, "m", input[pident<99.5&pident>=98.5],"\033[0m")
		input[pident<98.5&pident>=97.5] <- paste0("\033[0;", 34, "m", input[pident<98.5&pident>=97.5],"\033[0m")
		input[pident<97.5&pident>=96.5] <- paste0("\033[0;", 33, "m", input[pident<97.5&pident>=96.5],"\033[0m")
		input[pident<96.5] <-             paste0("\033[0;", 31, "m", input[pident<96.5],"\033[0m")
		return(input)}
	
	function_col <- function(x = character()) {
		vec_assert(x, character())
		new_vctr(x, class = "vctrs_function_col")}
	
	format.vctrs_function_col <- function(x,...) {
		gsub("function",crayon::red("function"),vec_data(x))}
	
	pr_df <- function(vector_var,wide=40){
		wide=wide
		temp <- wide-nchar(vector_var)
		for (i in 1:length(temp)) {
			if (temp[i]>=0) {
				temp[i] <- paste(rep(" ",temp[i]), collapse = "")
			}else{
				vector_var[i] <- substr(vector_var[i],0,wide)
				temp[i] <- ""
			}
		}
		return(paste0(vector_var,temp))
	}
	
	allduplicates <- function(vector){
		vector%in%unique(vector[duplicated(vector)])
	}
	
	unlist_uneven_list <- function(list){
		v <- vector("numeric")
		for (i in 1:length(list)) {
			v[i] <- length(list[[i]])
		}
		x <- as.data.frame(matrix(NA,length(list),max(v)))
		for (i in 1:length(list)) {
			x[i,1:length(list[[i]])] <- list[[i]]
		}
		return(x)
	}
	
	is.defined <- function(sym) {
		sym <- deparse(substitute(sym))
		env <- parent.frame()
		exists(sym, env)
	}
	
	rowpaste <- function(inn){
		vvv <- vector(length=length(nrow(inn)))
		for (k in 1:nrow(inn)) {
			vvv[k] <- paste(inn[k,],collapse = " ")
		}
		return(vvv)
	}
	
	input <- tibble::tibble(input)
	columns <- c("Ssciname","scommname","qseqid","sseqid","pident",
							 "length","mismatch","gapopen","qcovus","qstart","qend",
							 "sstart","send","evalue","staxids","qlen","qcovs")
	if (!"annotated_tax"%in%colnames(input)) input$annotated_tax <- NA
	if (!"pmatchsel"%in%colnames(input)) input$pmatchsel <- NA
	cat("The data frame should contain these columns \n",columns)
	cat("checking the names of the input column")
	col_missing <- sum(!columns%in%colnames(input))
	if(col_missing>0)cat(columns[!columns%in%colnames(input)]," column is missing")
	input$Ssciname[grep(";",input$Ssciname)] <- gsub(";", " | ", input$Ssciname[grep(";",input$Ssciname)])
	# if(!sum(is.na(query_vector))){input <- input[input$qseqid%in%query_vector,]}
	
	summarised_reads <- read_data %>%
		rename(hash = !!rlang::sym(sequence_id)) %>%
		group_by(hash) %>%
		summarise(sum_Reads = sum(.data[[reads]]))
	
	input <- left_join(input,summarised_reads,by=c('qseqid'='hash'))
	
	input <- input %>% arrange(desc(sum_Reads))
	
	if (method==1) {
		j=1 #build the loop
		while (j <= length(unique(input$qseqid))){ #Loop of method 1
			pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
			i <- unique(input$qseqid)[j] #select the unique query ID
			list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
			
			# if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 1
			# else { #Not skipping method 1
			#Esthetics of the list
			list_of_query_print <- tibble::tibble(list_of_query);
			list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
			list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
			list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
			list_of_query_print$pident <- function_col(list_of_query_print$pident)
			
			#Print
			cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
			print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
			
			#Create a summary list
			summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
			summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
			
			#Menu and choices
			chsl <- c(paste0(summary_of_list$Ssciname),"Undefined","Enter manually") #Create taxa selection options
			cat("\n");questions("Select the subject that your query belongs to ::");cat("\n")
			sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
															 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
															 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
												 "Unidentified","Enter manually"))
			
			if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa ")}
			#Writing on the fasta file
			input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") #Assigns the selected taxa to the input database
			input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
			
			#Esthetics
			setTxtProgressBar(pb,j);cat("\n") #progress bar
			j <- j+1 #loop mechanism
			# write.csv(input,"manually_annotated_function_outcome.csv",row.names = FALSE)
			# } #Not skipping method 1
		} #Loop of method 1
	} #Method 1
	
	else if (method==2) { #Method 2
		j=1 #build the loop
		while (j <= length(unique(input$qseqid))){ #Loop of method 2
			pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
			i <- unique(input$qseqid)[j] #select the unique query ID
			list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
			
			#Esthetics of the list
			# if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 2
			# else { #Not skipping method 2
			list_of_query_print <- tibble::tibble(list_of_query);
			list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
			list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
			list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
			list_of_query_print$pident <- function_col(list_of_query_print$pident)
			
			#Print
			printbold(unique(list_of_query_print$sum_Reads),wide = 60,col = 47)
			cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
			print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
			
			#Create a summary list
			summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
			summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
			
			chsl <- c(paste0(summary_of_list$Ssciname),"Undefined","Enter manually") #Create taxa selection options
			
			match_id <- 100
			l <- summary_of_list$pident_max==match_id
			if (sum(l)==0) {
				l <- summary_of_list$pident_max==max(summary_of_list$pident_max)
				match_id <- max(summary_of_list$pident_max)
			}
			if (sum(l)==1) {#If only 1 subject is 100%
				cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
				cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
										magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
										magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
				cat("\n");questions("The algorithm selected ::")
				cat("\n");cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname[l],30),summary_of_list$pident_max[l]),"::",
															magic.pident(pr_df(summary_of_list$pident_max[l],5),summary_of_list$pident_max[l]), 
															magic.pident(pr_df(summary_of_list$pident[l],5),summary_of_list$pident[l]))),sep="\n")
				cat("\n");questions("Do you agree ?")
				slyn <- menu(c("yes","no"))
				if (slyn==1) {#If yes agree to algorithm's suggestion on annotating to the only subject that is 100#
					input$annotated_tax[input$qseqid==i] <- paste(summary_of_list$Ssciname[l])
					input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%summary_of_list$Ssciname[l]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
					
				}#If yes agree to algorithm's suggestion on annotating to the only subject that is 100#
				else if (slyn==2) {#If no agree to algorithm's suggestion on annotating to the only subject that is 100#
					cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
					sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																	 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																	 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
														 "Unidentified","Enter manually"))
					
					if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa ")}
					
					#Writing on the fasta file
					input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") #!!!
					input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
					
				}#If no agree to algorithm's suggestion on annotating to the only subject that is 100#
			}#If only 1 subject is 100%
			
			
			else{#If multiple subjects are 100%
				multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
				colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
				multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus,multiple_sbj$species)),]
				multiple_sbj <- multiple_sbj[,1:2]
				multiple_sbj$comb_name <- rowpaste(multiple_sbj)
				multiple_sbj$pident_max <- NA
				
				#      multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
				#      colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
				if (prod(allduplicates(multiple_sbj$genus))==1) { #If multiple subjects are the same genus
					
					cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
											magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
											magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
					cat("\n");questions("The algorithm suggests to annotate to genus rank ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus)," spp."),30),mean(summary_of_list$pident_max[l])),"::",
											magic.pident(pr_df(mean(summary_of_list$pident_max[l]),5),mean(summary_of_list$pident_max[l])), 
											magic.pident(pr_df(mean(summary_of_list$pident[l]),5),mean(summary_of_list$pident[l])))),sep = "\n")
					cat("\n");questions("Do you agree ?")
					slyn <- menu(c("yes","no"))
					if (slyn==1) {#If yes agree to algorithm's suggestion on annotating to the genus rank
						sl_0 <- paste0(paste(unique(multiple_sbj$genus))," spp.",collapse = " | ")
						#         sl <- paste(sl_0," || ", paste(chsl[grep(strsplit(sl_0," ")[[1]][1],chsl)],collapse = " | "))
						sl <- paste(sl_0," || ", paste(multiple_sbj$comb_name,collapse = " | "))
						input$annotated_tax[input$qseqid==i] <- sl
						
						for (ii in length(multiple_sbj$comb_name)) {multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii],summary_of_list$Ssciname))]}
						input$pmatchsel[input$qseqid==i] <- paste(round((multiple_sbj$pident_max/100),3),collapse = " | ")
						
						#          input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[grep(strsplit(sl_0," ")[[1]][1],chsl)]]/100),3),collapse = " | ")
					}#If yes agree to algorithm's suggestion on annotating to the genus rank
					
					else if (slyn==2) { #If no agree to algorithm's suggestion on annotating to the genus rank
						cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
						sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																		 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																		 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
															 "Unidentified","Enter manually"))
						if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")}
						
						#Writing on the fasta file
						input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") 
						input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
						
					}#If no agree to algorithm's suggestion on annotating to the genus rank
				} #If multiple subjects are the same genus
				else if (prod(allduplicates(multiple_sbj$genus))==0) {#If multiple subjects are not the same genus
					cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
											magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
											magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
					cat("\n");questions("The algorithm suggests to annotate to all subjects with 100% match ::");cat("\n")
					sl <- paste(unique(paste(multiple_sbj$genus,multiple_sbj$species)),collapse = " | ")
					cat(sl);cat("\n")
					cat("\n");questions("Do you agree ?")
					slyn <- menu(c("yes","no"))#If yes agree to algorithm's suggestion on annotating all subjects that are 100
					if (slyn==1) {
						input$annotated_tax[input$qseqid==i] <- sl
						input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%str_trim(unlist(strsplit(sl, "[|]")))]/100),3),collapse = " | ")
					}#If yes agree to algorithm's suggestion on annotating all subjects that are 100
					
					else if (slyn==2) {#If no agree to algorithm's suggestion on annotating all subjects that are 100
						cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
						sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
																		 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
																		 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
															 "Unidentified","Enter manually"))
						if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")}
						
						#Writing on the fasta file
						input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") 
						input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
					}#If no agree to algorithm's suggestion on annotating all subjects that are 100
				}#If multiple subjects are not the same genus
			}#If multiple subjects are 100%
			
			#Esthetics
			setTxtProgressBar(pb,j);cat("\n") #progress bar
			j <- j+1 #loop mechanism
			# write.csv(input,"manually_annotated_function_outcome.csv",row.names = F)
			# } #Not skipping method 2
		} #Loop of method 2
	} #Method 2
	
	else if (method==3) { #Method 3
		j=1 #build the loop
		while (j <= length(unique(input$qseqid))){ #Loop of method 3
			pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
			i <- unique(input$qseqid)[j] #select the unique query ID
			list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
			
			#Esthetics of the list
			# if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 2
			# else { #Not skipping method 3
			list_of_query_print <- tibble::tibble(list_of_query);
			list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
			list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
			list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
			list_of_query_print$pident <- function_col(list_of_query_print$pident)
			
			#Print
			cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
			print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
			
			#Create a summary list
			summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
			summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
			
			chsl <- summary_of_list$Ssciname
			
			match_id <- 100
			l <- summary_of_list$pident_max==match_id
			if (sum(l)==0) {
				l <- summary_of_list$pident_max==max(summary_of_list$pident_max)
				match_id <- max(summary_of_list$pident_max)
			}
			if (sum(l)==1) {#If only 1 subject is 100%
				cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
				cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
										magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
										magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
				cat("\n");questions("The algorithm selected ::")
				cat("\n");cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname[l],30),summary_of_list$pident_max[l]),"::",
															magic.pident(pr_df(summary_of_list$pident_max[l],5),summary_of_list$pident_max[l]), 
															magic.pident(pr_df(summary_of_list$pident[l],5),summary_of_list$pident[l]))),sep="\n")
				
				input$annotated_tax[input$qseqid==i] <- paste(summary_of_list$Ssciname[l])
				input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%summary_of_list$Ssciname[l]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
				
			}#If only 1 subject is 100%
			
			else{#If multiple subjects are 100%
				multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
				colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
				multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus,multiple_sbj$species)),]
				multiple_sbj <- multiple_sbj[,1:2]
				multiple_sbj$comb_name <- rowpaste(multiple_sbj)
				multiple_sbj$pident_max <- NA
				if (prod(allduplicates(multiple_sbj$genus))==1) { #If multiple subjects are the same genus
					chsl <- multiple_sbj$comb_name
					cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
											magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
											magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
					cat("\n");questions("The algorithm suggests to annotate to genus rank ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus)," spp."),30),mean(summary_of_list$pident_max[l])),"::",
											magic.pident(pr_df(mean(summary_of_list$pident_max[l]),5),mean(summary_of_list$pident_max[l])), 
											magic.pident(pr_df(mean(summary_of_list$pident[l]),5),mean(summary_of_list$pident[l])))),sep = "\n")
					
					sl_0 <- paste0(paste(unique(multiple_sbj$genus))," spp.",collapse = " | ")
					sl <- paste(sl_0," || ", paste(multiple_sbj$comb_name,collapse = " | "))
					input$annotated_tax[input$qseqid==i] <- sl
					
					for (ii in length(multiple_sbj$comb_name)) {multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii],summary_of_list$Ssciname))]}
					input$pmatchsel[input$qseqid==i] <- paste(round((multiple_sbj$pident_max/100),3),collapse = " | ")
				} #If multiple subjects are the same genus
				
				else if (prod(allduplicates(multiple_sbj$genus))==0) {#If multiple subjects are not the same genus
					cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
					cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,40),summary_of_list$pident_max),"::",
											magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
											magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
					cat("\n");questions("The algorithm suggests to annotate to all subjects with 100% match ::");cat("\n")
					sl <- paste(unique(paste(multiple_sbj$genus,multiple_sbj$species)),collapse = " | ")
					input$annotated_tax[input$qseqid==i] <- sl
					input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%str_trim(unlist(strsplit(sl, "[|]")))]/100),3),collapse = " | ")
				}#If multiple subjects are not the same genus
			}#If multiple subjects are 100%
			
			#Esthetics
			setTxtProgressBar(pb,j);cat("\n") #progress bar
			j <- j+1 #loop mechanism
			# write.csv(input,"manually_annotated_function_outcome.csv",row.names = F)
			# } #Not skipping method 3
		} #Loop of method 3
	} #Method 3
	
	ret_df <- setNames(data.frame(matrix(NA,length(unique(input$qseqid)),3)),c("qseqid","ssciname","pident"))
	ret_df$qseqid <- c(unique(input$qseqid))
	for (ii in unique(input$qseqid)) {
		ret_df[ret_df$qseqid==ii,2] <- unique(input$annotated_tax[input$qseqid==ii])  
		ret_df[ret_df$qseqid==ii,3] <- unique(input$pmatchsel[input$qseqid==ii])  
	}
	write.csv(input,"3_manually_annotated_function_outcome.csv",row.names = F)
	return(ret_df)
}#function
