merge.snp.data <-
function (x, y, ..., error_amount = 1e+06, replacena = TRUE,
		forcestranduse = FALSE, sort = TRUE, intersected_snps_only = FALSE)
{
	if (error_amount <= 0) {
		stop("error_amount can not be <=0")
	}
	if (replacena != T & replacena != F) {
		stop("replacena can be \"TRUE\" or \"FALSE\" only")
	}
	if (forcestranduse != T & forcestranduse != F) {
		stop("forcestranduse can be TRUE or FALSE only")
	}
	if (!is(x, "snp.data")) {
		stop("Wrong data class: the first argument should be snp.data")
	}
	if (!is(x, "snp.data")) {
		stop("Wrong data class: the second argument should be snp.data")
	}
	which_id_intersect_in_x <- which(is.element(x@idnames, y@idnames) *
					c(1:x@nids) > 0)
	which_id_intersect_in_y <- match(x@idnames, y@idnames, nomatch = -1)
	which_id_intersect_in_y <- which_id_intersect_in_y[which_id_intersect_in_y >
					0]
	num_ids_intersected <- length(which_id_intersect_in_x)
	ids_intersected <- c(which_id_intersect_in_x, which_id_intersect_in_y)
	which_snp_intersect_in_x <- which(is.element(x@snpnames,
					y@snpnames) * c(1:x@nsnps) > 0)
	which_snp_intersect_in_y <- match(x@snpnames, y@snpnames,
			nomatch = -1)
	which_snp_intersect_in_y <- which_snp_intersect_in_y[which_snp_intersect_in_y >
					0]
	monos <- function(data) {
		ss <- strsplit(data, "")
		ff <- function(a) {
			if (a[1] == a[2])
				return(TRUE)
			else return(FALSE)
		}
		return(unlist(lapply(ss, FUN = ff)))
	}

	browser()

	if_mono_x_not_mono_y <- monos(coding(x)[which_snp_intersect_in_x]) &
			!monos(coding(y)[which_snp_intersect_in_y])
	which_mono_x_not_mono_y <- which_snp_intersect_in_x[if_mono_x_not_mono_y]
	if_mono_y_not_mono_x <- monos(coding(y)[which_snp_intersect_in_y]) &
			!monos(coding(x)[which_snp_intersect_in_x])
	which_mono_y_not_mono_x <- which_snp_intersect_in_y[if_mono_y_not_mono_x]
	cdng_x <- coding(x)[which_snp_intersect_in_x[if_mono_x_not_mono_y]]
	cdng_y <- coding(y)[which_snp_intersect_in_y[if_mono_x_not_mono_y]]
	strnd_x <- strand(x)[which_snp_intersect_in_x[if_mono_x_not_mono_y]]
	strnd_y <- strand(y)[which_snp_intersect_in_y[if_mono_x_not_mono_y]]
	jjj <- 1
	for (iii in which_mono_x_not_mono_y) {
		mono <- cdng_x[jjj]
		poly <- cdng_y[jjj]
		monoStrand <- strnd_x[jjj]
		polyStrand <- strnd_y[jjj]
		newC <- translate_mono_coding(mono, poly, monoStrand,
				polyStrand, forcestranduse = forcestranduse)
		jjj <- jjj + 1
		coding(x)[iii] <- newC
	}
	cdng_x <- coding(x)[which_snp_intersect_in_x[if_mono_y_not_mono_x]]
	cdng_y <- coding(y)[which_snp_intersect_in_y[if_mono_y_not_mono_x]]
	strnd_x <- strand(x)[which_snp_intersect_in_x[if_mono_y_not_mono_x]]
	strnd_y <- strand(y)[which_snp_intersect_in_y[if_mono_y_not_mono_x]]
	jjj <- 1
	for (iii in which_mono_y_not_mono_x) {
		mono <- cdng_y[jjj]
		poly <- cdng_x[jjj]
		monoStrand <- strnd_y[jjj]
		polyStrand <- strnd_x[jjj]
		newC <- translate_mono_coding(mono, poly, monoStrand,
				polyStrand, forcestranduse = forcestranduse)
		jjj <- jjj + 1
		coding(y)[iii] <- newC
	}
	if (intersected_snps_only) {
		x_new <- x[, which_snp_intersect_in_x]
		y_new <- y[, which_snp_intersect_in_y]
		results <- merge.snp.data(x = x_new, y = y_new, error_amount,
				replacena, forcestranduse, sort, intersected_snps_only)
		return(results)
	}
	chromosome_logic_vec <- c(as.character(x@chromosome[which_snp_intersect_in_x])) ==
			c(as.character(y@chromosome[which_snp_intersect_in_y]))
	chromosome_logic_vec_factor <- factor(chromosome_logic_vec)
	if (length(levels(chromosome_logic_vec_factor)) > 1) {
		position_of_error_snps_in_x <- which_snp_intersect_in_x[!chromosome_logic_vec]
		position_of_error_snps_in_y <- which_snp_intersect_in_y[!chromosome_logic_vec]
		error_output <- data.frame(snp_positon_in_x = position_of_error_snps_in_x,
				snp_positon_in_y = position_of_error_snps_in_y, snpname_in_x = x@snpnames[position_of_error_snps_in_x],
				snpname_in_y = y@snpnames[position_of_error_snps_in_y],
				chromosome_in_x = as.character(x@chromosome)[position_of_error_snps_in_x],
				chromosome_in_y = as.character(y@chromosome)[position_of_error_snps_in_y],
				codding_in_x = as.character(x@coding[position_of_error_snps_in_x]),
				codding_in_y = as.character(y@coding[position_of_error_snps_in_y]))
		write.table(error_output, file = "merge_snp_data_chromosome_error.txt",
				row.names = F, quote = F)
		non_overlaped_chr_num <- table(chromosome_logic_vec)[1]
		stop(paste("There is difference found in chromosome names for ",
						non_overlaped_chr_num, " intersected SNPs. It is expected that all snps have exactly the same names of chromosome. Capital/lower-case letters are considered as different. The file merge_snp_data_chromosome_error.txt contains information about the SNPs having wrong chromosome number. The column snp_positon_in_x and snp_positon_in_y contains order position in x (first) and y (second) input data set. Check column chromosome_in_x and chromosome_in_y and replace the chromosome names in input data set.\n"),
				sep = "")
	}
	if (length(levels(chromosome_logic_vec_factor)) == 0) {
		cat("There are no intersecting SNPs\n")
		if (intersected_snps_only)
			return(list(data = NULL, id = NULL, snp = NULL))
	}
	num_snps_intersected <- length(which_snp_intersect_in_x)
	snp_intersected <- c(which_snp_intersect_in_x, which_snp_intersect_in_y)
	snp_intersected_names <- x@snpnames[which_snp_intersect_in_x]
	id_amount_in_new_array <- x@nids + y@nids - num_ids_intersected
	byte_amount_for_one_snp_in_new_array <- ceiling((id_amount_in_new_array)/4)
	snps_amount_in_new_array <- x@nsnps + y@nsnps - num_snps_intersected
	byte_number_in_new_array <- byte_amount_for_one_snp_in_new_array *
			(snps_amount_in_new_array)
	alleleID_raw <- alleleID.char2raw()
	alleleID_names_char <- names(alleleID_raw)
	names(alleleID_raw) <- NULL
	alleleID_amount <- length(alleleID_raw)
	alleleID_reverse_raw <- alleleID.revstrand()
	names(alleleID_reverse_raw) <- NULL
	return_val <- .C("fast_merge_C_", as.raw(x@gtps), as.integer(x@nids),
			as.integer(x@nsnps), as.raw(y@gtps), as.integer(y@nids),
			as.integer(y@nsnps), as.integer(num_ids_intersected),
			as.integer(num_snps_intersected), as.integer(snp_intersected),
			as.integer(ids_intersected), as.logical(replacena), x@strand@.Data,
			y@strand@.Data, x@coding@.Data, y@coding@.Data, alleleID_raw,
			alleleID_names_char, as.integer(alleleID_amount), alleleID_reverse_raw,
			as.integer(error_amount), found_error_amount_snp = integer(1),
			found_id_error_amount_id = integer(1), id_position_error = integer(error_amount),
			id_snpposition_error = integer(error_amount), val_x_error = raw(error_amount),
			val_y_error = raw(error_amount), snp_position_error = integer(error_amount),
			snp_x_codding_error = raw(error_amount), snp_y_codding_error = raw(error_amount),
			as.logical(forcestranduse), merged_gtps = raw(byte_number_in_new_array))
	new_array_raw <- return_val$merged_gtps
	found_id_error_amount_id <- return_val$found_id_error_amount_id
	id_position_error <- return_val$id_position_error[1:found_id_error_amount_id]
	id_snpposition_error <- return_val$id_snpposition_error[1:found_id_error_amount_id]
	val_x_error <- as.numeric(return_val$val_x_error[1:found_id_error_amount_id])
	val_y_error <- as.numeric(return_val$val_y_error[1:found_id_error_amount_id])
	val_x_error[val_x_error == 0] = NA
	val_x_error[val_x_error == 1] = 0
	val_x_error[val_x_error == 2] = 1
	val_x_error[val_x_error == 3] = 2
	val_y_error[val_y_error == 0] = NA
	val_y_error[val_y_error == 1] = 0
	val_y_error[val_y_error == 2] = 1
	val_y_error[val_y_error == 3] = 2
	found_error_amount_snp <- return_val$found_error_amount_snp
	snp_position_error <- return_val$snp_position_error[1:found_error_amount_snp]
	snp_x_codding_error <- return_val$snp_x_codding_error[1:found_error_amount_snp]
	snp_y_codding_error <- return_val$snp_y_codding_error[1:found_error_amount_snp]
	dim(new_array_raw) <- c(byte_amount_for_one_snp_in_new_array,
			snps_amount_in_new_array)
	new_array_raw <- new("snp.mx", new_array_raw)
	gc(verbose = FALSE)
	tmp_logic_array <- c()
	tmp_logic_array[c(1:length(y@idnames))] = TRUE
	tmp_logic_array[which_id_intersect_in_y] = FALSE
	y_idnames_without_intersected <- y@idnames[tmp_logic_array]
	ids_array <- c(x@idnames, y_idnames_without_intersected)
	y_male_without_intersected <- y@male[tmp_logic_array]
	male_array <- c(x@male, y_male_without_intersected)
	tmp_logic_array <- c()
	tmp_logic_array[c(1:length(y@snpnames))] = TRUE
	tmp_logic_array[which_snp_intersect_in_y] = FALSE
	y_snps_without_intersected <- y@snpnames[tmp_logic_array]
	snp_array <- c(x@snpnames, y_snps_without_intersected)
	y_chromosome_without_intersected <- as.character(y@chromosome[tmp_logic_array])
	chrom_array <- factor(c(as.character(x@chromosome), y_chromosome_without_intersected))
	y_map_without_intersected <- y@map[tmp_logic_array]
	map_array <- c(x@map, y_map_without_intersected)
	y_coding_without_intersected <- y@coding[tmp_logic_array]
	coding_array <- new("snp.coding", c(x@coding, y_coding_without_intersected))
	y_strand_without_intersected <- y@strand[tmp_logic_array]
	strand_array <- new("snp.strand", c(x@strand, y_strand_without_intersected))
	mearged_set_snp_data <- snp.data(nids = id_amount_in_new_array,
			rawdata = new_array_raw, idnames = ids_array, snpnames = snp_array,
			chromosome = chrom_array, map = map_array, coding = coding_array,
			strand = strand_array, male = male_array)
	id_snpposition_error <- x@snpnames[id_snpposition_error]
	snp_position_error <- y@snpnames[snp_position_error]
	cod2index <- 1:length(alleleID.char2raw())
	names(cod2index) <- alleleID.char2raw()
	alleleID <- alleleID.char2raw()
	snp_y_codding_error <- names(alleleID[cod2index[as.character(snp_y_codding_error)]])
	snp_x_codding_error <- names(alleleID[cod2index[as.character(snp_x_codding_error)]])
	id_position_error <- x@idnames[as.numeric(id_position_error)]
	if (found_error_amount_snp != 0) {
		snp_error_data.frame <- data.frame(snpnames = snp_position_error,
				x = snp_y_codding_error, y = snp_x_codding_error,
				stringsAsFactors = F)
	}
	else {
		snp_error_data.frame <- data.frame(snpnames = "...",
				x = "...", y = "...", stringsAsFactors = F)
	}
	if (found_id_error_amount_id != 0) {
		id_error_data.frame <- data.frame(id = id_position_error,
				snpnames = id_snpposition_error, x = val_x_error,
				y = val_y_error, stringsAsFactors = F)
	}
	else {
		id_error_data.frame <- data.frame(id = "...", snpnames = "...",
				x = "...", y = "...", stringsAsFactors = F)
	}
	output_list <- list(data = mearged_set_snp_data, id = id_error_data.frame,
			snp = snp_error_data.frame)
	cat("...merging finished...\n")
	return(output_list)
}