package My_xgrm::Controller::Action;
use My_xgrm::Controller::Utils;
use Mojo::Base 'Mojolicious::Controller';
use JSON;
use LWP::Simple;
use List::Util qw( min max sum );
use POSIX qw(strftime);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

# Render template "dcGO_sitemap.html.ep"
# Render template "dcGO_hie.html.ep"
sub index {
	my $c = shift;
  	
  	
  	##########################
	## Ontology Information
	my %task_name = (
		"EAgene" => 'Enrichment Analyser (Genes)',
		"EAsnp" => 'Enrichment Analyser (SNPs)',
		"RElyser" => 'Enrichment Analyser (Genomic Regions)',
		"EAdomain" => 'Enrichment Analyser (Protein Domains)',
		"SAgene" => 'Subnetwork Analyser (Genes)',
		"SAsnp" => 'Subnetwork Analyser (SNPs)',
		"SAregion" => 'Subnetwork Analyser (Genomic Regions)',
	);
	my %task_des = (
		"EAgene" => 'Enrichment analysis for genes using ontologies',
		"EAsnp" => 'SNPs linked to genes for enrichment analysis',
		"RElyser" => 'Genomic regions linked to genes for enrichment analysis',
		"EAdomain" => 'Enrichment analysis for protein domains using ontologies',
		"SAgene" => 'Subnetwork analysis for gene-level summary data',
		"SAsnp" => 'SNPs linked to genes for subnetwork analysis',
		"SAregion" => 'Genomic regions linked to genes for subnetwork analysis',
	);
	
	# json_task
	my @data_task;
	foreach my $id (keys %task_name) {
		my $rec;
		$rec->{id}=$id;
		$rec->{name}=$task_name{$id};
		$rec->{description}=$task_des{$id};
			
		push @data_task,$rec;
	}
	print STDERR scalar(@data_task)."\n";
	$c->stash(json_task => encode_json(\@data_task));
  	##########################
  	
  	
  	$c->render();
}


# Render template "EAdomain.html.ep"
sub EAdomain {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
	my $domain_type = $c->param('domain_type') || 'pfam'; # by default: pfam
  	my $domainlist = $c->param('domainlist');
  	my $obo = $c->param('obo') || 'GOMF'; # by default: GOMF
  	
	my $FDR_cutoff = $c->param('FDR_cutoff') || 0.05;
	my $min_overlap = $c->param('min_overlap') || 5;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($domainlist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;

		my $input_filename=$tmpFolder.'/'.'data.Domains.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'EAdomain.Domains.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'EAdomain.Domains.'.$digit16.'.r';
	
		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $domainlist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_openxgr'){
			# mac
			#$placeholder="/Users/hfang/Sites/SVN/github/bigdata_fdb";
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_openxgr";
		}elsif(-e '/var/www/html/bigdata_openxgr'){
			# huawei
			#$placeholder="/var/www/bigdata_fdb";
			$placeholder="/var/www/html/bigdata_openxgr";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# www.genomicsummary.com
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/00*
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", domain.type="", obo="", FDR.cutoff="", min.overlap="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRplus-site
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_fdb"
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_openxgr"
		
		library(tidyverse)
		
		domain.type <- "pfam"
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAdomain_Pfam.txt"
		
		domain.type <- "sf"
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAdomain_SF.txt"
		
		data <- read_delim(input.file, delim="\t", col_names=F) %>% as.data.frame() %>% pull(1)
		FDR.cutoff <- 1
		min.overlap <- 3
		obo <- "GOMF"
	}
	
	# read input file
	data <- read_delim(input.file, delim="\t", col_names=F) %>% as.data.frame() %>% pull(1)
	
	if(FDR.cutoff == "NULL"){
		FDR.cutoff <- 1
	}else{
		FDR.cutoff <- as.numeric(FDR.cutoff)
	}
	
	min.overlap <- as.numeric(min.overlap)

	set <- oRDS(str_c("dcGOdb.SET.", domain.type, "2", obo), placeholder=placeholder)
	background <- set$domain_info %>% pull(id)
	eset <- oSEA(data, set, background, test="fisher", min.overlap=min.overlap)

	if(class(eset)=="eSET"){
		# *_enrichment.txt
		df_eTerm <- eset %>% oSEAextract() %>% filter(adjp < FDR.cutoff)
		df_eTerm %>% write_delim(output.file, delim="\t")
		
		# *_enrichment.xlsx
		output.file.enrichment <- gsub(".txt$", ".xlsx", output.file, perl=T)
		df_eTerm %>% openxlsx::write.xlsx(output.file.enrichment)
		#df_eTerm %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/XGRplus-site/app/examples/EAdomain_enrichment.xlsx")
		
		# Dotplot
		message(sprintf("Drawing dotplot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		gp_dotplot <- df_eTerm %>% mutate(name=str_c(id)) %>% oSEAdotplot(FDR.cutoff=0.05, label.top=5, size.title="Number of domains", label.direction.y=c("left","right","none")[3], colors=c("#95c11f","#026634"))
		output.file.dotplot.pdf <- gsub(".txt$", "_dotplot.pdf", output.file, perl=T)
		#output.file.dotplot.pdf <-  "/Users/hfang/Sites/XGR/XGRplus-site/app/examples/EAdomain_enrichment_dotplot.pdf"
		ggsave(output.file.dotplot.pdf, gp_dotplot, device=cairo_pdf, width=5, height=4)
		output.file.dotplot.png <- gsub(".txt$", "_dotplot.png", output.file, perl=T)
		ggsave(output.file.dotplot.png, gp_dotplot, type="cairo", width=5, height=4)
		
		# Forest plot
		if(1){
			message(sprintf("Drawing forest (%s) ...", as.character(Sys.time())), appendLF=TRUE)
			#zlim <- c(0, -log10(df_eTerm$adjp) %>% max() %>% ceiling())
			zlim <- c(0, -log10(df_eTerm$adjp) %>% quantile(0.95) %>% ceiling())
			gp_forest <- df_eTerm %>% mutate(name=str_c(id)) %>% oSEAforest(top=10, colormap="spectral.top", color.title=expression(-log[10]("FDR")), zlim=zlim, legend.direction=c("auto","horizontal","vertical")[3], sortBy=c("or","none")[1], size.title="Number\nof domains", wrap.width=50)		
			output.file.forestplot.pdf <- gsub(".txt$", "_forest.pdf", output.file, perl=T)
			#output.file.forest.pdf <-  "/Users/hfang/Sites/XGR/XGRplus-site/app/examples/EAdomain_enrichment_forestplot.pdf"
			ggsave(output.file.forestplot.pdf, gp_forest, device=cairo_pdf, width=5, height=3.5)
			output.file.forestplot.png <- gsub(".txt$", "_forest.png", output.file, perl=T)
			ggsave(output.file.forestplot.png, gp_forest, type="cairo", width=5, height=3.5)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRplus-site/pier_app/public
		## but outputs at public/tmp/eV2CG.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		if(1){
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- set$domain_info %>% select(id,level,description) %>% semi_join(tibble(id=data), by="id") %>% set_names(c("Identifier","Level","Description"))
		ls_rmd$min_overlap <- min.overlap
		ls_rmd$xlsx_enrichment <- gsub("public/", "", output.file.enrichment, perl=T)
		ls_rmd$pdf_dotplot <- gsub("public/", "", output.file.dotplot.pdf, perl=T)
		ls_rmd$png_dotplot <- gsub("public/", "", output.file.dotplot.png, perl=T)
		ls_rmd$pdf_forestplot <- gsub("public/", "", output.file.forestplot.pdf, perl=T)
		ls_rmd$png_forestplot <- gsub("public/", "", output.file.forestplot.png, perl=T)
		
		output_dir <- gsub("EAdomain.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_EAdomain.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[1], hightlight="default"), output_dir=output_dir)

		}
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)

# huawei
vec <- list.files(path='/root/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
#vec <- list.files(path='/Users/hfang/Sites/XGR/Fang/R', pattern='.r', full.names=T)
vec <- list.files(path='/Users/hfang/Sites/XGR/OpenXGR/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", domain.type=\"$domain_type\", obo=\"$obo\", FDR.cutoff=\"$FDR_cutoff\", min.overlap=\"$min_overlap\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- EAdomain: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_EAdomain.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_EAdomain.html";
			#public/tmp/RMD_eV2CG.html	
			print STDERR "RMD_EAdomain (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_EAdomain.html";
			#public/tmp/digit16_RMD_EAdomain.html
			print STDERR "RMD_EAdomain (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_EAdomain.html
				print STDERR "RMD_EAdomain (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "GElyser.html.ep"
sub GElyser {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $genelist = $c->param('genelist');
  	my $obo = $c->param('obo') || 'GOMF'; # by default: GOMF
  	
	my $FDR_cutoff = $c->param('FDR_cutoff') || 0.05;
	my $min_overlap = $c->param('min_overlap') || 5;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;

		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'GElyser.Genes.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'GElyser.Genes.'.$digit16.'.r';
	
		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_xgrm'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_xgrm";
		}elsif(-e '/var/www/html/bigdata_xgrm'){
			# huawei
			$placeholder="/var/www/html/bigdata_xgrm";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", obo="", FDR.cutoff="", min.overlap="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRm-site
		#placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_fdb"
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_xgrm"
		
		library(tidyverse)
		
		input.file <- "~/Sites/XGR/XGRm-site/app/examples/eg_GElyser.txt"
		
		data <- read_delim(input.file, delim="\t", col_names=F) %>% as.data.frame() %>% pull(1)
		FDR.cutoff <- 1
		min.overlap <- 3
		obo <- "GOMF"
		
		obo <- "MDODD"
		
	}
	
	# read input file
	data <- read_delim(input.file, delim="\t", col_names=F) %>% as.data.frame() %>% pull(1)
	
	if(FDR.cutoff == "NULL"){
		FDR.cutoff <- 1
	}else{
		FDR.cutoff <- as.numeric(FDR.cutoff)
	}
	
	min.overlap <- as.numeric(min.overlap)

	set <- oRDS(str_c("org.Mm.eg", obo), placeholder=placeholder)
	background <- set$info %>% pull(member) %>% unlist() %>% unique()
	eset <- oSEA(data, set, background, test="fisher", min.overlap=min.overlap)

	if(class(eset)=="eSET"){
		# *_enrichment.txt
		df_eTerm <- eset %>% oSEAextract() %>% filter(adjp < FDR.cutoff)		
		####################
		if(nrow(df_eTerm)==0){
			return(NULL)
		}else{
			df_eTerm %>% write_delim(output.file, delim="\t")
		}
		####################

		# *_enrichment.xlsx
		output.file.enrichment <- gsub(".txt$", ".xlsx", output.file, perl=T)
		df_eTerm %>% openxlsx::write.xlsx(output.file.enrichment)
		#df_eTerm %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/XGRm-site/app/examples/GElyser_enrichment.xlsx")
		
		# Dotplot
		message(sprintf("Drawing dotplot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		gp_dotplot <- df_eTerm %>% mutate(name=str_c(name)) %>% oSEAdotplot(FDR.cutoff=0.05, label.top=5, size.title="Number of genes", label.direction.y=c("left","right","none")[3], colors=c("skyblue","orangered"))
		output.file.dotplot.pdf <- gsub(".txt$", "_dotplot.pdf", output.file, perl=T)
		#output.file.dotplot.pdf <-  "/Users/hfang/Sites/XGR/XGRm-site/app/examples/GElyser_enrichment_dotplot.pdf"
		ggsave(output.file.dotplot.pdf, gp_dotplot, device=cairo_pdf, width=5, height=4)
		output.file.dotplot.png <- gsub(".txt$", "_dotplot.png", output.file, perl=T)
		ggsave(output.file.dotplot.png, gp_dotplot, type="cairo", width=5, height=4)
		
		# Forest plot
		if(1){
			message(sprintf("Drawing forest (%s) ...", as.character(Sys.time())), appendLF=TRUE)
			#zlim <- c(0, -log10(df_eTerm$adjp) %>% max() %>% ceiling())
			zlim <- c(0, -log10(df_eTerm$adjp) %>% quantile(0.95) %>% ceiling())
			gp_forest <- df_eTerm %>% mutate(name=str_c(name)) %>% oSEAforest(top=10, colormap="spectral.top", color.title=expression(-log[10]("FDR")), zlim=zlim, legend.direction=c("auto","horizontal","vertical")[3], sortBy=c("or","none")[1], size.title="Number\nof genes", wrap.width=50)
			
			output.file.forestplot.pdf <- gsub(".txt$", "_forest.pdf", output.file, perl=T)
			#output.file.forest.pdf <-  "/Users/hfang/Sites/XGR/XGRm-site/app/examples/GElyser_enrichment_forestplot.pdf"
			ggsave(output.file.forestplot.pdf, gp_forest, device=cairo_pdf, width=5, height=3.5)
			output.file.forestplot.png <- gsub(".txt$", "_forest.png", output.file, perl=T)
			ggsave(output.file.forestplot.png, gp_forest, type="cairo", width=5, height=3.5)
			
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRm-site/my_xgrm/public
		## but outputs at public/tmp/
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		if(1){
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		gene_info <- oRDS("org.Mm.eg", placeholder=placeholder)
		ls_rmd$data_input <- gene_info$info %>% select(Symbol,description) %>% semi_join(tibble(Symbol=data), by="Symbol") %>% transmute(Genes=Symbol, Description=description)
		ls_rmd$min_overlap <- min.overlap
		ls_rmd$xlsx_enrichment <- gsub("public/", "", output.file.enrichment, perl=T)
		ls_rmd$pdf_dotplot <- gsub("public/", "", output.file.dotplot.pdf, perl=T)
		ls_rmd$png_dotplot <- gsub("public/", "", output.file.dotplot.png, perl=T)
		ls_rmd$pdf_forestplot <- gsub("public/", "", output.file.forestplot.pdf, perl=T)
		ls_rmd$png_forestplot <- gsub("public/", "", output.file.forestplot.png, perl=T)
		
		output_dir <- gsub("GElyser.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_GElyser.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[2], hightlight="default"), output_dir=output_dir)

		}
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)

# huawei
vec <- list.files(path='/root/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
vec <- list.files(path='/Users/hfang/Sites/XGR/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", obo=\"$obo\", FDR.cutoff=\"$FDR_cutoff\", min.overlap=\"$min_overlap\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- EAgene: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_GElyser.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_GElyser.html";
			#public/tmp/RMD_GElyser.html	
			print STDERR "RMD_GElyser (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_GElyser.html";
			#public/tmp/digit16_RMD_GElyser.html
			print STDERR "RMD_GElyser (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_GElyser.html
				print STDERR "RMD_GElyser (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "GSlyser.html.ep"
sub GSlyser {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $genelist = $c->param('genelist');
  	my $network = $c->param('network') || 'STRING_high'; # by default: STRING_highest
	my $subnet_size = $c->param('subnet_size') || 30;
	my $subnet_sig = $c->param('subnet_sig') || 'yes';
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;
		
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'GSlyser.Genes.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'GSlyser.Genes.'.$digit16.'.r';
	
		my $my_input;
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_xgrm'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_xgrm";
		}elsif(-e '/var/www/html/bigdata_xgrm'){
			# huawei
			$placeholder="/var/www/html/bigdata_xgrm";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", network="", subnet.size="", subnet.sig="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRm-site
		#placeholder <- "~/Sites/SVN/github/bigdata_fdb"
		placeholder <- "~/Sites/SVN/github/bigdata_xgrm"
		
		library(tidyverse)
		library(igraph)
		
		input.file <- "~/Sites/XGR/XGRm-site/app/examples/eg_GSlyser.txt"
		data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
		network <- "STRINGmm_fin"
		network <- "KEGGmm_pin"
		subnet.size <- 30
		
		ig <- oRDS("ig.KEGGmm.merged", placeholder=placeholder)
		ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
		
	}
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
	
	# subnet.size
	subnet.size <- as.numeric(subnet.size)
	
	message(sprintf("Performing subnetwork analysis restricted to %d network genes (%s) ...", subnet.size, as.character(Sys.time())), appendLF=TRUE)
	#ig <- oRDS("KEGGmm_pin", placeholder=placeholder)
	ig <- oRDS(network, placeholder=placeholder)
	ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
	
	df_data <- tibble(name=data[,1], pvalue=data[,2]) %>% as.data.frame()
	subg <- oSubneterGenes(df_data, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder)

	if(vcount(subg)>0){
	
		vec <- V(subg)$significance %>% as.numeric()
		vec[vec==0] <- min(vec[vec!=0])
		V(subg)$logP <- -log10(vec)
	
		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
		
		df_subg <- subg %>% oIG2TB("nodes") %>% transmute(Genes=name, Pvalue=as.numeric(significance), Description=description) %>% arrange(Pvalue)
		
		vec <- df_subg$Pvalue
		vec[vec==0] <- min(vec[vec!=0])
		vec <- -log10(vec)
		if(max(vec)<20){
			zlim <- c(0, ceiling(max(vec)))
		}else{
			zlim <- c(0, floor(max(vec)/10)*10)
		}
		
		gp_rating <- oGGnetwork(g=subg, node.label="name", node.label.size=3, node.label.color="black", node.label.alpha=0.95, node.label.padding=0.5, node.label.arrow=0, node.label.force=0.4, node.shape=19, node.xcoord="xcoord", node.ycoord="ycoord", node.color="logP", node.color.title=expression(-log[10]("pvalue")), colormap="spectral.top", zlim=zlim, node.size.range=5, title="", edge.color="steelblue4", edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05)
		
		
		# *_crosstalk.txt
		df_subg %>% write_delim(output.file, delim="\t")
		# *_crosstalk.xlsx
		output.file.crosstalk <- gsub(".txt$", "_crosstalk.xlsx", output.file, perl=T)
		df_subg %>% openxlsx::write.xlsx(output.file.crosstalk)

		# *_crosstalk.pdf *_crosstalk.png
		output.file.crosstalk.pdf <- gsub(".txt$", "_crosstalk.pdf", output.file, perl=T)
		ggsave(output.file.crosstalk.pdf, gp_rating, device=cairo_pdf, width=6, height=6)
		output.file.crosstalk.png <- gsub(".txt$", "_crosstalk.png", output.file, perl=T)
		ggsave(output.file.crosstalk.png, gp_rating, type="cairo", width=6, height=6)
		
		combinedP <- 1
		if(subnet.sig=="yes"){
			subg.sig <- oSubneterGenes(df_data, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder, test.permutation=T, num.permutation=10, respect=c("none","degree")[2], aggregateBy="fishers")
			combinedP <- signif(subg.sig$combinedP, digits=2)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRm-site/my_xgrm/public
		## but outputs at public/tmp/
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		gene_info <- oRDS("org.Mm.eg", placeholder=placeholder)
		ls_rmd$data_input <- gene_info$info %>% select(Symbol,description) %>% inner_join(df_data, by=c("Symbol"="name")) %>% transmute(Genes=Symbol, Pvalue=pvalue, Description=description) %>% arrange(Genes)
		ls_rmd$vcount <- nrow(df_subg)
		ls_rmd$combinedP <- combinedP
		ls_rmd$xlsx_crosstalk <- gsub("public/", "", output.file.crosstalk, perl=T)
		ls_rmd$pdf_crosstalk <- gsub("public/", "", output.file.crosstalk.pdf, perl=T)
		ls_rmd$png_crosstalk <- gsub("public/", "", output.file.crosstalk.png, perl=T)
		
		output_dir <- gsub("GSlyser.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_GSlyser.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[2], hightlight="default"), output_dir=output_dir)
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(igraph)

# huawei
vec <- list.files(path='/root/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
vec <- list.files(path='~/Sites/XGR/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", network=\"$network\", subnet.size=\"$subnet_size\", subnet.sig=\"$subnet_sig\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- GSlyser: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_GSlyser.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_GSlyser.html";
			#public/tmp/RMD_GSlyser.html	
			print STDERR "RMD_GSlyser (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_GSlyser.html";
			#public/tmp/digit16_RMD_GSlyser.html
			print STDERR "RMD_GSlyser (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_GSlyser.html
				print STDERR "RMD_GSlyser (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "SAsnp.html.ep"
sub SAsnp {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $crosslink = $c->param('crosslink') || 'proximity_10000';
  	my $network = $c->param('network') || 'STRING_high'; # by default: STRING_highest
	my $subnet_size = $c->param('subnet_size') || 30;
	my $subnet_sig = $c->param('subnet_sig') || 'yes';
  	
	my $significance_threshold = $c->param('significance_threshold') || 0.05;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;

		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'SAsnp.SNPs.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'SAsnp.SNPs.'.$digit16.'.r';
	
		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_openxgr'){
			# mac
			#$placeholder="/Users/hfang/Sites/SVN/github/bigdata_fdb";
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_openxgr";
		}elsif(-e '/var/www/html/bigdata_openxgr'){
			# huawei
			#$placeholder="/var/www/bigdata_fdb";
			$placeholder="/var/www/html/bigdata_openxgr";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# www.genomicsummary.com
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", population="", crosslink="", significance.threshold="", network="", subnet.size="", subnet.sig="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRplus-site
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_fdb"
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_openxgr"
		
		library(tidyverse)
		library(GenomicRanges)
		library(igraph)
		
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAsnp.txt"
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAsnp_IND.txt"
		
		#####################
		# for eg_SAregion.txt
		data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
		GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)
		ind <- match(data$snp, names(GR.SNP))
		gr <- GR.SNP[ind[!is.na(ind)]]
		df_gr <- gr %>% as.data.frame() %>% as_tibble(rownames="snp") %>% transmute(snp,region=str_c(seqnames,":",start,"-",end))
		df_gr %>% inner_join(data, by="snp") %>% select(region,pvalue) %>% write_delim("~/Sites/XGR/XGRplus-site/app/examples/eg_SAregion.txt", delim="\t")
		#####################
		
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAsnp_IND.txt"
		data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
		
		LD.customised <- oRDS("GWAS_LD.EUR", placeholder=placeholder) %>% as.data.frame()
		significance.threshold=5e-5
		distance.max=20000
		relative.importance=c(1/3,1/3,1/3)
		
		crosslink <- "pQTL_Plasma"
		crosslink <- "PCHiC_PMID27863249_Activated_total_CD4_T_cells"
		crosslink <- "proximity_20000"
		
		network <- "STRING_high"
		network <- "KEGG"
		subnet.size <- 30
	}
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2) %>% set_names("snp","pvalue")
	
	if(significance.threshold == "NULL"){
		significance.threshold <- NULL
	}else{
		significance.threshold <- as.numeric(significance.threshold)
	}
	
	if(population=="NA"){
		LD.customised <- NULL
	}else{
		LD.customised <- oRDS(str_c("GWAS_LD.", population), placeholder=placeholder) %>% as.data.frame()
	}
	
	relative.importance <- c(0,0,0)
	include.QTL <- NULL
	include.RGB <- NULL
	if(str_detect(crosslink, "proximity")){
		relative.importance <- c(1,0,0)
		distance.max <- str_replace_all(crosslink, "proximity_", "") %>% as.numeric()
	}else if(str_detect(crosslink, "QTL")){
		relative.importance <- c(0,1,0)
		include.QTL <- crosslink
	}else if(str_detect(crosslink, "PCHiC")){
		relative.importance <- c(0,0,1)
		include.RGB <- crosslink
	}
	
	#GR.SNP <- oRDS("dbSNP_GWAS", placeholder=placeholder)
	GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)
	GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)
	
	xGene <- oSNP2xGenes(data, score.cap=100, LD.customised=LD.customised, significance.threshold=significance.threshold, distance.max=distance.max, decay.kernel="constant", GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=include.QTL, include.RGB=include.RGB, relative.importance=relative.importance, placeholder=placeholder)
	
	subnet.size <- as.numeric(subnet.size)
	
	message(sprintf("Performing subnetwork analysis restricted to %d network genes (%s) ...", subnet.size, as.character(Sys.time())), appendLF=TRUE)
	ig <- oDefineNet(network=network, STRING.only=c("experimental_score","database_score"), placeholder=placeholder)
	ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
	
	df_data <- tibble(name=xGene$xGene$Gene, pvalue=10^(-xGene$xGene$LScore)) %>% as.data.frame()
	subg <- oSubneterGenes(df_data, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder)

	if(vcount(subg)>0){
		# *_LG.xlsx
		output.file.LG <- gsub(".txt$", "_LG.xlsx", output.file, perl=T)
		df_LG <- xGene$xGene
		df_LG %>% openxlsx::write.xlsx(output.file.LG)
		
		# *_LG_evidence.xlsx
		output.file.LG_evidence <- gsub(".txt$", "_LG_evidence.xlsx", output.file, perl=T)
		df_evidence <- xGene$Evidence %>% transmute(Gene,SNP,SNP_type=ifelse(SNP_Flag=="Lead","Input","LD"),Evidence=Context)
		df_evidence %>% openxlsx::write.xlsx(output.file.LG_evidence)

		vec <- V(subg)$significance %>% as.numeric()
		vec[vec==0] <- min(vec[vec!=0])
		V(subg)$logP <- -log10(vec)

		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
		
		df_subg <- subg %>% oIG2TB("nodes") %>% transmute(Genes=name, Pvalue=as.numeric(significance), Description=description) %>% arrange(Pvalue)
		
		gp_rating <- oGGnetwork(g=subg, node.label="name", node.label.size=3, node.label.color="black", node.label.alpha=0.95, node.label.padding=0.5, node.label.arrow=0, node.label.force=0.4, node.shape=19, node.xcoord="xcoord", node.ycoord="ycoord", node.color="logP", node.color.title="Linked\ngene\nscores", colormap="spectral.top", zlim=c(0,10), node.size.range=5, title="", edge.color="steelblue4", edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05)
		
		
		# *_crosstalk.txt
		df_subg %>% write_delim(output.file, delim="\t")
		# *_crosstalk.xlsx
		output.file.crosstalk <- gsub(".txt$", "_crosstalk.xlsx", output.file, perl=T)
		df_subg %>% openxlsx::write.xlsx(output.file.crosstalk)

		# *_crosstalk.pdf *_crosstalk.png
		output.file.crosstalk.pdf <- gsub(".txt$", "_crosstalk.pdf", output.file, perl=T)
		ggsave(output.file.crosstalk.pdf, gp_rating, device=cairo_pdf, width=6, height=6)
		output.file.crosstalk.png <- gsub(".txt$", "_crosstalk.png", output.file, perl=T)
		ggsave(output.file.crosstalk.png, gp_rating, type="cairo", width=6, height=6)
		
		combinedP <- 1
		if(subnet.sig=="yes"){
			subg.sig <- oSubneterGenes(df_data, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder, test.permutation=T, num.permutation=10, respect=c("none","degree")[2], aggregateBy="fishers")
			combinedP <- signif(subg.sig$combinedP, digits=2)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRplus-site/my_xgrplus/public
		## but outputs at public/tmp/eV2CG.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD %s %f (%s) ...", subnet.sig, combinedP, as.character(Sys.time())), appendLF=TRUE)
		
		if(1){
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		ind <- match(data$snp, names(GR.SNP))
		ls_rmd$data_input <- tibble(SNPs=data$snp[!is.na(ind)], Pvalue=data$pvalue[!is.na(ind)], GR.SNP[ind[!is.na(ind)]] %>% as.data.frame() %>% as_tibble() %>% transmute(Locus=str_c(seqnames,":",start,"-",end)))
		ls_rmd$xlsx_LG <- gsub("public/", "", output.file.LG, perl=T)
		ls_rmd$xlsx_LG_evidence <- gsub("public/", "", output.file.LG_evidence, perl=T)
		ls_rmd$num_LG <- nrow(df_LG)
		ls_rmd$vcount <- nrow(df_subg)
		ls_rmd$combinedP <- combinedP
		ls_rmd$xlsx_crosstalk <- gsub("public/", "", output.file.crosstalk, perl=T)
		ls_rmd$pdf_crosstalk <- gsub("public/", "", output.file.crosstalk.pdf, perl=T)
		ls_rmd$png_crosstalk <- gsub("public/", "", output.file.crosstalk.png, perl=T)
		
		output_dir <- gsub("SAsnp.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_SAsnp.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[1], hightlight="default"), output_dir=output_dir)

		}
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(GenomicRanges)
library(igraph)

# huawei
vec <- list.files(path='/root/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
#vec <- list.files(path='/Users/hfang/Sites/XGR/Fang/R', pattern='.r', full.names=T)
vec <- list.files(path='/Users/hfang/Sites/XGR/OpenXGR/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", crosslink=\"$crosslink\", significance.threshold=\"$significance_threshold\", network=\"$network\", subnet.size=\"$subnet_size\", subnet.sig=\"$subnet_sig\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- SAsnp: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_SAsnp.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_SAsnp.html";
			#public/tmp/RMD_eV2CG.html	
			print STDERR "RMD_SAsnp (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_SAsnp.html";
			#public/tmp/digit16_RMD_SAsnp.html
			print STDERR "RMD_SAsnp (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_SAsnp.html
				print STDERR "RMD_SAsnp (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "EAsnp.html.ep"
sub EAsnp {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $crosslink = $c->param('crosslink') || 'proximity_10000';
  	my $obo = $c->param('obo') || 'GOMF'; # by default: GOMF
  	
	my $significance_threshold = $c->param('significance_threshold') || 0.05;
	my $FDR_cutoff = $c->param('FDR_cutoff') || 0.05;
	my $min_overlap = $c->param('min_overlap') || 5;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;

		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'EAsnp.SNPs.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'EAsnp.SNPs.'.$digit16.'.r';
	
		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_openxgr'){
			# mac
			#$placeholder="/Users/hfang/Sites/SVN/github/bigdata_fdb";
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_openxgr";
		}elsif(-e '/var/www/html/bigdata_openxgr'){
			# huawei
			#$placeholder="/var/www/bigdata_fdb";
			$placeholder="/var/www/html/bigdata_openxgr";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# www.genomicsummary.com
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", population="", crosslink="", significance.threshold="", obo="", FDR.cutoff="", min.overlap="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRplus-site
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_fdb"
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_openxgr"
		
		library(tidyverse)
		library(GenomicRanges)
		
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAsnp.txt"
		input.file <- "~/Sites/XGR/XGRplus-site/app/examples/eg_EAsnp_IND.txt"
		
		data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
		
		LD.customised <- oRDS("GWAS_LD.EUR", placeholder=placeholder) %>% as.data.frame()
		significance.threshold <- 5e-5
		distance.max <- 2000
		relative.importance <- c(1/3,1/3,1/3)
		
		crosslink <- "pQTL_Plasma"
		crosslink <- "PCHiC_PMID27863249_Activated_total_CD4_T_cells"
		
		FDR.cutoff <- 0.05
		min.overlap <- 3
		obo <- "GOMF"
		
		obo <- "KEGGEnvironmentalOrganismal"
	}
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2) %>% set_names("snp","pvalue")
	
	if(significance.threshold == "NULL"){
		significance.threshold <- 1
	}else{
		significance.threshold <- as.numeric(significance.threshold)
	}
	
	if(FDR.cutoff == "NULL"){
		FDR.cutoff <- 1
	}else{
		FDR.cutoff <- as.numeric(FDR.cutoff)
	}
	
	if(population=="NA"){
		LD.customised <- NULL
	}else{
		LD.customised <- oRDS(str_c("GWAS_LD.", population), placeholder=placeholder) %>% as.data.frame()
	}
	
	relative.importance <- c(0,0,0)
	include.QTL <- NULL
	include.RGB <- NULL
	if(str_detect(crosslink, "proximity")){
		relative.importance <- c(1,0,0)
		distance.max <- str_replace_all(crosslink, "proximity_", "") %>% as.numeric()
	}else if(str_detect(crosslink, "QTL")){
		relative.importance <- c(0,1,0)
		include.QTL <- crosslink
	}else if(str_detect(crosslink, "PCHiC")){
		relative.importance <- c(0,0,1)
		include.RGB <- crosslink
	}
	
	GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)
	GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)
	
	xGene <- oSNP2xGenes(data, score.cap=100, LD.customised=LD.customised, significance.threshold=significance.threshold, distance.max=distance.max, decay.kernel="constant", GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=include.QTL, include.RGB=include.RGB, relative.importance=relative.importance, placeholder=placeholder)
	
	min.overlap <- as.numeric(min.overlap)
	set <- oRDS(str_c("org.Hs.eg", obo), placeholder=placeholder)
	background <- set$info %>% pull(member) %>% unlist() %>% unique()
	eset <- oSEA(xGene$xGene %>% pull(Gene), set, background, test="fisher", min.overlap=min.overlap)

	if(class(eset)=="eSET"){
		# *_LG.xlsx
		output.file.LG <- gsub(".txt$", "_LG.xlsx", output.file, perl=T)
		df_LG <- xGene$xGene
		df_LG %>% openxlsx::write.xlsx(output.file.LG)
		
		# *_LG_evidence.xlsx
		output.file.LG_evidence <- gsub(".txt$", "_LG_evidence.xlsx", output.file, perl=T)
		df_evidence <- xGene$Evidence %>% transmute(SNP,SNP_type=ifelse(SNP_Flag=="Lead","Input","LD"),Gene,Evidence=Context)
		df_evidence %>% openxlsx::write.xlsx(output.file.LG_evidence)

		# *_enrichment.txt
		df_eTerm <- eset %>% oSEAextract() %>% filter(adjp < FDR.cutoff)
		####################
		if(nrow(df_eTerm)==0){
			return(NULL)
		}else{
			df_eTerm %>% write_delim(output.file, delim="\t")
		}
		####################
		
		# *_enrichment.xlsx
		output.file.enrichment <- gsub(".txt$", ".xlsx", output.file, perl=T)
		df_eTerm %>% openxlsx::write.xlsx(output.file.enrichment)
		#df_eTerm %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/XGRplus-site/app/examples/EAsnp_enrichment.xlsx")
		
		# Dotplot
		message(sprintf("Drawing dotplot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		gp_dotplot <- df_eTerm %>% mutate(name=str_c(name)) %>% oSEAdotplot(FDR.cutoff=0.05, label.top=5, size.title="Number of genes", label.direction.y=c("left","right","none")[3], colors=c("#95c11f","#026634"))
		output.file.dotplot.pdf <- gsub(".txt$", "_dotplot.pdf", output.file, perl=T)
		#output.file.dotplot.pdf <-  "/Users/hfang/Sites/XGR/XGRplus-site/app/examples/EAsnp_enrichment_dotplot.pdf"
		ggsave(output.file.dotplot.pdf, gp_dotplot, device=cairo_pdf, width=5, height=4)
		output.file.dotplot.png <- gsub(".txt$", "_dotplot.png", output.file, perl=T)
		ggsave(output.file.dotplot.png, gp_dotplot, type="cairo", width=5, height=4)
		
		# Forest plot
		if(1){
			message(sprintf("Drawing forest (%s) ...", as.character(Sys.time())), appendLF=TRUE)
			#zlim <- c(0, -log10(df_eTerm$adjp) %>% max() %>% ceiling())
			zlim <- c(0, -log10(df_eTerm$adjp) %>% quantile(0.95) %>% ceiling())
			gp_forest <- df_eTerm %>% mutate(name=str_c(name)) %>% oSEAforest(top=10, colormap="spectral.top", color.title=expression(-log[10]("FDR")), zlim=zlim, legend.direction=c("auto","horizontal","vertical")[3], sortBy=c("or","none")[1], size.title="Number\nof genes", wrap.width=50)
			output.file.forestplot.pdf <- gsub(".txt$", "_forest.pdf", output.file, perl=T)
			#output.file.forest.pdf <-  "/Users/hfang/Sites/XGR/XGRplus-site/app/examples/EAsnp_enrichment_forestplot.pdf"
			ggsave(output.file.forestplot.pdf, gp_forest, device=cairo_pdf, width=5, height=3.5)
			output.file.forestplot.png <- gsub(".txt$", "_forest.png", output.file, perl=T)
			ggsave(output.file.forestplot.png, gp_forest, type="cairo", width=5, height=3.5)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRplus-site/my_xgrplus/public
		## but outputs at public/tmp/eV2CG.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		if(1){
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		ind <- match(data$snp, names(GR.SNP))
		ls_rmd$data_input <- tibble(SNPs=data$snp[!is.na(ind)], Pvalue=data$pvalue[!is.na(ind)], GR.SNP[ind[!is.na(ind)]] %>% as.data.frame() %>% as_tibble() %>% transmute(Locus=str_c(seqnames,":",start,"-",end)))
		ls_rmd$xlsx_LG <- gsub("public/", "", output.file.LG, perl=T)
		ls_rmd$xlsx_LG_evidence <- gsub("public/", "", output.file.LG_evidence, perl=T)
		ls_rmd$num_LG <- nrow(df_LG)
		ls_rmd$min_overlap <- min.overlap
		ls_rmd$xlsx_enrichment <- gsub("public/", "", output.file.enrichment, perl=T)
		ls_rmd$pdf_dotplot <- gsub("public/", "", output.file.dotplot.pdf, perl=T)
		ls_rmd$png_dotplot <- gsub("public/", "", output.file.dotplot.png, perl=T)
		ls_rmd$pdf_forestplot <- gsub("public/", "", output.file.forestplot.pdf, perl=T)
		ls_rmd$png_forestplot <- gsub("public/", "", output.file.forestplot.png, perl=T)
		
		output_dir <- gsub("EAsnp.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_EAsnp.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[1], hightlight="default"), output_dir=output_dir)

		}
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(GenomicRanges)

# huawei
vec <- list.files(path='/root/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
#vec <- list.files(path='/Users/hfang/Sites/XGR/Fang/R', pattern='.r', full.names=T)
vec <- list.files(path='/Users/hfang/Sites/XGR/OpenXGR/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", crosslink=\"$crosslink\", significance.threshold=\"$significance_threshold\", obo=\"$obo\", FDR.cutoff=\"$FDR_cutoff\", min.overlap=\"$min_overlap\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- EAsnp: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_EAsnp.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_EAsnp.html";
			#public/tmp/RMD_eV2CG.html	
			print STDERR "RMD_EAsnp (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_EAsnp.html";
			#public/tmp/digit16_RMD_EAsnp.html
			print STDERR "RMD_EAsnp (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_EAsnp.html
				print STDERR "RMD_EAsnp (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "RElyser.html.ep"
sub RElyser {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $regionlist = $c->param('regionlist');
	my $crosslink = $c->param('crosslink') || 'proximity_10000';
  	my $obo = $c->param('obo') || 'GOMF'; # by default: GOMF
  	
	my $FDR_cutoff = $c->param('FDR_cutoff') || 0.05;
	my $min_overlap = $c->param('min_overlap') || 5;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($regionlist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;

		my $input_filename=$tmpFolder.'/'.'data.Regions.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'RElyser.Regions.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'RElyser.Regions.'.$digit16.'.r';

		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $regionlist)){
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_xgrm'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_xgrm";
		}elsif(-e '/var/www/html/bigdata_xgrm'){
			# huawei
			$placeholder="/var/www/html/bigdata_xgrm";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", crosslink="", obo="", FDR.cutoff="", min.overlap="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRm-site
		library(tidyverse)
		#placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_fdb"
		placeholder <- "/Users/hfang/Sites/SVN/github/bigdata_xgrm"

		input.file <- "~/Sites/XGR/XGRm-site/app/examples/eg_RElyser.txt"
		data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1)
		
		crosslink <- "proximity_50000"
		
		crosslink <- "RGBmm_EP_PMID32728240"
		
		FDR.cutoff <- 0.05
		min.overlap <- 3
		obo <- "GOMF"
	}
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1) %>% set_names("region")
	
	if(FDR.cutoff == "NULL"){
		FDR.cutoff <- 1
	}else{
		FDR.cutoff <- as.numeric(FDR.cutoff)
	}
	
	if(str_detect(crosslink, "proximity")){
		nearby.distance.max <- str_replace_all(crosslink, "proximity_", "") %>% as.numeric()
		crosslink <- "nearby"
		crosslink.customised  <- NULL
	}else{
		# replace RGB_PCHiC_ with RGB.PCHiC_
		crosslink <- str_replace_all(crosslink, "RGBmm_", "RGBmm.")
		crosslink.customised <- oRDS(crosslink, placeholder=placeholder)
	}
	
	xGene <- oGR2xGenes(data, format="chr:start-end", crosslink=crosslink, crosslink.customised=crosslink.customised, nearby.distance.max=nearby.distance.max, nearby.decay.kernel="constant", GR.Gene="UCSCmm_knownGene", placeholder=placeholder)
	
	min.overlap <- as.numeric(min.overlap)
	set <- oRDS(str_c("org.Mm.eg", obo), placeholder=placeholder)
	background <- set$info %>% pull(member) %>% unlist() %>% unique()
	eset <- xGene$xGene %>% pull(Gene) %>% oSEA(set, background, test="fisher", min.overlap=min.overlap)

	if(class(eset)=="eSET"){
		# *_LG.xlsx
		output.file.LG <- gsub(".txt$", "_LG.xlsx", output.file, perl=T)
		df_LG <- xGene$xGene
		df_LG %>% openxlsx::write.xlsx(output.file.LG)
		
		# *_LG_evidence.xlsx
		output.file.LG_evidence <- gsub(".txt$", "_LG_evidence.xlsx", output.file, perl=T)
		df_evidence <- xGene$Evidence %>% transmute(dGR,GR,Gene,Evidence=Context)
		df_evidence %>% openxlsx::write.xlsx(output.file.LG_evidence)

		# *_enrichment.txt
		df_eTerm <- eset %>% oSEAextract() %>% filter(adjp < FDR.cutoff)
		df_eTerm %>% write_delim(output.file, delim="\t")
		
		# *_enrichment.xlsx
		output.file.enrichment <- gsub(".txt$", ".xlsx", output.file, perl=T)
		df_eTerm %>% openxlsx::write.xlsx(output.file.enrichment)
		#df_eTerm %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/XGRm/app/examples/RElyser_enrichment.xlsx")
		
		# Dotplot
		message(sprintf("Drawing dotplot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		gp_dotplot <- df_eTerm %>% mutate(name=str_c(name)) %>% oSEAdotplot(FDR.cutoff=0.05, label.top=5, size.title="Number of genes", label.direction.y=c("left","right","none")[3], colors=c("skyblue","orangered"))
		output.file.dotplot.pdf <- gsub(".txt$", "_dotplot.pdf", output.file, perl=T)
		#output.file.dotplot.pdf <-  "/Users/hfang/Sites/XGR/XGRm-site/app/examples/RElyser_enrichment_dotplot.pdf"
		ggsave(output.file.dotplot.pdf, gp_dotplot, device=cairo_pdf, width=5, height=4)
		output.file.dotplot.png <- gsub(".txt$", "_dotplot.png", output.file, perl=T)
		ggsave(output.file.dotplot.png, gp_dotplot, type="cairo", width=5, height=4)
		
		# Forest plot
		if(1){
			message(sprintf("Drawing forest (%s) ...", as.character(Sys.time())), appendLF=TRUE)
			#zlim <- c(0, -log10(df_eTerm$adjp) %>% max() %>% ceiling())
			zlim <- c(0, -log10(df_eTerm$adjp) %>% quantile(0.95) %>% ceiling())
			gp_forest <- df_eTerm %>% mutate(name=str_c(name)) %>% oSEAforest(top=10, colormap="spectral.top", color.title=expression(-log[10]("FDR")), zlim=zlim, legend.direction=c("auto","horizontal","vertical")[3], sortBy=c("or","none")[1], size.title="Number\nof genes", wrap.width=50)		
			output.file.forestplot.pdf <- gsub(".txt$", "_forest.pdf", output.file, perl=T)
			#output.file.forest.pdf <-  "/Users/hfang/Sites/XGR/XGRm-site/app/examples/RElyser_enrichment_forestplot.pdf"
			ggsave(output.file.forestplot.pdf, gp_forest, device=cairo_pdf, width=5, height=3.5)
			output.file.forestplot.png <- gsub(".txt$", "_forest.png", output.file, perl=T)
			ggsave(output.file.forestplot.png, gp_forest, type="cairo", width=5, height=3.5)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRm-site/my_xgrm/public
		## but outputs at public/tmp/
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		if(1){
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- xGene$GR
		ls_rmd$xlsx_LG <- gsub("public/", "", output.file.LG, perl=T)
		ls_rmd$xlsx_LG_evidence <- gsub("public/", "", output.file.LG_evidence, perl=T)
		ls_rmd$num_LG <- nrow(df_LG)
		ls_rmd$min_overlap <- min.overlap
		ls_rmd$xlsx_enrichment <- gsub("public/", "", output.file.enrichment, perl=T)
		ls_rmd$pdf_dotplot <- gsub("public/", "", output.file.dotplot.pdf, perl=T)
		ls_rmd$png_dotplot <- gsub("public/", "", output.file.dotplot.png, perl=T)
		ls_rmd$pdf_forestplot <- gsub("public/", "", output.file.forestplot.pdf, perl=T)
		ls_rmd$png_forestplot <- gsub("public/", "", output.file.forestplot.png, perl=T)
		
		output_dir <- gsub("RElyser.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_RElyser.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[2], hightlight="default"), output_dir=output_dir)

		}
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(GenomicRanges)

# huawei
vec <- list.files(path='/root/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
vec <- list.files(path='/Users/hfang/Sites/XGR/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", crosslink=\"$crosslink\", obo=\"$obo\", FDR.cutoff=\"$FDR_cutoff\", min.overlap=\"$min_overlap\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- RElyser: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_RElyser.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_RElyser.html";
			#public/tmp/RMD_RElyser.html	
			print STDERR "RMD_RElyser (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_RElyser.html";
			#public/tmp/digit16_RMD_RElyser.html
			print STDERR "RMD_RElyser (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_RElyser.html
				print STDERR "RMD_RElyser (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "RSlyser.html.ep"
sub RSlyser {
  	my $c = shift;
	
	my $ip = $c->tx->remote_address;
	print STDERR "Remote IP address: $ip\n";
	
	my $host = $c->req->url->to_abs->host;
	my $port = $c->req->url->to_abs->port;
	my $host_port = "http://".$host.":".$port."/";
	print STDERR "Server available at ".$host_port."\n";
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $regionlist = $c->param('regionlist');
	my $crosslink = $c->param('crosslink') || 'proximity_10000';
  	my $network = $c->param('network') || 'STRING_high'; # by default: STRING_highest
	my $subnet_size = $c->param('subnet_size') || 30;
	my $subnet_sig = $c->param('subnet_sig') || 'yes';
  	
	my $significance_threshold = $c->param('significance_threshold') || 0.05;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($regionlist)){
		my $tmpFolder = $My_xgrm::Controller::Utils::tmpFolder; # public/tmp
		
		# 14 digits: year+month+day+hour+minute+second
		my $datestring = strftime "%Y%m%d%H%M%S", localtime;
		# 2 randomly generated digits
		my $rand_number = int rand 99;
		my $digit16 =$datestring.$rand_number."_".$ip;

		my $input_filename=$tmpFolder.'/'.'data.Regions.'.$digit16.'.txt';
		my $output_filename=$tmpFolder.'/'.'RSlyser.Regions.'.$digit16.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'RSlyser.Regions.'.$digit16.'.r';
	
		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $regionlist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			My_xgrm::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_xgrm'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_xgrm";
		}elsif(-e '/var/www/html/bigdata_xgrm'){
			# huawei
			$placeholder="/var/www/html/bigdata_xgrm";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", significance.threshold="", crosslink="", network="", subnet.size="", subnet.sig="", placeholder="", host.port="", ...){
	
	sT <- Sys.time()
	
	# for test
	if(0){
		#cd ~/Sites/XGR/XGRm-site
		#placeholder <- "~/Sites/SVN/github/bigdata_fdb"
		placeholder <- "~/Sites/SVN/github/bigdata_xgrm"
		
		input.file <- "~/Sites/XGR/XGRm-site/app/examples/eg_RSlyser.txt"
		data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
		
		format <- "data.frame"
		build.conversion <- c(NA,"hg38.to.hg19","hg18.to.hg19")
		
		significance.threshold=5e-2
		significance.threshold=NULL
		
		crosslink <- "proximity_50000"
		crosslink <- "RGBmm_EP_PMID32728240"
		
		network <- "STRINGmm_fin"
		network <- "KEGGmm_pin"
		subnet.size <- 30

	}
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2) %>% set_names("region","pvalue]")
	
	if(significance.threshold == "NULL"){
		significance.threshold <- NULL
	}else{
		significance.threshold <- as.numeric(significance.threshold)
	}
	
	if(str_detect(crosslink, "proximity")){
		nearby.distance.max <- str_replace_all(crosslink, "proximity_", "") %>% as.numeric()
		crosslink <- "nearby"
		crosslink.customised  <- NULL
	}else{
		# replace RGB_PCHiC_ with RGB.PCHiC_
		crosslink <- str_replace_all(crosslink, "RGBmm_", "RGBmm.")
		crosslink.customised <- oRDS(crosslink, placeholder=placeholder)
	}
	
	xGene <- oGR2xGeneScores(data, significance.threshold=significance.threshold, crosslink=crosslink, crosslink.customised=crosslink.customised, nearby.distance.max=nearby.distance.max, nearby.decay.kernel="constant", GR.Gene="UCSCmm_knownGene", placeholder=placeholder)
	
	subnet.size <- as.numeric(subnet.size)
	
	message(sprintf("Performing subnetwork analysis restricted to %d network genes (%s) ...", subnet.size, as.character(Sys.time())), appendLF=TRUE)
	#ig <- oRDS("KEGGmm_pin", placeholder=placeholder)
	ig <- oRDS(network, placeholder=placeholder)
	ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
	
	df_data <- tibble(name=xGene$xGene$Gene, pvalue=10^(-xGene$xGene$GScore)) %>% as.data.frame()
	subg <- oSubneterGenes(df_data, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder)

	if(vcount(subg)>0){
		# *_LG.xlsx
		output.file.LG <- gsub(".txt$", "_LG.xlsx", output.file, perl=T)
		df_LG <- xGene$xGene
		df_LG %>% openxlsx::write.xlsx(output.file.LG)
		
		# *_LG_evidence.xlsx
		output.file.LG_evidence <- gsub(".txt$", "_LG_evidence.xlsx", output.file, perl=T)
		df_evidence <- xGene$Evidence %>% transmute(dGR,GR,Gene,Evidence=Context)
		df_evidence %>% openxlsx::write.xlsx(output.file.LG_evidence)

		vec <- V(subg)$significance %>% as.numeric()
		vec[vec==0] <- min(vec[vec!=0])
		V(subg)$logP <- -log10(vec)
	
		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
		
		df_subg <- subg %>% oIG2TB("nodes") %>% transmute(Genes=name, Pvalue=as.numeric(significance), Description=description) %>% arrange(Pvalue)
		
		gp_rating <- oGGnetwork(g=subg, node.label="name", node.label.size=3, node.label.color="black", node.label.alpha=0.95, node.label.padding=0.5, node.label.arrow=0, node.label.force=0.4, node.shape=19, node.xcoord="xcoord", node.ycoord="ycoord", node.color="logP", node.color.title="Linked\ngene\nscores", colormap="spectral.top", zlim=c(0,10), node.size.range=5, title="", edge.color="steelblue4", edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05)
		
		
		# *_crosstalk.txt
		df_subg %>% write_delim(output.file, delim="\t")
		# *_crosstalk.xlsx
		output.file.crosstalk <- gsub(".txt$", "_crosstalk.xlsx", output.file, perl=T)
		df_subg %>% openxlsx::write.xlsx(output.file.crosstalk)

		# *_crosstalk.pdf *_crosstalk.png
		output.file.crosstalk.pdf <- gsub(".txt$", "_crosstalk.pdf", output.file, perl=T)
		ggsave(output.file.crosstalk.pdf, gp_rating, device=cairo_pdf, width=6, height=6)
		output.file.crosstalk.png <- gsub(".txt$", "_crosstalk.png", output.file, perl=T)
		ggsave(output.file.crosstalk.png, gp_rating, type="cairo", width=6, height=6)
		
		combinedP <- 1
		if(subnet.sig=="yes"){
			subg.sig <- oSubneterGenes(df_data, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder, test.permutation=T, num.permutation=10, respect=c("none","degree")[2], aggregateBy="fishers")
			combinedP <- signif(subg.sig$combinedP, digits=2)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/XGRm-site/my_xgrm/public
		## but outputs at public/tmp/
		######################################
		message(sprintf("RMD %s %f (%s) ...", subnet.sig, combinedP, as.character(Sys.time())), appendLF=TRUE)
		
		if(1){
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$host_port <- host.port
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- xGene$GR
		ls_rmd$xlsx_LG <- gsub("public/", "", output.file.LG, perl=T)
		ls_rmd$xlsx_LG_evidence <- gsub("public/", "", output.file.LG_evidence, perl=T)
		ls_rmd$num_LG <- nrow(df_LG)
		ls_rmd$vcount <- nrow(df_subg)
		ls_rmd$combinedP <- combinedP
		ls_rmd$xlsx_crosstalk <- gsub("public/", "", output.file.crosstalk, perl=T)
		ls_rmd$pdf_crosstalk <- gsub("public/", "", output.file.crosstalk.pdf, perl=T)
		ls_rmd$png_crosstalk <- gsub("public/", "", output.file.crosstalk.png, perl=T)
		
		output_dir <- gsub("RSlyser.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_RSlyser.Rmd", bookdown::html_document2(number_sections=F,theme=c("readable","united")[2], hightlight="default"), output_dir=output_dir)

		}
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(GenomicRanges)
library(igraph)

# huawei
vec <- list.files(path='/root/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
#vec <- list.files(path='/Users/hfang/Sites/XGR/Fang/R', pattern='.r', full.names=T)
vec <- list.files(path='/Users/hfang/Sites/XGR/XGRm/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", significance.threshold=\"$significance_threshold\", crosslink=\"$crosslink\", network=\"$network\", subnet.size=\"$subnet_size\", subnet.sig=\"$subnet_sig\", placeholder=\"$placeholder\", host.port=\"$host_port\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- RSlyser: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
My_xgrm::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			##########################
			### for RMD_RSlyser.html
			##########################
			$tmp_file=$tmpFolder."/"."RMD_RSlyser.html";
			#public/tmp/RMD_RSlyser.html	
			print STDERR "RMD_RSlyser (local & original) locates at $tmp_file\n";
			$ajax_rmd_html_file=$tmpFolder."/".$digit16."_RMD_RSlyser.html";
			#public/tmp/digit16_RMD_RSlyser.html
			print STDERR "RMD_RSlyser (local & new) locates at $ajax_rmd_html_file\n";
			if(-e $tmp_file){
				# do replacing
    			$command="mv $tmp_file $ajax_rmd_html_file";
				if(system($command)==1){
					print STDERR "Cannot execute: $command\n";
				}
				$ajax_rmd_html_file=~s/^public//g;
				#/tmp/digit16_RMD_RSlyser.html
				print STDERR "RMD_RSlyser (server) locates at $ajax_rmd_html_file\n";
			}
			
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


1;
