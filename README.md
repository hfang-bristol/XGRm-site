# [Source code for the XGRm web server](https://github.com/hfang-bristol/XGRm-site)

## @ Overview

> The [XGRm](http://www.xgrmod.com) web server offers real-time enrichment and subnetwork analyses for a user-input list of `genes`, `SNPs`, or `genomic regions`, through leveraging ontologies, networks, and many others (such as e/pQTL, promoter capture Hi-C, and enhancer-gene maps).

> The web server provides users with a more integrated and user-friendly experience, featuring: [ENRICHMENT ANALYSER (GENES) - EAG](http://www.xgrmod.com/XGRm/EAgene) identifying enriched ontology terms from input gene list; [ENRICHMENT ANALYSER (REGIONS) - EAR](http://www.xgrmod.com/XGRm/EAregion) identifying enriched ontology terms for genes linked from input genomic region list; [SUBNETWORK ANALYSER (GENES) - SAG](http://www.xgrmod.com/XGRm/SAgene) identifying a gene subnetwork based on input gene-level summary data; and [SUBNETWORK ANALYSER (REGIONS) - SAR](http://www.xgrmod.com/XGRm/SAregion) identifying a gene subnetwork based on genes linked from input genomic region-level summary data.

> To learn how to use the XGRm web server, a user manual has been made available [here](http://www.xgrmod.com/XGRmodbooklet/index.html) with step-by-step instructions.

## @ Development

> The XGRm web server was developed using [Mojolicious](https://www.mojolicious.org) and [Bootstrap](https://getbootstrap.com), supporting a mobile-first and responsive webserver across all major platform browsers.

> The folder `my_xgrm` has a tree-like directory structure with three levels:
```ruby
my_xgrm
├── lib
│   └── My_xgrm
│       └── Controller
├── public
│   ├── XGRmbooklet
│   │   ├── index_files
│   │   └── libs
│   ├── app
│   │   ├── ajex
│   │   ├── css
│   │   ├── examples
│   │   ├── img
│   │   └── js
│   ├── dep
│   │   ├── HighCharts
│   │   ├── Select2
│   │   ├── bootstrap
│   │   ├── bootstrapselect
│   │   ├── bootstraptoggle
│   │   ├── dataTables
│   │   ├── fontawesome
│   │   ├── jcloud
│   │   ├── jqcloud
│   │   ├── jquery
│   │   ├── tabber
│   │   └── typeahead
│   └── tmp
├── script
├── t
└── templates
    └── layouts
```


## @ Installation

Assume you have a `ROOT (sudo)` privilege on `Ubuntu`

### 1. Install Mojolicious and other perl modules

```ruby
sudo su
# here enter your password

curl -L cpanmin.us | perl - Mojolicious
perl -e "use Mojolicious::Plugin::PODRenderer"
perl -MCPAN -e "install Mojolicious::Plugin::PODRenderer"
perl -MCPAN -e "install DBI"
perl -MCPAN -e "install Mojo::DOM"
perl -MCPAN -e "install Mojo::Base"
perl -MCPAN -e "install LWP::Simple"
perl -MCPAN -e "install JSON::Parse"
perl -MCPAN -e "install local::lib"
perl -MCPAN -Mlocal::lib -e "install JSON::Parse"
```

### 2. Install R

```ruby
sudo su
# here enter your password

# install R
wget http://www.stats.bris.ac.uk/R/src/base/R-4/R-4.3.1.tar.gz
tar xvfz R-4.3.1.tar.gz
cd ~/R-4.3.1
./configure
make
make check
make install
R # start R
```

### 3. Install pandoc

```ruby
sudo su
# here enter your password

# install pandoc
wget https://github.com/jgm/pandoc/releases/download/3.1/pandoc-3.1-linux-amd64.tar.gz
tar xvzf pandoc-3.1-linux-amd64.tar.gz --strip-components 1 -C /usr/local/

# use pandoc to render R markdown
R
rmarkdown::pandoc_available()
Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
rmarkdown::render(YOUR_RMD_FILE, bookdown::html_document2(number_sections=F, theme="readable", hightlight="default"))
```


## @ Deployment

Assume you place `my_xgrm` under your `home` directory

```ruby
cd ~/my_xgrm
systemctl stop apache2.service
morbo -l 'http://*:80/' script/my_xgrm
```

## @ Contact

> For any bug reports, please drop [email](mailto:fh12355@rjh.com.cn).


