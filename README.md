# lastz-pipeline

Pairwise genome alignment pipeline

The following R script makes use of the CNEr package to pipeline the pairwise alignment of two genomes. The pipeline performs pairwise alignment of all specified chromosomes or scaffolds of target and query genomes using [lastz](https://www.bx.psu.edu/~rsharris/lastz/). The remainder of the pipeline uses some binaries from UCSC [kent tools](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/) for filtering and chaining alignments. The number of pairwise combinations grows quadratically, and the runtime becomes unruly with too many scaffolds or chromosomes. When working with draft genomes, often the most reasonable approach is to limit the input to the top 30-40 longest scaffolds. The same goes for reference genomes, for example limit the chrom.sizes file to contain only cannonical chromosomes, no `chrUn*`. 

Two primate genomes and 50 threads can take 3-5 hours, so make use of as many CPUs as your can find.

# Setup

Install all dependencies in a fresh conda environment using the following command, assuming you have a working installation of conda. To install miniconda on an linux system follow these [instructions](https://docs.conda.io/en/latest/miniconda.html).

```
conda env create -n lastz-test -f environment.yaml
```

# Usage

```
Rscript lastz-pipe.R -h
Options:
	-t CHARACTER, --target=CHARACTER
		Target genome 2bit file (absolute path).

	-q CHARACTER, --query=CHARACTER
		Query genome 2bit file (absolute path).

	-T CHARACTER, --targetSizes=CHARACTER
		Target chrom.sizes file specifies chroms to include. Uses known chroms from twoBit file if not given.

	-Q CHARACTER, --querySizes=CHARACTER
		Query chrom.sizes file specifies chroms to include. Uses known chros from twoBit file if not given.

	-d CHARACTER, --outdir=CHARACTER
		Output directory (absolute path).

	-p NUMERIC, --threads=NUMERIC
		Number of threads to use for lastz

	-h, --help
		Show this help message and exit
```


The output file of primary importance has the form {target}.{query}.net.axt. This file can be ingested by our synteny breakpoint detection tool [axtToSyn](https://github.com/carbonelab/axtToSyn).

To create a useful visualization of this pairwise genome alignment, you can convert it to `maf` format using the kent utility `axtToMaf` as follows:

```
axtToMaf in.axt tSizes qSizes out.maf
```

Then create pariwise alignment dotplot using the tool [last-dotplot](https://github.com/lpryszcz/last/blob/master/scripts/last-dotplot). 

```
last-dotplot [options] maf-or-tab-alignments dotplot.png
```

If the dotplot is too busy with smaller noisy alignments, you can filter the lower scoreing alignments out of the maf file with `mafFiler`. I have found 10000 can be a noise-reducing threshold, but can vary based on genome assembly quality and completeness.

```
mafFiler -minScore=10000 in.maf out.filt.maf
```
