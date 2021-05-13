# lastz-pipeline

Pairwise genome alignment pipeline

The following R script makes use of the CNEr package to pipeline the pairwise alignment of two genomes.

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
