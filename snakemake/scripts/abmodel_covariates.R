#!/usr/bin/env Rscript

# Debug: Print environment info to stderr immediately
cat("=== DEBUG START ===\n", file=stderr())
cat("TMPDIR env:", Sys.getenv("TMPDIR"), "\n", file=stderr())
cat("tempdir():", tempdir(), "\n", file=stderr())
cat("R_HOME:", Sys.getenv("R_HOME"), "\n", file=stderr())
cat("getwd():", getwd(), "\n", file=stderr())

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    cat("DEBUG: Running via snakemake\n", file=stderr())
    logfile <- snakemake@log[[1]]
    cat("DEBUG: logfile =", logfile, "\n", file=stderr())
    cat("DEBUG: Opening log file connection...\n", file=stderr())
    con <- file(logfile, 'w')
    cat("DEBUG: Connection opened, setting up sink for output...\n", file=stderr())
    sink(con, type='output')
    cat("DEBUG: Output sink set, setting up sink for message...\n", file=stderr())
    sink(con, type='message')
    cat("DEBUG: Both sinks set up successfully\n")

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['integrated_table'],
            snakemake@input['abfits'],
            snakemake@input['bedgz'],
            snakemake@params['sc_sample'],
            snakemake@params['genome'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

cat("DEBUG: Parsing command line arguments...\n")
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
    stop("usage: spatial_sens_abmodel.R integrated_table.tab.gz abfits.rda regions.bed.gz singlecell_id genome_string out.tab out.tab.gz n_cores")
}

cat("DEBUG: Assigning variables from args...\n")
integrated.table.path <- args[1]
abfits.path <- args[2]
bed.path <- args[3]
sc.sample <- args[4]
genome.string <- args[5]
out.tab <- args[6]
out.tab.gz <- args[7]
n.cores <- as.integer(args[8])

cat("DEBUG: integrated.table.path =", integrated.table.path, "\n")
cat("DEBUG: abfits.path =", abfits.path, "\n")
cat("DEBUG: bed.path =", bed.path, "\n")
cat("DEBUG: sc.sample =", sc.sample, "\n")
cat("DEBUG: genome.string =", genome.string, "\n")
cat("DEBUG: out.tab =", out.tab, "\n")
cat("DEBUG: out.tab.gz =", out.tab.gz, "\n")
cat("DEBUG: n.cores =", n.cores, "\n")

cat("DEBUG: Checking if output files already exist...\n")
for (f in c(out.tab, out.tab.gz))
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))

cat("DEBUG: Loading library scan2...\n")
suppressMessages(library(scan2))
cat("DEBUG: Loading library Rsamtools...\n")
suppressMessages(library(Rsamtools))
cat("DEBUG: Loading library future...\n")
suppressMessages(library(future))
cat("DEBUG: Loading library progressr...\n")
suppressMessages(library(progressr))
cat("DEBUG: Setting up multicore plan with", n.cores, "workers...\n")
plan(multicore, workers=n.cores)

cat("DEBUG: Loading ab.fits from", abfits.path, "...\n")
ab.fits <- get(load(abfits.path))
cat("DEBUG: ab.fits loaded successfully\n")

cat("DEBUG: Reading sens.regions from", bed.path, "...\n")
sens.regions <- data.table::fread(cmd = paste("zcat", bed.path))
colnames(sens.regions) <- c('chr', 'start', 'end')
sens.regions$end <- sens.regions$end - 1
cat("DEBUG: sens.regions has", nrow(sens.regions), "rows\n")

cat("DEBUG: Creating GRanges object...\n")
gr.sens.regions <- GenomicRanges::GRanges(seqnames=sens.regions$chr,
    ranges=IRanges::IRanges(start=sens.regions$start, end=sens.regions$end),
    seqinfo=genome.string.to.seqinfo.object(genome.string))
cat("DEBUG: GRanges object created\n")

cat("DEBUG: Creating parallelization tiles...\n")
# Don't want that many tiles here; each one causes two extra tabix reads.
grs.for.para <- analysis.set.tiling.for.parallelization.helper(
    regions=GenomicRanges::reduce(gr.sens.regions),
    total.tiles=10)
cat("DEBUG: Parallelization tiles created\n")

cat("DEBUG: Starting compute.spatial.sensitivity.abmodel...\n")
with_progress({
    handlers(handler_newline())
    abmodel <- compute.spatial.sensitivity.abmodel(
        single.cell.id=sc.sample,
        ab.fits=ab.fits,
        integrated.table.path=integrated.table.path,
        grs.for.sens=gr.sens.regions,
        grs.for.parallelization=grs.for.para,
        genome.string=genome.string)
}, enable=TRUE)
cat("DEBUG: compute.spatial.sensitivity.abmodel completed\n")

cat("DEBUG: Writing results to", out.tab, "\n")
colnames(abmodel)[1] <- "#chr"
data.table::fwrite(abmodel, file=out.tab, sep="\t", col.names=TRUE)
cat("DEBUG: fwrite completed\n")

cat("DEBUG: Running bgzip on", out.tab, "->", out.tab.gz, "\n")
Rsamtools::bgzip(out.tab, out.tab.gz)
cat("DEBUG: bgzip completed\n")

cat("DEBUG: Running indexTabix on", out.tab.gz, "\n")
Rsamtools::indexTabix(file=out.tab.gz, format='bed', comment='#')
cat("DEBUG: indexTabix completed\n")

cat("DEBUG: Script completed successfully\n")

if ('snakemake' %in% ls()) {
    sink()
}
