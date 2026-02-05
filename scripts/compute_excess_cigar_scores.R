#!/usr/bin/env Rscript

# Debug: Print environment info immediately
cat("=== DEBUG START ===\n", file=stderr())
cat("TMPDIR env:", Sys.getenv("TMPDIR"), "\n", file=stderr())
cat("tempdir():", tempdir(), "\n", file=stderr())
cat("getwd():", getwd(), "\n", file=stderr())

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    cat("DEBUG: Running via snakemake\n", file=stderr())
    logfile <- snakemake@log[[1]]
    cat("DEBUG: logfile =", logfile, "\n", file=stderr())
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')
    cat("DEBUG: Sinks set up successfully\n")

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['config_yaml'],
            snakemake@params['single_cell'],
            snakemake@input['inttab'],
            snakemake@input['sccigars'],
            snakemake@input['bulkcigars'],
            snakemake@input['trainingdata'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@threads,
            snakemake@params['chroms']
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

cat("DEBUG: Parsing command line arguments...\n")
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 9) {
    stop("usage: compute_ab_ests_and_models.R config.yaml single.cell.ID integrated_table.tab.gz sccigars.tab.gz bulkcigars.tab.gz cigardata.tab.gz out.tab out.tab.gz n.cores [ chromosomes ]")
}

cat("DEBUG: Assigning variables from args...\n")
config.yaml <- args[1]
single.cell <- args[2]
int.tab <- args[3]
sccigars.path <- args[4]
bulkcigars.path <- args[5]
trainingcigars.path <- args[6]
out.tab <- args[7]
out.tab.gz <- args[8]
n.cores <- as.integer(args[9])

cat("DEBUG: config.yaml =", config.yaml, "\n")
cat("DEBUG: single.cell =", single.cell, "\n")
cat("DEBUG: int.tab =", int.tab, "\n")
cat("DEBUG: sccigars.path =", sccigars.path, "\n")
cat("DEBUG: bulkcigars.path =", bulkcigars.path, "\n")
cat("DEBUG: trainingcigars.path =", trainingcigars.path, "\n")
cat("DEBUG: out.tab =", out.tab, "\n")
cat("DEBUG: out.tab.gz =", out.tab.gz, "\n")
cat("DEBUG: n.cores =", n.cores, "\n")

chroms <- c()
if (length(args) > 9) {
    chroms <- args[-(1:9)]
    cat("DEBUG: chroms =", chroms, "\n")
}

cat("DEBUG: Checking if output files already exist...\n")
for (f in c(out.tab, out.tab.gz, paste0(out.tab.gz, '.tbi'))) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

cat("DEBUG: Loading library scan2...\n")
suppressMessages(library(scan2))
cat("DEBUG: Loading library future...\n")
suppressMessages(library(future))
cat("DEBUG: Loading library progressr...\n")
suppressMessages(library(progressr))

cat("DEBUG: Setting up plan with", n.cores, "cores...\n")
# Force sequential to debug the error
cat("DEBUG: FORCING SEQUENTIAL PLAN FOR DEBUGGING\n")
plan(sequential)

cat("DEBUG: Creating scan object with make.scan()...\n")
object <- make.scan(config.path=config.yaml, single.cell=single.cell)
cat("DEBUG: make.scan() completed\n")

# Restrict the analysis set to a set of chromosomes, if specified
if (length(chroms) > 0) {
    cat("DEBUG: Restricting analysis regions to chromosomes:", chroms, "\n")
    object@analysis.regions <- object@analysis.regions[seqnames(object@analysis.regions) %in% chroms,]
}
print(object@analysis.regions)

# Use fewer tiles when run on individual chromosomes to avoid unnecessary
# overhead from futures.
target.tile.number <- 300
if (length(chroms) > 0) {
    target.tile.number <- 100
}

cat('DEBUG: target.tile.number =', target.tile.number, '\n')

cat("DEBUG: Computing parallelization tiles...\n")
grs.for.para <- analysis.set.tiling.for.parallelization(object, total.tiles=target.tile.number)
cat("DEBUG: Parallelization tiles computed\n")

cat("DEBUG: Starting run.chunked.pipeline with what.to.compute='excess.cigar'...\n")
cat("DEBUG: int.tab =", int.tab, "\n")
cat("DEBUG: sccigars.path =", sccigars.path, "\n")
cat("DEBUG: bulkcigars.path =", bulkcigars.path, "\n")
cat("DEBUG: trainingcigars.path =", trainingcigars.path, "\n")

with_progress({
    handlers(handler_newline())
    object <- run.chunked.pipeline(object=object, int.tab=int.tab,
        sccigars=sccigars.path, bulkcigars=bulkcigars.path, trainingcigars=trainingcigars.path,
        grs.for.parallelization=grs.for.para,
        what.to.compute='excess.cigar',
        verbose=FALSE, report.mem=TRUE)
}, enable=TRUE)
cat("DEBUG: run.chunked.pipeline completed\n")

print('hi')
print(object@gatk)

cat("DEBUG: Writing results to", out.tab, "\n")
colnames(object@gatk)[1] <- '#chr'
data.table::fwrite(object@gatk[,.(`#chr`,pos,refnt,altnt,muttype,id.score,hs.score)],
    file=out.tab, sep='\t', quote=FALSE)
cat("DEBUG: fwrite completed\n")

cat("DEBUG: Running bgzip on", out.tab, "->", out.tab.gz, "\n")
Rsamtools::bgzip(out.tab, out.tab.gz)
cat("DEBUG: bgzip completed\n")

cat("DEBUG: Running indexTabix on", out.tab.gz, "\n")
Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')
cat("DEBUG: indexTabix completed\n")

cat("DEBUG: Script completed successfully\n")

if ('snakemake' %in% ls()) {
    sink()
}
