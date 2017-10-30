z = read.csv2("gwas-association-downloaded_2017-08-15-schizophrenia.tsv", sep="\t")

z$INITIAL.SAMPLE.SIZE = gsub(",", "", z$INITIAL.SAMPLE.SIZE)
z$SAMPLE.SIZE <- regmatches(z$INITIAL.SAMPLE.SIZE, gregexpr("[[:digit:]]+", z$INITIAL.SAMPLE.SIZE))
