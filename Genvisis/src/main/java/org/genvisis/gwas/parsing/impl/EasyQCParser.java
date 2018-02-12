package org.genvisis.gwas.parsing.impl;

import java.io.IOException;
import org.genvisis.common.ext;
import org.genvisis.gwas.parsing.AliasedFileColumn;
import org.genvisis.gwas.parsing.DoubleWrapperColumn;
import org.genvisis.gwas.parsing.FileColumn;
import org.genvisis.gwas.parsing.FileParser;
import org.genvisis.gwas.parsing.FileParserFactory;
import org.genvisis.gwas.parsing.StandardFileColumns;

public class EasyQCParser {

  public void run(String inputFile, String outputFile) {
    String outFile = outputFile;
    if (outputFile == null) {
      outFile = ext.rootOf(inputFile, false) + "_out.txt";
    }

    FileColumn<String> snp = StandardFileColumns.snp("SNP");
    FileColumn<Byte> chr = StandardFileColumns.chr("CHR");
    FileColumn<Integer> pos = StandardFileColumns.pos("POS");
    FileColumn<String> strand = new AliasedFileColumn("STRAND", "strand");
    FileColumn<String> eff = StandardFileColumns.a1("EFFECT_ALLELE");
    FileColumn<String> oth = new AliasedFileColumn("OTHER_ALLELE", "base_allele");
    FileColumn<Integer> n = StandardFileColumns.n("N");
    FileColumn<Double> beta = new DoubleWrapperColumn(new AliasedFileColumn("BETA", "beta_int"));
    FileColumn<Double> se = new DoubleWrapperColumn(new AliasedFileColumn("SE", "se_int"));
    FileColumn<Double> p = new DoubleWrapperColumn(new AliasedFileColumn("PVAL", "chi_P_2df"));
    FileColumn<String> imp = new AliasedFileColumn("IMPUTATION", "Imputation_value");

    FileParser parser = FileParserFactory.setup(inputFile, snp, chr, pos, strand, eff, oth, n, beta,
                                                se, p, imp)
                                         .build();
    try {
      parser.parseToFile(outFile, " ");
    } catch (IOException e) {
      e.printStackTrace();
    }

  }

}
