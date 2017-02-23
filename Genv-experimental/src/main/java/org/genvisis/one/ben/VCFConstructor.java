package org.genvisis.one.ben;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeader.HEADER_FIELDS;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.seq.manage.ReferenceGenome;

public class VCFConstructor {
	
	Logger log;
	String inputFile;
	String dosageFile;
	String outputFile;
	
	String[] markers;
	int[][] locations;
	List<List<String>> alleles;
	String[] ids;
	
	private void readInputFile() {
		// TODO check inputFile is set and exists and is readable
		markers = HashVec.loadFileToStringArray(inputFile, false, false, new int[] {0}, true, false, "\t");
		SnpMarkerSet markerSet = new SnpMarkerSet(markers);
		markerSet.parseSNPlocations(log);
		locations = markerSet.getChrAndPositionsAsInts();
		parseAlleles(markerSet.getMarkerNames());
	}
	
	private void readDosageFile() {
		// TODO check dosageFile is set and exists and is readable
		String[][] matr = HashVec.loadFileToStringMatrix(dosageFile, false, null, false);
		ids = ArrayUtils.subArray(Matrix.extractColumn(matr, 1), 1);
		String[] snps = ArrayUtils.subArray(matr[0], 2);
		
	}
	
	private void parseAlleles(String[] markerNames) {
		alleles = new ArrayList<List<String>>();
		for (int i = 0; i < markerNames.length; i++) {
			if (markerNames[i].startsWith("rs")) {
				log.reportError("Cannot parse alleles from RS marker");
				alleles.add(new ArrayList<String>());
			} else {
				String[] pts = markerNames[i].split(":");
				ArrayList<String> a = new ArrayList<String>();
				a.add(pts[2]);
				a.add(pts[3]);
				alleles.add(a);
			}
		}
	}
	
	public void build() {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputFile);
		builder.clearOptions();
		builder.setOption(Options.INDEX_ON_THE_FLY);
		VCFHeader vcfHeader = new VCFHeader();
		vcfHeader.addMetaDataLine(new VCFFormatHeaderLine("FORMAT=<ID=GT,Number=1,Type=String,Description=\"A blank genotype\">", VCFHeaderVersion.VCF4_2));
		
		SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(Resources.genome(GENOME_BUILD.HG19, log).getFASTA().getAbsolute(), log)
				.getIndexedFastaSequenceFile().getSequenceDictionary();

		builder.setReferenceDictionary(samSequenceDictionary);
		vcfHeader.setSequenceDictionary(samSequenceDictionary);
		VariantContextWriter writer = builder.build();
		writer.writeHeader(vcfHeader);

		for (int i = 0; i < markers.length; i++) {/// your data
			VariantContextBuilder builderVc = new VariantContextBuilder();
			builderVc.chr("chr" + locations[i][0]);
			int len = alleles.get(i).get(0).replaceAll("\\*","").length() - 1;
			List<Allele> l = new ArrayList<Allele>();
			for (int a = 0; a < alleles.get(i).size(); a++) {
				l.add(Allele.create(alleles.get(i).get(a), a == 0));
			}
			builderVc.alleles(l);
			builderVc.start(locations[i][1]);
			builderVc.stop(locations[i][1] + len);
			builderVc.id(markers[i]);
			Collection<Genotype> col = new ArrayList<Genotype>();
			for (String id : ids) {
				col.add(GenotypeBuilder.create(id, l));
			}
			builderVc.genotypes(col);

			writer.add(builderVc.make());
		}
		writer.close();
	}
	
	public static void main(String[] args) {
		VCFConstructor c = new VCFConstructor();
		c.inputFile = "F:/temp/variantviewer/inputData.txt";
		c.outputFile = "F:/temp/variantviewer/data.vcf";
		c.dosageFile = "F:/temp/variantviewer/FRZ.db.xln";
		c.log = new Logger();
		c.readInputFile();
		c.readDosageFile();
		c.build();
	}
	
}
