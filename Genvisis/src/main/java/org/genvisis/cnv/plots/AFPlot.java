package org.genvisis.cnv.plots;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.CHROMOSOME;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.Logger;
import org.genvisis.gwas.parsing.AliasedFileColumn;
import org.genvisis.gwas.parsing.ColumnFilter;
import org.genvisis.gwas.parsing.DoubleWrapperColumn;
import org.genvisis.gwas.parsing.FileColumn;
import org.genvisis.gwas.parsing.FileParser;
import org.genvisis.gwas.parsing.FileParserFactory;
import org.genvisis.gwas.parsing.StandardFileColumns;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.Allele;

public class AFPlot {

	AbstractPanel AFPanel = new AbstractPanel() {

		@Override
		public void highlightPoints() {

		}

		@Override
		public void generatePoints() {

		}

		@Override
		public void assignAxisLabels() {
			displayXAxis = displayYAxis = true;
			xAxisLabel = "Observed";
			yAxisLabel = "Expected";
		}
	};

	Project proj;
	Logger log;
	Set<Integer> chrs;

	Set<String> sharedMarkers;

	Map<String, Marker> g1KMarkers;
	Map<String, Map<POPULATION, Double>> g1KData;

	Map<String, Double> observeds;

	enum POPULATION {
		ALL,
		EAS,
		EUR,
		AFR,
		AMR,
		SAS
	};

	private void load1000G() {
		// TODO progress bar / monitor?
		// TODO use Task api?
		g1KMarkers = new HashMap<>();
		g1KData = new HashMap<>();
		Map<POPULATION, Double> dataMap;
		for (Integer chr : chrs) {
			FileColumn<String> snpCol = StandardFileColumns.snp("SNP");
			FileColumn<Integer> chrCol = StandardFileColumns.chr("CHROM");
			FileColumn<Integer> posCol = StandardFileColumns.pos("POS");
			FileColumn<String> refCol = StandardFileColumns.a1("REF");
			FileColumn<String> altCol = StandardFileColumns.a2("ALT");
			FileColumn<Double> afAll = new DoubleWrapperColumn(new AliasedFileColumn("AF", "AF"));
			FileColumn<Double> afEas = new DoubleWrapperColumn(new AliasedFileColumn("EAS", "EAS_AF"));
			FileColumn<Double> afEur = new DoubleWrapperColumn(new AliasedFileColumn("EUR", "EUR_AF"));
			FileColumn<Double> afAfr = new DoubleWrapperColumn(new AliasedFileColumn("AFR", "AFR_AF"));
			FileColumn<Double> afAmr = new DoubleWrapperColumn(new AliasedFileColumn("AMR", "AMR_AF"));
			FileColumn<Double> afSas = new DoubleWrapperColumn(new AliasedFileColumn("SAS", "SAS_AF"));

			String G1KFile = Resources.genome(proj == null ? GENOME_BUILD.HG19
																										 : proj.GENOME_BUILD_VERSION.getValue(),
																				log)
																.chr(CHROMOSOME.valueOf("C" + Integer.toString(chr)))
																.getG1Kphase3v5AlleleFreq().get();
			FileParser parser = FileParserFactory.setup(G1KFile, snpCol, chrCol, posCol, refCol, altCol,
																									afAll, afEas, afEur, afAfr, afAmr, afSas)
																					 .build();
			for (Map<FileColumn<?>, String> line : parser) {
				g1KMarkers.put(line.get(snpCol), new Marker(line.get(snpCol),
																										new GenomicPosition(Byte.parseByte(line.get(chrCol)),
																																				Integer.parseInt(line.get(posCol))),
																										Allele.create(line.get(refCol), true),
																										Allele.create(line.get(altCol), false)));
				g1KData.put(line.get(snpCol), dataMap = new HashMap<>());
				dataMap.put(POPULATION.ALL, Double.parseDouble(line.get(afAll)));
				dataMap.put(POPULATION.EAS, Double.parseDouble(line.get(afEas)));
				dataMap.put(POPULATION.EUR, Double.parseDouble(line.get(afEur)));
				dataMap.put(POPULATION.AFR, Double.parseDouble(line.get(afAfr)));
				dataMap.put(POPULATION.AMR, Double.parseDouble(line.get(afAmr)));
				dataMap.put(POPULATION.SAS, Double.parseDouble(line.get(afSas)));
			}
			try {
				parser.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		dataMap = null;
	}

	private void setChrs(Set<Integer> chrsToLoad) {
		chrs = new HashSet<Integer>();
		if (chrsToLoad == null) {
			for (int i = 1; i < 23; i++) { // TODO add X/Y/MT?
				chrs.add(i);
			}
		} else {
			chrs.addAll(chrsToLoad);
		}
	}

	private void loadFromMarkerMetrics(Project proj, Set<Integer> chrsToLoad) {
		loadFromFile(proj.MARKER_METRICS_FILENAME.getValue(), chrsToLoad);

		// Map<String, Marker> markerNmMap = proj.getMarkerSet().getMarkerNameMap();
		// could iterate through alleles and check if correct / flipped ... or just create the plot and
		// let the user find out if they're flipped
		// - maybe include a menu option to check for / flip alleles
	}

	private void loadFromFile(String file, Set<Integer> chrsToLoad) {
		setChrs(chrsToLoad);
		load1000G();

		FileColumn<String> snpCol = StandardFileColumns.snp("SNP");
		FileColumn<Double> mafCol = new DoubleWrapperColumn(new AliasedFileColumn("MAF", "MAF"));
		// TODO add Load-If-Present functionality to FileParser
		FileParser parser = FileParserFactory.setup(file, snpCol, mafCol)
																				 .filter(new ColumnFilter() {
																					 @Override
																					 public List<FileColumn<?>> getFilterColumns() {
																						 return Lists.newArrayList(snpCol);
																					 }

																					 @Override
																					 public boolean filter(Map<FileColumn<?>, String> values) {
																						 return g1KMarkers.containsKey(values.get(snpCol));
																					 }
																				 }).build();
		for (Map<FileColumn<?>, String> line : parser) {
			observeds.put(line.get(snpCol), Double.parseDouble(line.get(mafCol)));
		}
		// iterate through 1000G markers, pull existing out of markerNmMap and MAF data
		// - check alleles (if present) and respond appropriately to differences
	}

}
