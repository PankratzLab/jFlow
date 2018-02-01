package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import javax.swing.AbstractAction;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.FileChooser;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.CHROMOSOME;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Aliases;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.common.gui.JProgressBarListener;
import org.genvisis.common.gui.SimpleIndeterminateTask;
import org.genvisis.common.gui.Task;
import org.genvisis.gwas.parsing.AliasedFileColumn;
import org.genvisis.gwas.parsing.ColumnFilter;
import org.genvisis.gwas.parsing.DoubleWrapperColumn;
import org.genvisis.gwas.parsing.FileColumn;
import org.genvisis.gwas.parsing.FileParser;
import org.genvisis.gwas.parsing.FileParserFactory;
import org.genvisis.gwas.parsing.StandardFileColumns;
import org.genvisis.seq.manage.StrandOps.CONFIG;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.Allele;
import net.miginfocom.swing.MigLayout;

public class AFPlot {

	private JFrame frame;
	private AFPanel afPanel;
	private JPanel alleleInfoPanel;
	private Project proj;
	private Logger log;

	/**
	 * Loaded file
	 */
	private String file;
	/**
	 * Set of chrs to load
	 */
	private Set<Integer> chrs;
	/**
	 * Map of Marker name to marker objects. If null, attempts will be made to load location info from
	 * the data file.
	 */
	private Map<String, Marker> mkrMap;

	private Map<String, Marker> g1KMarkers;
	private Map<String, String> g1KChrPosLookup;
	private Map<String, Map<POPULATION, Double>> g1KData;

	private Map<String, String[]> obsAlleles;
	private Map<String, String[]> obsChrPosAlleles;
	private Map<String, Double> observeds;
	private Map<String, Double> observedsChrPos;
	private Map<CONFIG, Integer> alleleInfo;

	private volatile boolean forceRedraw = false;
	private volatile POPULATION selectedPop = POPULATION.ALL;
	private volatile boolean maskCenter = false;
	private volatile boolean chrPosLookup = false;
	private volatile boolean loading = false;
	private volatile boolean trim = true;
	private volatile boolean colorByConfig = true;

	/**
	 * Number of chromosomes constituting "all". Doesn't currently include X/Y/XY/MY.
	 */
	private static final int ALL_CHRS_COUNT = 22;
	private static final boolean DEFAULT_SHOW_ALLELE_INFO = false;
	private static final String TASK_CHANNEL = "LOAD_ALLELE_FREQS";

	private JCheckBoxMenuItem chrPosItem;

	/**
	 * 1KG allele freq populations
	 */
	enum POPULATION {
		ALL,
		EAS,
		EUR,
		AFR,
		AMR,
		SAS
	};

	private void loadReferencePanel(Set<String> dataSnps) {
		setG1KMarkers(new HashMap<>());
		setG1KData(new HashMap<>());
		setG1KChrPosLookup(new HashMap<>());
		Map<POPULATION, Double> dataMap;
		for (Integer chr : chrs) {
			FileColumn<String> snpCol = new AliasedFileColumn("SNP", "ID");
			FileColumn<Integer> chrCol = StandardFileColumns.chr("CHROM");
			FileColumn<Integer> posCol = StandardFileColumns.pos("POS");
			FileColumn<String> refCol = StandardFileColumns.a2("REF");
			FileColumn<String> altCol = StandardFileColumns.a1("ALT");
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
																					 .filter(new ColumnFilter() {
																						 @Override
																						 public List<FileColumn<?>> getFilterColumns() {
																							 return Lists.newArrayList(snpCol);
																						 }

																						 @Override
																						 public boolean filter(Map<FileColumn<?>, String> values) {
																							 String snpV = values.get(snpCol);
																							 String chrPos = values.get(chrCol)
																															 + ":"
																															 + values.get(posCol);
																							 boolean hasSnp = dataSnps.contains(snpV)
																																|| dataSnps.contains(chrPos);
																							 if (!hasSnp)
																								 return false;
																							 boolean goodAll = !values.get(refCol).contains(",")
																																 && !values.get(altCol)
																																					 .contains(",");
																							 if (!goodAll)
																								 return false;
																							 boolean goodChrPos = !ext.isMissingValue(values.get(chrCol))
																																		&& !ext.isMissingValue(values.get(posCol));
																							 return goodChrPos;
																						 }
																					 })
																					 .build();
			for (Map<FileColumn<?>, String> line : parser) {
				try {
					Marker m = new Marker(line.get(snpCol),
																new GenomicPosition(Byte.parseByte(line.get(chrCol)),
																										Integer.parseInt(line.get(posCol))),
																Allele.create(line.get(refCol), true),
																Allele.create(line.get(altCol), false));
					double[] afs = {
													Double.parseDouble(line.get(afAll)),
													Double.parseDouble(line.get(afEas)),
													Double.parseDouble(line.get(afEur)),
													Double.parseDouble(line.get(afAfr)),
													Double.parseDouble(line.get(afAmr)),
													Double.parseDouble(line.get(afSas))
					};
					getG1KMarkers().put(line.get(snpCol), m);
					getG1KChrPosLookup().put(line.get(chrCol) + ":" + line.get(posCol), line.get(snpCol));
					getG1KData().put(line.get(snpCol), dataMap = new HashMap<>());
					dataMap.put(POPULATION.ALL, afs[0]);
					dataMap.put(POPULATION.EAS, afs[1]);
					dataMap.put(POPULATION.EUR, afs[2]);
					dataMap.put(POPULATION.AFR, afs[3]);
					dataMap.put(POPULATION.AMR, afs[4]);
					dataMap.put(POPULATION.SAS, afs[5]);
				} catch (Exception e) {
				}
			}
			try {
				parser.close();
			} catch (IOException e) {
			}
		}
		dataMap = null;
	}

	private void setFile(String f) {
		this.file = f;
	}

	private void setChrs(Set<Integer> chrsToLoad) {
		chrs = new HashSet<Integer>();
		if (chrsToLoad == null) {
			for (int i = 1; i <= ALL_CHRS_COUNT; i++) {
				chrs.add(i);
			}
		} else {
			chrs.addAll(chrsToLoad);
		}
	}

	private void setMarkerMap(Map<String, Marker> map) {
		mkrMap = map == null ? null : Collections.unmodifiableMap(map);
	}

	private Set<String> loadDataSnps() throws IOException {
		FileColumn<String> snpCol = StandardFileColumns.snp("SNP");
		FileColumn<Integer> chrCol = StandardFileColumns.chr("CHR");

		FileParserFactory factory = FileParserFactory.setup(file, snpCol);
		if (mkrMap == null) {
			// if we don't have a marker map, try to load the chr column
			factory.optionalColumns(chrCol);
		}
		if (chrs.size() < ALL_CHRS_COUNT) {
			// filter the list of snps down to the selected chrs
			factory.filter(new ColumnFilter() {
				@Override
				public List<FileColumn<?>> getFilterColumns() {
					return mkrMap == null ? Lists.newArrayList(chrCol) : new ArrayList<>();
				}

				@Override
				public boolean filter(Map<FileColumn<?>, String> values) {
					boolean chr = false;
					if (mkrMap != null && mkrMap.containsKey(values.get(snpCol))) {
						chr = chrs.contains((int) mkrMap.get(values.get(snpCol)).getChr());
					} else if (values.containsKey(chrCol)) {
						chr = chrs.contains(Integer.parseInt(values.get(chrCol)));
					}
					return chr;
				}
			});
		}
		long lines = -1;
		if (trim) {
			lines = Files.countLines(file, 1);
			if (lines > 1000000) {
				// if a huge file, trim
				int lineMod = (int) (lines / 1000000);
				factory.filter(new ColumnFilter() {

					@Override
					public List<FileColumn<?>> getFilterColumns() {
						return new ArrayList<>();
					}

					long count = 0;

					@Override
					public boolean filter(Map<FileColumn<?>, String> values) {
						return count++ % lineMod == 0;
					}
				});
			}
		}

		Set<String> dataSnps = factory.build().load(snpCol).keySet();
		long rsCount = countRSIds(dataSnps);
		if (rsCount == 0) {
			log.report("No RS-IDs found in data - switching to CHR:POS mode.");
			setChrPosLookup(true);
			if (chrPosItem != null) {
				chrPosItem.setSelected(true);
			}
			rsCount = dataSnps.size();
		}
		if (trim) {
			if (lines != dataSnps.size()) {
				log.report("Trimmed data snps to " + dataSnps.size() + " of "
									 + lines);
			}
		}
		return dataSnps;
	}

	private void loadDataFreqs() throws IOException {
		FileColumn<String> snpCol = StandardFileColumns.snp("SNP");
		FileColumn<Double> mafCol = new DoubleWrapperColumn(new AliasedFileColumn("MAF",
																																							Aliases.ALLELE_FREQS));
		FileColumn<String> a1 = StandardFileColumns.a1("REF");
		FileColumn<String> a2 = StandardFileColumns.a2("ALT");
		FileColumn<Integer> chr = StandardFileColumns.chr("CHR");
		FileColumn<Integer> pos = StandardFileColumns.pos("POS");
		FileParserFactory factory = FileParserFactory.setup(file, snpCol, mafCol);
		if (mkrMap == null) {
			factory.optionalColumns(chr, pos, a1, a2);
		}
		factory.filter(new ColumnFilter() {
			@Override
			public List<FileColumn<?>> getFilterColumns() {
				return Lists.newArrayList(snpCol);
			}

			@Override
			public boolean filter(Map<FileColumn<?>, String> values) {
				boolean key = getG1KMarkers().containsKey(values.get(snpCol));
				if (!key) {
					if (values.containsKey(chr) && values.containsKey(pos)) {
						key = getG1KChrPosLookup().containsKey(values.get(chr)
																									 + ":"
																									 + values.get(pos));
					}
				}
				boolean miss = !ext.isMissingValue(values.get(mafCol));
				return key && miss;
			}
		});
		FileParser parser = factory.build();

		for (Map<FileColumn<?>, String> line : parser) {
			String key = line.get(snpCol);
			double maf = Double.parseDouble(line.get(mafCol));
			getObservedData().put(key, maf);
			if (line.containsKey(a1) && line.containsKey(a2)) {
				getObservedAlleles().put(key, new String[] {line.get(a1), line.get(a2)});
			}
			if (mkrMap != null) {
				Marker mkr = mkrMap.get(key);
				if (mkr != null) {
					String[] alleles = new String[] {mkr.getRef().getBaseString(),
																					 mkr.getAlt().getBaseString()};
					getObservedAlleles().put(key, alleles);
					key = mkr.getGenomicPosition().getChr() + ":"
								+ mkr.getGenomicPosition().getPosition();
					getObservedDataChrPos().put(key, maf);
					getObservedChrPosAlleles().put(key, alleles);
				}
			}
			if (line.containsKey(chr) && line.containsKey(pos)) {
				key = line.get(chr) + ":" + line.get(pos);
				getObservedDataChrPos().put(key, maf);
				if (line.containsKey(a1) && line.containsKey(a2)) {
					getObservedChrPosAlleles().put(key, new String[] {line.get(a1), line.get(a2)});
				}
			}
		}
		parser.close();
	}

	/**
	 * Load data from a file. Assumes the file has a header.<br />
	 * If this {@link AFPlot} was constructed without a Project, assumes the file has chr/pos info.
	 * <br />
	 * Loads allele info if present in file.
	 * 
	 * @param file Path to file
	 * @param chrsToLoad Set of chromosome numbers to load
	 */
	public void loadFromFile(String file, Set<Integer> chrsToLoad) {
		setLoading(true);
		setFile(file);
		setChrs(chrsToLoad);

		JProgressBarListener listener = new JProgressBarListener(TASK_CHANNEL);
		setProgressBar(listener.getBar());

		Task<String, String> loadTask = new SimpleIndeterminateTask(TASK_CHANNEL) {

			@Override
			protected String doInBackground() throws Exception {
				String msg;
				msg = "Loading list of data snps...";
				log.reportTime(msg);
				afPanel.setNullMessage(msg);
				afPanel.paintAgain();

				try {
					Set<String> dataSnps = loadDataSnps();

					msg = "Loading reference panel allele freq data...";
					afPanel.setNullMessage(msg);
					afPanel.paintAgain();
					log.reportTime(msg);

					loadReferencePanel(dataSnps);

					msg = ("Loading data snp allele freqs...");
					afPanel.setNullMessage(msg);
					afPanel.paintAgain();
					log.reportTime(msg);

					loadDataFreqs();

					msg = "Drawing plot - "
								+ (isChrPosLookup() ? getObservedDataChrPos().size() : getObservedData().size())
								+ " matched...";
					afPanel.setNullMessage(msg);
					log.reportTime(msg);
				} catch (IOException e) {
					setError(e);
				}
				return null;
			}
		};

		loadTask.execute();
		try {
			loadTask.get();
			setLoading(false);
		} catch (InterruptedException | ExecutionException e) {
			setError(e);
		}

		removeProgressBar(listener.getBar());
		afPanel.paintAgain();
		// iterate through 1000G markers, pull existing out of markerNmMap and MAF data
		// - check alleles (if present) and respond appropriately to differences
	}

	private void setError(Exception e) {
		String msg = e.getMessage();
		afPanel.setNullMessage(msg);
		log.reportException(e);
	}

	private void setProgressBar(JProgressBar bar) {
		frame.add(bar, "cell 0 1, growx");
		frame.revalidate();
		frame.repaint();
	}

	private void removeProgressBar(JProgressBar bar) {
		frame.remove(bar);
		frame.revalidate();
		frame.repaint();
	}

	private long countRSIds(Set<String> dataSnps) {
		return dataSnps.stream().filter(s -> s.startsWith("rs")).count();
	}

	public void screenshot(String file) {
		if (this.afPanel.getSize().getWidth() <= 350 || this.afPanel.getSize().getHeight() <= 100) {
			this.afPanel.setSize(800, 800);
		}
		this.afPanel.screenCapture(file);
		log.reportTime("Screenshot to " + file);
	}

	private void createMenuBar() {
		JMenuBar bar = new JMenuBar();
		JMenu options = new JMenu("Options");
		bar.add(options);

		/**
		 * Load file menu item
		 */
		JMenuItem loadFileItem = new JMenuItem();
		loadFileItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				FileChooser fc = new FileChooser(frame, "", false, false, "Load File", log);
				if (!fc.isSelected())
					return;
				final String file = fc.getSelectedFile().getAbsolutePath();
				new Thread(() -> {
					loadFromFile(file, null);
				}).start();
				setForceRedraw(true);
				afPanel.paintAgain();
			}
		});
		loadFileItem.setText("Load File");
		loadFileItem.setMnemonic('F');
		options.add(loadFileItem);

		/**
		 * Create screenshot
		 */
		JMenuItem screenshotItem = new JMenuItem();
		screenshotItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				FileChooser fc = new FileChooser(frame, "", false, false, "Save Screenshot", log);
				if (!fc.isSelected())
					return;
				final String file = fc.getSelectedFile().getAbsolutePath();
				new Thread(() -> {
					screenshot(file);
				}).start();
			}
		});
		screenshotItem.setText("Screenshot");
		screenshotItem.setMnemonic('S');
		options.add(screenshotItem);

		options.addSeparator();

		/**
		 * Mask center line of plot
		 */
		JCheckBoxMenuItem maskCenterItem = new JCheckBoxMenuItem();
		maskCenterItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				AFPlot.this.setMaskCenter(!AFPlot.this.isMaskCenter());
				setForceRedraw(true);
				afPanel.paintAgain();
			}
		});
		maskCenterItem.setText("Mask Center");
		maskCenterItem.setSelected(AFPlot.this.isMaskCenter());
		maskCenterItem.setMnemonic('M');
		options.add(maskCenterItem);

		/**
		 * Use Chr:Pos lookup for markers; this will change automatically if, when loading a file, no
		 * rsIDs are found.
		 */
		chrPosItem = new JCheckBoxMenuItem();
		chrPosItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				AFPlot.this.setChrPosLookup(!AFPlot.this.isChrPosLookup());
				setForceRedraw(true);
				afPanel.paintAgain();
			}
		});
		chrPosItem.setText("Match Chr/Pos for Snps");
		chrPosItem.setSelected(AFPlot.this.isChrPosLookup());
		chrPosItem.setMnemonic('C');
		options.add(chrPosItem);

		/**
		 * Show panel with allele config info
		 */
		JCheckBoxMenuItem showAlleleItem = new JCheckBoxMenuItem();
		showAlleleItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				AFPlot.this.showHideAlleleInfo();
				afPanel.paintAgain();
			}
		});
		showAlleleItem.setText("Show Allele Stats");
		showAlleleItem.setSelected(DEFAULT_SHOW_ALLELE_INFO);
		showAlleleItem.setMnemonic('A');
		options.add(showAlleleItem);

		/**
		 * Trim data if large (negated option)
		 */
		JCheckBoxMenuItem dontTrimItem = new JCheckBoxMenuItem();
		dontTrimItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				AFPlot.this.trim = !AFPlot.this.trim;
				afPanel.paintAgain();
			}
		});
		dontTrimItem.setText("Don't Trim Data (slow)");
		dontTrimItem.setSelected(!trim);
		dontTrimItem.setMnemonic('T');
		options.add(dontTrimItem);

		/**
		 * Color points based on allele config info (and info is present)
		 */
		JCheckBoxMenuItem colorConfigItem = new JCheckBoxMenuItem();
		colorConfigItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				AFPlot.this.setColorByConfig(!AFPlot.this.isColorByConfig());
				afPanel.paintAgain();
			}
		});
		colorConfigItem.setText("Color by Allele Config");
		colorConfigItem.setSelected(isColorByConfig());
		colorConfigItem.setMnemonic('C');
		options.add(colorConfigItem);

		options.addSeparator();

		/**
		 * Select population to compare freqs against
		 */
		JMenu popMenu = new JMenu("Population");
		ButtonGroup bg = new ButtonGroup();
		for (POPULATION p : POPULATION.values()) {
			JRadioButtonMenuItem pI = new JRadioButtonMenuItem();
			pI.setAction(new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent e) {
					setSelectedPop(p);
					setForceRedraw(true);
					afPanel.paintAgain();
				}
			});
			pI.setText(ext.capitalizeFirst(p.name()));
			bg.add(pI);
			pI.setSelected(p == getSelectedPop());
			popMenu.add(pI);
		}
		options.add(popMenu);

		frame.setJMenuBar(bar);
	}

	void updateAlleleInfoPanel() {
		if (alleleInfoPanel == null)
			return;
		StringBuilder sb = new StringBuilder();
		sb.append("<html><body><table border=1>");

		for (CONFIG c : CONFIG.values()) {
			sb.append("<tr><td>").append(c.getDescription()).append("</td><td>")
				.append(getAlleleInfo().get(c) == null ? 0 : getAlleleInfo().get(c))
				.append("</td></tr>");
		}

		sb.append("</table></body></html>");
		SwingUtilities.invokeLater(() -> {
			alleleInfoPanel.removeAll();
			JLabel html = new JLabel(sb.toString());
			alleleInfoPanel.add(html);
		});
	}

	protected void showHideAlleleInfo() {
		if (alleleInfoPanel == null)
			return;
		updateAlleleInfoPanel();
		alleleInfoPanel.setVisible(!alleleInfoPanel.isVisible());
	}

	/**
	 * Create an AlleleFrequency plot. If not headless, will create a {@link JFrame}.
	 * 
	 * @param proj {@link Project}; can be null.
	 */
	public AFPlot(Project proj) {
		afPanel = new AFPanel(this);
		if (!GraphicsEnvironment.isHeadless()) {
			frame = new JFrame();
			frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
			frame.setSize(800, 700);
			frame.setLayout(new MigLayout("hidemode 3", "0px[grow]0px", "0px[grow, fill]0px[]0px"));
			frame.add(afPanel, "cell 0 0, grow");
			frame.setTitle("Genvisis - AFPlot"
										 + (proj == null ? "" : " - " + proj.PROJECT_NAME.getValue()));
			createMenuBar();

			alleleInfoPanel = new JPanel(new MigLayout("debug", "[center, grow]", "[center, grow]"));
			alleleInfoPanel.setBackground(Color.WHITE);
			alleleInfoPanel.setVisible(false);
			frame.add(alleleInfoPanel, "west, center");
			frame.revalidate();
		}

		log = proj == null ? new Logger() : proj.getLog();
		setMarkerMap(proj == null ? null : proj.getMarkerSet().getMarkerNameMap());

		setObservedData(new HashMap<>());
		setObsAlleles(new HashMap<>());
		setObservedDataChrPos(new HashMap<>());
		setObsChrPosAlleles(new HashMap<>());
		setAlleleInfo(new HashMap<>());
	}

	public void setVisible(boolean vis) {
		frame.setVisible(vis);
	}

	/**
	 * Wait for {@link AFPlot#isLoading()} to return false.
	 */
	public void waitForData() {
		while (isLoading()) {
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
			}
		}
	}

	public static void main(String[] args) {
		AFPlot plot = new AFPlot(null);
		plot.setVisible(true);

		// Set<Integer> chrs = null;// Sets.newHashSet(4);
		// plot.loadFromFile("F:\\BCX2\\results\\fixed\\CARDIA_BAS_EA_141217_BC.txt.gz", chrs);
		// plot.waitForData();
	}

	private void setObsChrPosAlleles(Map<String, String[]> obsChrPosAlleles) {
		this.obsChrPosAlleles = obsChrPosAlleles;
	}

	private void setObsAlleles(Map<String, String[]> obsAlleles) {
		this.obsAlleles = obsAlleles;
	}

	private void setObservedDataChrPos(Map<String, Double> observedsChrPos) {
		this.observedsChrPos = observedsChrPos;
	}

	private void setObservedData(Map<String, Double> observeds) {
		this.observeds = observeds;
	}

	private void setG1KChrPosLookup(Map<String, String> g1kChrPosLookup) {
		g1KChrPosLookup = g1kChrPosLookup;
	}

	private void setG1KData(Map<String, Map<POPULATION, Double>> g1kData) {
		g1KData = g1kData;
	}

	private void setLoading(boolean loading) {
		this.loading = loading;
	}

	private void setMaskCenter(boolean maskCenter) {
		this.maskCenter = maskCenter;
	}

	private void setChrPosLookup(boolean chrPosLookup) {
		this.chrPosLookup = chrPosLookup;
	}

	private void setSelectedPop(POPULATION selectedPop) {
		this.selectedPop = selectedPop;
	}

	private void setColorByConfig(boolean colorByConfig) {
		this.colorByConfig = colorByConfig;
	}

	private void setG1KMarkers(Map<String, Marker> g1kMarkers) {
		g1KMarkers = g1kMarkers;
	}

	private void setAlleleInfo(Map<CONFIG, Integer> alleleInfo) {
		this.alleleInfo = alleleInfo;
	}

	/**
	 * @return Map of Chr:Pos to String[2] alleles
	 */
	public Map<String, String[]> getObservedChrPosAlleles() {
		return obsChrPosAlleles;
	}

	/**
	 * @return Map of marker name to String[2] alleles
	 */
	public Map<String, String[]> getObservedAlleles() {
		return obsAlleles;
	}

	/**
	 * @return Map of Chr:Pos to AF value
	 */
	public Map<String, Double> getObservedDataChrPos() {
		return observedsChrPos;
	}

	/**
	 * @return Map of marker name to AF value
	 */
	public Map<String, Double> getObservedData() {
		return observeds;
	}

	/**
	 * @return Map of 1KG Chr:Pos to Markername
	 */
	public Map<String, String> getG1KChrPosLookup() {
		return g1KChrPosLookup;
	}

	/**
	 * @return Map of 1KG Markername to map of population to population allele frequency
	 */
	public Map<String, Map<POPULATION, Double>> getG1KData() {
		return g1KData;
	}

	/**
	 * @return Should force a redraw
	 */
	public boolean isForceRedraw() {
		return forceRedraw;
	}

	/**
	 * @param forceRedraw boolean force a redraw
	 */
	public void setForceRedraw(boolean forceRedraw) {
		this.forceRedraw = forceRedraw;
	}

	/**
	 * @return is loading a file?
	 */
	public boolean isLoading() {
		return loading;
	}

	/**
	 * @return should mask the center of the plot?
	 */
	public boolean isMaskCenter() {
		return maskCenter;
	}

	/**
	 * @return is using a chr:pos lookup rather than marker name?
	 */
	public boolean isChrPosLookup() {
		return chrPosLookup;
	}

	public POPULATION getSelectedPop() {
		return selectedPop;
	}

	/**
	 * @return Map of 1KG String to {@link Marker}
	 */
	public Map<String, Marker> getG1KMarkers() {
		return g1KMarkers;
	}

	/**
	 * @return Map of allele configuration info
	 */
	public Map<CONFIG, Integer> getAlleleInfo() {
		return alleleInfo;
	}

	/**
	 * @return should color points by allele {@link CONFIG}?
	 */
	public boolean isColorByConfig() {
		return colorByConfig;
	}

}

