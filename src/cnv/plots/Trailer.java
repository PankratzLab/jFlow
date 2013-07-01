// to run from within this class, the argument -notJar must be passed to it
package cnv.plots;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import common.*;
import cnv.filesys.*;
import cnv.gui.SingleClick;
import cnv.gui.ClickListener;
import cnv.manage.Transforms;
import cnv.var.CNVariant;
import cnv.var.IndiPheno;
import cnv.var.SampleData;
import filesys.*;

public class Trailer extends JFrame implements ActionListener, ClickListener, MouseListener, MouseMotionListener, MouseWheelListener {
	public static final long serialVersionUID = 1L;

//	public static final String DEFAULT_LOCATION = "chr6:161,739,001-163,119,245"; // parkin
//	public static final String DEFAULT_LOCATION = "chr5:151497267-151499003"; // hit
//	public static final String DEFAULT_LOCATION = "chr5:151,175,382-151,203,885"; // 1 gene
	public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32
//	public static final String DEFAULT_LOCATION = "chr12"; // USP32

	public static final String DEFAULT_SAMPLE = null;
//	public static final String DEFAULT_SAMPLE = "ND12014";
	
	// public static final String DEFAULT_LOCATION = "chr10:96,040,187-97,106,178";
	// public static final String DEFAULT_LOCATION = "chr23";

	// public static final String[] DEFAULT_CNV_FILES = {"data/penncnv.cnv", "data/quantisnp_windows.cnv", "data/quantisnp_linux.cnv", "data/dnacopy10.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"data/conf.cnv", "data/allMarkers.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"data/conf_100kb_5SNP_10.0.cnv", "data/allMarkers_100kb_5SNP_10.0.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"data/conf_100bp_5SNP_10.0.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"conf_100bp_5SNP_10.0.cnv"};
//	public static final String[] DEFAULT_CNV_FILES = {"data/penncnv_1SNP.cnv", "data/quantisnp_1SNP.cnv"};
//	public static final String DEFAULT_GENE_TRACK = GeneSet.DIRECTORY+GeneSet.REFSEQ_TRACK;
//	public static final String DEFAULT_GENE_TRACK = GeneSet.REFSEQ_TRACK;
	public static final boolean SHOW_MIDLINE = true;

	public static final double MOUSE_WHEEL_MULTIPLIER = 0.5;
	public static final int WIDTH_BUFFER = 25;
	public static final int HEIGHT_BUFFER = 10;
//	public static final int DYNAMIC_HEIGHT_LIMIT = 2000;
	public static final int DYNAMIC_HEIGHT_LIMIT = 0;
	public static final int DOUBLE_CLICK_INTERVAL = 500;
	public static final int SIZE = 4;
	public static final int MIN_BUFFER = 10000;
	// public static final String DEFAULT_SAMPLE_DATA = "data/SampleData.dat";
	public static final int DEFAULT_STARTX = 20;
	public static final int DEFAULT_STARTY = 20;
	public static final int DEFAULT_WIDTH = 1000;
	public static final int DEFAULT_HEIGHT = 720;

	// private static final String ALT_UP = "ALT UP";
	// private static final String ALT_DOWN = "ALT DOWN";
	// private static final String ALT_LEFT = "ALT LEFT";
	// private static final String ALT_RIGHT = "ALT RIGHT";

	private static final String FIRST_CHR = "First chr";
	private static final String PREVIOUS_CHR = "Previous chr";
	private static final String NEXT_CHR = "Next chr";
	private static final String LAST_CHR = "Last chr";
	private static final String FIRST_REGION = "First region";
	private static final String PREVIOUS_REGION = "Previous region";
	private static final String NEXT_REGION = "Next region";
	private static final String LAST_REGION = "Last region";

	private JComboBox<String> sampleList;
	private String[] samplesPresent;
	private JTextField navigationField;
	private JButton firstChr, previousChr, nextChr, lastChr, previousRegion, nextRegion;
	private Project proj;
	private String sample;
	private IndiPheno indiPheno;
	private boolean jar;
	private int numMarkers;
	private MarkerSet markerSet;
	private long fingerprint;
	private int[] positions;
	private boolean[] dropped;
	private int[][] chrBoundaries;
	private float[] lrrs, lrrValues;
	private float[] bafs;
	private byte[] genotypes;
	private byte chr;
	private int start, startMarker;
	private int stop, stopMarker;
//	private float lrrMin, lrrMax;
	private boolean inDrag;
	private int startX;
	private String[] cnvFilenames;
	private String[] cnvLabels;
	private CNVariant[][] cnvs;
	private String[] regionsList;
	private int regionsListIndex;
	private String[][] regions;
	private int regionIndex;
	private JPanel lrrPanel;
	private SingleClick leftClick;
	private SingleClick rightClick;
	private GeneTrack track;
	private SampleData sampleData;
	private JLabel commentLabel;
	private int transformation_type;
	private boolean transformSeparatelyByChromosome;
	
	// private Color[] colorScheme = {new Color(33, 31, 53), // dark dark
	// new Color(23, 58, 172), // dark blue
	// new Color(201, 30, 10), // deep red
	// new Color(140, 20, 180), // deep purple
	// new Color(33, 87, 0), // dark green
	// new Color(55, 129, 252), // light blue
	// new Color(217, 109, 194), // pink
	// new Color(94, 88, 214), // light purple
	// new Color(189, 243, 61), // light green
	// };

	private Color[][] colorScheme = {
			{new Color(23, 58, 172), new Color(55, 129, 252)}, // dark/light blue
			{new Color(140, 20, 180), new Color(94, 88, 214)}, // deep/light purple
			{new Color(33, 87, 0), new Color(189, 243, 61)}, // dark green
			{new Color(201, 30, 10), new Color(217, 109, 194)}, // deep red/pink
			{new Color(33, 31, 53), new Color(255, 255, 255)} // dark dark/ light light
	};
	
	public Trailer(Project proj, String selectedSample, String[] filenames, String location) {
		this(proj, selectedSample, filenames, location, DEFAULT_STARTX, DEFAULT_STARTX, DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}

	// TODO Trailer should have a createAndShowGUI, same as all the other plots, as opposed to being its own frame 
	public Trailer(Project proj, String selectedSample, String[] filenames, String location, int startX, int startY, int width, int height) {
		super("CNVis - Trailer");
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		System.out.println("startX: "+startX+"\t startY: "+startY+"\t width: "+width+"\t height: "+height);
		System.out.println(Array.toStr(filenames, "; "));

		long time;

		this.proj = proj;
		jar = proj.getJarStatus();
		cnvFilenames = filenames;
		
		time = new Date().getTime();

		chr = (byte)Positions.parseUCSClocation(location)[0]; 
//		parseLocation(location, false);

		loadMarkers();
		generateComponents();
		
		sample = selectedSample==null?samplesPresent[0]:selectedSample;
		sampleData = proj.getSampleData(2, cnvFilenames);
		if (sampleData.failedToLoad()) {
			return;
		}
		cnvLabels = sampleData.getCnvClasses();
		System.err.println("Error - ");
		System.err.println("Error - there are "+cnvLabels.length+" cnvLabels: "+Array.toStr(cnvFilenames));
		System.err.println("Error - ");
//		System.exit(1);
//		loadCNVfiles(proj, cnvFilenames);

		regionsList = proj.getFilenames(Project.REGION_LIST_FILENAMES);
		regionsListIndex = 0;
		if (regionsList.length > 0) {
			if (Files.exists(regionsList[regionsListIndex], jar)) {
				loadRegions();
			} else {
				System.err.println("Error - couldn't find '"+regionsList[regionsListIndex]+"' in data directory; populating with CNVs of current subject");
			}
		}
		regionIndex = -1;

        time = new Date().getTime();
//		track = GeneTrack.load(proj.getDir(Project.DATA_DIRECTORY)+GeneSet.REFSEQ_TRACK, jar);
        if (new File(proj.getFilename(Project.GENETRACK_FILENAME, false, false)).exists()) {
        	track = GeneTrack.load(proj.getFilename(Project.GENETRACK_FILENAME), jar);
        } else if (new File(GeneSet.DIRECTORY+GeneSet.REFSEQ_TRACK).exists()) {
            track = GeneTrack.load(GeneSet.DIRECTORY+GeneSet.REFSEQ_TRACK, jar);
        } else if (new File(GeneSet.REFSEQ_TRACK).exists()) {
        	track = GeneTrack.load(GeneSet.REFSEQ_TRACK, jar);
        } else {
//			JOptionPane.showMessageDialog(this, "Gene track is not installed. Gene boundaries will not be displayed.", "FYI", JOptionPane.INFORMATION_MESSAGE);
        	track = null;
        }
		System.out.println("Loaded track in "+ext.getTimeElapsed(time));
		
		
		updateSample(sample);
		
		System.out.println("All in "+ext.getTimeElapsed(time));

		parseLocation(location);
//		setBounds(20, 20, 1000, 720);
		setBounds(startX, startY, width, height);
		setVisible(true);
	}
	
	public void generateComponents() {
		
		JPanel dataPanel = new JPanel();
		dataPanel.setLayout(new GridLayout(3, 1));

		lrrPanel = new JPanel() {
			public static final long serialVersionUID = 2L;

			public void paintComponent(Graphics g) {
				float min, max;

				if (lrrValues != null) {
					min = Array.min(lrrValues);
					max = Array.max(lrrValues);

					// System.out.println("Displaying "+(stopMarker-startMarker)+"
					// markers");
					if (stopMarker-startMarker<DYNAMIC_HEIGHT_LIMIT) {
						min = Float.POSITIVE_INFINITY;
						max = Float.NEGATIVE_INFINITY;
						for (int i = startMarker; i<=stopMarker; i++) {
							if (lrrValues[i]<min) {
								min = lrrValues[i];
							}
							if (lrrValues[i]>max) {
								max = lrrValues[i];
							}
						}
						min = (float)Math.floor(min);
						max = (float)Math.ceil(max);
					}
	
					// min = -1.5;
					// max = 1.5;

					switch (transformation_type) {
					case 0:
						min = -3;
						max = 3;
						break;
					case 1:
						min = 0;
						max = 1;
						break;
					case 2:
						min = -8;
						max = 8;
						break;
					case 3:
						min = -12;
						max = 12;
						break;
						
					}
					
					g.setFont(new Font("Arial", 0, 20));
					g.drawString("Log R Ratio", WIDTH_BUFFER, 20);
					
					if (SHOW_MIDLINE) {
						g.setColor(Color.LIGHT_GRAY);
						g.drawLine(WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER, getWidth()-WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
					}
	
					g.setFont(new Font("Arial", 0, 12));
					for (int i = startMarker; i<=stopMarker; i++) {
						// if (genotypes[i] == 1) {
						if (bafs[i]>0.2&&bafs[i]<0.8) {
							g.setColor(Color.RED);
							// colorScheme[2]
						} else {
							g.setColor(Color.BLACK);
						}
						if (!Float.isNaN(lrrValues[i])) {
							if (dropped[i]) {
//								g.drawString("X", getX(positions[i]), getHeight()-(int)((double)(lrrValues[i]-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
							} else if (lrrValues[i] < min){
								g.drawString("v", getX(positions[i]), getHeight()-(int)((double)(min-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
							} else if (lrrValues[i] > max){
								g.drawString("^", getX(positions[i]), getHeight()-(int)((double)(max-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
							} else {
								g.fillOval(getX(positions[i]), getHeight()-(int)((double)(lrrValues[i]-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER, SIZE, SIZE);
							}
						}
					}
				}
			}

			public int getX(int pos) {
				return (int)((double)(pos-start)/(double)(stop-start)*(double)(getWidth()-2*WIDTH_BUFFER))+WIDTH_BUFFER;
			}
		};
		lrrPanel.addMouseListener(this);
		lrrPanel.addMouseMotionListener(this);
		lrrPanel.addMouseWheelListener(this);
		dataPanel.add(lrrPanel);

		JPanel panel = new JPanel() {
			public static final long serialVersionUID = 8L;

			public void paintComponent(Graphics g) {
				GeneData[] genes;
				int[][] exons;
				Vector<Segment> v = new Vector<Segment>();
				Segment[] segs;
				int width, begin, end, source;
				Segment currentView;
				String text;
				
//				g.drawRect(0, 0, this.getWidth()-1, this.getHeight()-1);
				
				
				if (track == null) {
					text = "Gene track is not installed";
					width = g.getFontMetrics(g.getFont()).stringWidth(text);
					g.drawString(text, this.getWidth()/2-width/2, 10);
				} else {
					if (stop-start > 10000000) {
						g.drawString("Zoom in to see genes", 10, 10);
					} else {
						genes = track.getBetween(chr, start, stop, 30);
//						System.out.println(ext.getUCSCformat(new int[] {chr, start, stop}));
						g.setColor(Color.BLACK);
						for (int i = 0; i<genes.length; i++) {
							begin = getX(genes[i].getStart());
							end = getX(genes[i].getStop());
							g.drawRoundRect(begin, 0*15, end-begin, 10, 2, 2);
							v.add(new Segment(begin, end));
							exons = genes[i].getExonBoundaries();
							for (int j = 0; j<exons.length; j++) {
								begin = getX(exons[j][0]);
								end = getX(exons[j][1]);
								if (j==0 || j==exons.length-1) {
									g.fillRoundRect(begin, 0*15, end-begin+1, 10, 2, 2);
								} else {
									g.fillRect(begin, 0*15, end-begin+1, 10);
								}
								
							}
//							System.out.println(genes[i].getGeneName()+"\t"+genes[i].getStart()+"\t"+genes[i].getStop());
	                    }
//						System.out.println();
						Segment.mergeOverlapsAndSort(v);
						segs = Segment.toArray(v);
						g.setFont(new Font("Arial", 0, 14));
						
						for (int i = 0; i<genes.length; i++) {
							begin = getX(genes[i].getStart());
							width = g.getFontMetrics(g.getFont()).stringWidth(genes[i].getGeneName());
							if (!Segment.overlapsAny(new Segment(begin-width-5, begin-1), segs)) {
								g.drawString(genes[i].getGeneName(), begin-width-3, 0*15+10);
							}
						}
					}
				}
				currentView = new Segment(chr, start, stop);
				int firstBegin;
				for (int i = 0; cnvs != null && i<cnvs.length; i++) {
					source = i;
					firstBegin = Integer.MAX_VALUE;
					for (int j = 0; j<cnvs[i].length; j++) {
						if (cnvs[i][j].overlaps(currentView)) {
							begin = getX(cnvs[i][j].getStart());
							if (begin < firstBegin) {
								firstBegin = begin;
							}
							end = getX(cnvs[i][j].getStop());
							g.setColor(colorScheme[source][cnvs[i][j].getCN()<2?0:1]);
							g.fillRoundRect(begin, (source+2)*15, end-begin+1, 10, 2, 2);
//							g.drawString(ext.rootOf(cnvLabels[source]), begin+2, (source+1)*15+10);
						}
					}
					if (firstBegin != Integer.MAX_VALUE) {
						g.drawString(ext.rootOf(cnvLabels[source]), firstBegin-Grafik.getTextWidth(cnvLabels[source], g)-3, (source+2)*15+10);
					}
				}
			}

			public int getX(int pos) {
				return (int)((double)(pos-start)/(double)(stop-start)*(double)(getWidth()-2*WIDTH_BUFFER))+WIDTH_BUFFER;
			}
		};
		panel.setMaximumSize(new Dimension(getWidth(), 20));
		dataPanel.add(panel);

		panel = new JPanel() {
			public static final long serialVersionUID = 7L;

			public void paintComponent(Graphics g) {
				g.setFont(new Font("Arial", 0, 20));
				g.drawString("B Allele Frequency", WIDTH_BUFFER, 20);

				g.setFont(new Font("Arial", 0, 10));
				if (lrrs != null) {
					for (int i = startMarker; i<=stopMarker; i++) {
						if (!Float.isNaN(lrrs[i])) {
							if (dropped[i]) {
								g.drawString("X", getX(positions[i]), getHeight()-(int)(bafs[i]*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER+5);
							} else if (genotypes != null && genotypes[i]==-1) {
								g.drawString("+", getX(positions[i]), getHeight()-(int)(bafs[i]*(double)(getHeight()-4*HEIGHT_BUFFER))-HEIGHT_BUFFER/2);
							} else {
								g.fillOval(getX(positions[i]), getHeight()-(int)(bafs[i]*(double)(getHeight()-4*HEIGHT_BUFFER))-HEIGHT_BUFFER, SIZE, SIZE);
							}
						}
					}
				}
			}

			public int getX(int pos) {
				return (int)((double)(pos-start)/(double)(stop-start)*(double)(getWidth()-2*WIDTH_BUFFER))+WIDTH_BUFFER;
			}
		};
		dataPanel.add(panel);

		// getContentPane().setBackground(Color.WHITE);
		getContentPane().add(dataPanel, BorderLayout.CENTER);

		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new BoxLayout(descrPanel, BoxLayout.Y_AXIS));
//		descrPanel.setLayout(new GridLayout(2, 1));

		// JLabel sampleName = new JLabel(sample, JLabel.CENTER);
		// sampleName.setFont(new Font("Arial", 0, 20));
		// descrPanel.add(sampleName);

		JPanel sampPanel = new JPanel();
		previousRegion = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previousRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previousRegion.addActionListener(this);
		previousRegion.setActionCommand(PREVIOUS_REGION);
		previousRegion.setPreferredSize(new Dimension(25, 25));
		sampPanel.add(previousRegion);

		sampleList = new JComboBox<String>();
		sampleList.setFont(new Font("Arial", 0, 20));
		createSampleList();
		
		DefaultListCellRenderer dlcr = new DefaultListCellRenderer();
	    dlcr.setHorizontalAlignment(DefaultListCellRenderer.CENTER);
	    sampleList.setRenderer(dlcr);
		sampleList.setBorder(BorderFactory.createEtchedBorder());
		sampleList.setEditable(false);
		sampleList.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("unchecked")
				JComboBox<String> jcb = (JComboBox<String>)e.getSource();
				int index = jcb.getSelectedIndex();
//				System.out.println("Selected Index = "+index);
				if (index == samplesPresent.length-1) {
					createSampleList();
				} else if (sample != samplesPresent[index]) {
					updateSample(samplesPresent[index]);
				}
			}
		});
		sampPanel.add(sampleList);
		
		nextRegion = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		nextRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		nextRegion.addActionListener(this);
		nextRegion.setActionCommand(NEXT_REGION);
		nextRegion.setPreferredSize(new Dimension(25, 25));
		sampPanel.add(nextRegion);

		descrPanel.add(sampPanel);
		

		descrPanel.add(Box.createHorizontalGlue());
		commentLabel = new JLabel("", JLabel.CENTER);
		commentLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		commentLabel.setFont(new Font("Arial", 0, 14));
		descrPanel.add(commentLabel);
		descrPanel.add(Box.createHorizontalGlue());

		JPanel navigationPanel = new JPanel();
		firstChr = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
		firstChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
		firstChr.addActionListener(this);
		firstChr.setActionCommand(FIRST_CHR);
		firstChr.setPreferredSize(new Dimension(20, 20));
		previousChr = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previousChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previousChr.addActionListener(this);
		previousChr.setActionCommand(PREVIOUS_CHR);
		previousChr.setPreferredSize(new Dimension(20, 20));
		navigationField = new JTextField("", 20);
		navigationField.setHorizontalAlignment(JTextField.CENTER);
		navigationField.setFont(new Font("Arial", 0, 14));
		navigationField.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	        	JTextField navField;
	        	
	        	navField = (JTextField)e.getSource();
	        	navField.setText(navField.getText().trim());
				parseLocation(navField.getText());
	        }
		});

//		navigationField.addFocusListener(new FocusListener() {
//			public void focusGained(FocusEvent focusevent) {}
//
//			public void focusLost(FocusEvent fe) {
//				parseLocation(((JTextField)fe.getSource()).getText());
//			}
//		});

		nextChr = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		nextChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		nextChr.addActionListener(this);
		nextChr.setActionCommand(NEXT_CHR);
		nextChr.setPreferredSize(new Dimension(20, 20));
		lastChr = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
		lastChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		lastChr.addActionListener(this);
		lastChr.setActionCommand(LAST_CHR);
		lastChr.setPreferredSize(new Dimension(20, 20));
		navigationPanel.add(firstChr);
		navigationPanel.add(previousChr);
		navigationPanel.add(navigationField);
		navigationPanel.add(nextChr);
		navigationPanel.add(lastChr);
		descrPanel.add(navigationPanel);
		
		JPanel transformationPanel = new JPanel();
		JLabel label = new JLabel("Log R Ratio transformation: ");
		label.setFont(new Font("Arial", 0, 14));
		transformationPanel.add(label);
		
		ButtonGroup typeRadio = new ButtonGroup();
		JRadioButton[] transformationRadioButtons = new JRadioButton[Transforms.TRANFORMATIONS.length];
		ItemListener typeListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					for (int i = 0; i<Transforms.TRANFORMATIONS.length; i++) {
						if (jrb.getText().equals(Transforms.TRANFORMATIONS[i])) {
							transformation_type = i;
							System.out.println("Transformation type: "+transformation_type);
							if (transformation_type > 0) {
								lrrValues = Transforms.transform(lrrs, transformation_type, transformSeparatelyByChromosome, markerSet);
							} else {
								lrrValues = lrrs;
							}
							updateGUI();
						}
					}
				}
			}
		};
		for (int i = 0; i<Transforms.TRANFORMATIONS.length; i++) {
			transformationRadioButtons[i] = new JRadioButton(Transforms.TRANFORMATIONS[i], false);
			transformationRadioButtons[i].setFont(new Font("Arial", 0, 14));
			typeRadio.add(transformationRadioButtons[i]);
			transformationRadioButtons[i].addItemListener(typeListener);
//			transformationRadioButtons[i].setBackground(BACKGROUND_COLOR);
			transformationPanel.add(transformationRadioButtons[i]);
		}
		transformationRadioButtons[0].setSelected(true);
		descrPanel.add(transformationPanel);
		
		JPanel scopePanel = new JPanel();
		label = new JLabel("Transform by: ");
		label.setFont(new Font("Arial", 0, 14));
		scopePanel.add(label);
		ButtonGroup scopeRadio = new ButtonGroup();
		JRadioButton[] scopeRadioButtons = new JRadioButton[2];
		ItemListener scopeListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButton jrb = (JRadioButton)ie.getItem();
				if (jrb.isSelected()) {
					transformSeparatelyByChromosome = jrb.getText().equals(Transforms.SCOPES[1]); 
					if (transformation_type > 0) {
						lrrValues = Transforms.transform(lrrs, transformation_type, transformSeparatelyByChromosome, markerSet);
					} else {
						lrrValues = lrrs;
					}
					updateGUI();
				}
			}
		};
		for (int i = 0; i<Transforms.SCOPES.length; i++) {
			scopeRadioButtons[i] = new JRadioButton(Transforms.SCOPES[i], false);
			scopeRadioButtons[i].setFont(new Font("Arial", 0, 14));
			scopeRadio.add(scopeRadioButtons[i]);
			scopeRadioButtons[i].addItemListener(scopeListener);
//			scopeRadioButtons[i].setBackground(BACKGROUND_COLOR);
			scopePanel.add(scopeRadioButtons[i]);
		}
		scopeRadioButtons[0].setSelected(true);
		descrPanel.add(scopePanel);

		getContentPane().add(descrPanel, BorderLayout.NORTH);

		InputMap inputMap = panel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		// inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP,
		// InputEvent.ALT_MASK), ALT_UP);
		// inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN,
		// InputEvent.ALT_MASK), ALT_DOWN);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS_REGION);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT_REGION);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, InputEvent.CTRL_MASK), PREVIOUS_CHR);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, InputEvent.CTRL_MASK), NEXT_CHR);

		ActionMap actionMap = panel.getActionMap();
		actionMap.put(FIRST_CHR, new AbstractAction() {
			public static final long serialVersionUID = 5L;

			public void actionPerformed(ActionEvent e) {
				parseLocation("chr1");
			}
		});
		actionMap.put(PREVIOUS_CHR, new AbstractAction() {
			public static final long serialVersionUID = 6L;

			public void actionPerformed(ActionEvent e) {
				parseLocation("chr"+Math.max(chr-1, 1));
			}
		});
		actionMap.put(NEXT_CHR, new AbstractAction() {
			public static final long serialVersionUID = 7L;

			public void actionPerformed(ActionEvent e) {
				parseLocation("chr"+Math.min(chr+1, 25));
			}
		});
		actionMap.put(LAST_CHR, new AbstractAction() {
			public static final long serialVersionUID = 8L;

			public void actionPerformed(ActionEvent e) {
				parseLocation("chr25");
			}
		});
		actionMap.put(FIRST_REGION, new AbstractAction() {
			public static final long serialVersionUID = 9L;

			public void actionPerformed(ActionEvent e) {
				regionIndex = 0;
				showRegion();
			}
		});
		actionMap.put(PREVIOUS_REGION, new AbstractAction() {
			public static final long serialVersionUID = 10L;

			public void actionPerformed(ActionEvent e) {
				regionIndex = Math.max(regionIndex-1, 0);
				showRegion();
			}
		});
		actionMap.put(NEXT_REGION, new AbstractAction() {
			public static final long serialVersionUID = 11L;

			public void actionPerformed(ActionEvent e) {
				regionIndex = Math.min(regionIndex+1, regions.length-1);
				showRegion();
			}
		});
		actionMap.put(LAST_REGION, new AbstractAction() {
			public static final long serialVersionUID = 12L;

			public void actionPerformed(ActionEvent e) {
				regionIndex = regions.length-1;
				showRegion();
			}
		});
		panel.setActionMap(actionMap);

		// doesn't seem to get captured properly...
		nextChr.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT_REGION);
		nextChr.setActionMap(actionMap);
		previousChr.setActionMap(actionMap);
		
		
		
	}
	
	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();

		if (command.equals(FIRST_CHR)) {
			parseLocation("chr1");
		} else if (command.equals(PREVIOUS_CHR)) {
			parseLocation("chr"+Math.max(chr-1, 1));
		} else if (command.equals(NEXT_CHR)) {
			parseLocation("chr"+Math.min(chr+1, 25));
		} else if (command.equals(LAST_CHR)) {
			parseLocation("chr25");
		} else if (command.equals(FIRST_REGION)) {
			regionIndex = 0;
			showRegion();
		} else if (command.equals(PREVIOUS_REGION)) {
			regionIndex = Math.max(regionIndex-1, 0);
			showRegion();
		} else if (command.equals(NEXT_REGION)) {
			System.out.println("next");
			if (regions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded; files include: "+Array.toStr(proj.getFilenames(Project.REGION_LIST_FILENAMES), ", "), "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			regionIndex = Math.min(regionIndex+1, regions.length-1);
			showRegion();
		} else if (command.equals(LAST_REGION)) {
			regionIndex = regions.length-1;
			showRegion();
		} else {
			System.err.println("Error - unknown command '"+command+"'");
		}
	}

	public void mouseClicked(MouseEvent e) {
		if (e.getButton()==MouseEvent.BUTTON1) {
			if (leftClick!=null && leftClick.isAlive()) {
				leftClick.cancel();
				leftClick = null;
				zoomProportionally(false, e.getPoint(), true);
			} else {
				new Thread(leftClick = new SingleClick(this, MouseEvent.BUTTON1, DOUBLE_CLICK_INTERVAL)).start();
			}

		} else if (e.getButton()==MouseEvent.BUTTON3) {
			if (rightClick!=null && rightClick.isAlive()) {
				rightClick.cancel();
				rightClick = null;
				zoom(1, 1);
			} else {
				new Thread(rightClick = new SingleClick(this, MouseEvent.BUTTON3, DOUBLE_CLICK_INTERVAL)).start();
			}
		}
	}

	public void singleLeftClick() {
	// System.out.println("Single left click");
	}

	public void singleRightClick() {
	// System.out.println("Single right click");
	}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void mousePressed(MouseEvent e) {
		startX = e.getPoint().x;
		inDrag = true;
	}

	public void mouseReleased(MouseEvent e) {
		inDrag = false;
	}

	public void mouseDragged(MouseEvent e) {
		int curX = e.getPoint().x;
		int distance = startX-curX;

		distance *= (stop-start)/(getWidth()-2*WIDTH_BUFFER);

		if (distance<0) {
			distance = Math.max(distance, 1-start);
		} else {
			distance = Math.min(distance, positions[chrBoundaries[chr][1]]-stop);
		}

		if ((start<=1&&distance<0)||(stop>=positions[chrBoundaries[chr][1]]&&distance>0)) {

		} else {
			start += distance;
			stop += distance;
		}

		if (inDrag) {
			updateGUI();
			startX = curX;
		}
	}

	public void mouseMoved(MouseEvent e) {}

	public void mouseWheelMoved(MouseWheelEvent e) {
		zoomProportionally(e.getWheelRotation()>0, e.getPoint(), false);
	}

	public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
		int width = lrrPanel.getWidth()-2*WIDTH_BUFFER;
		double x = p.getX()-WIDTH_BUFFER;
		double xHat = width-x;
		double left = x/width;
		double right = xHat/width;
		double multiplier = MOUSE_WHEEL_MULTIPLIER/(outNotIn?1:-2);

		if (!outNotIn&&center) {
			right = 0.25-right;
			left = 0.25-left;
			multiplier = MOUSE_WHEEL_MULTIPLIER;
		}

		zoom(left*multiplier, right*multiplier);

	}

	public void zoom(double leftProportion, double rightProportion) {
		int dist = stop-start;
		start = start-(int)(leftProportion*dist);
		stop = stop+(int)(rightProportion*dist);
		updateGUI();
	}

	public void createSampleList() {
		long time = new Date().getTime();
		String[] filesPresent = Files.list(proj.getDir(Project.SAMPLE_DIRECTORY), Sample.SAMPLE_DATA_FILE_EXTENSION, jar);
		FontMetrics fontMetrics = sampleList.getFontMetrics(sampleList.getFont());
		String refresh = "refresh list";
		int maxWidth = fontMetrics.stringWidth(refresh);

		System.out.println("Determined sample list in "+ext.getTimeElapsed(time));

		if (filesPresent == null || filesPresent.length == 0) {
			samplesPresent = new String[] {proj.get(Project.SAMPLE_DIRECTORY)+" directory is empty", refresh};
			maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(samplesPresent[0]));
		} else {
			samplesPresent = new String[filesPresent.length+1];
			for (int i = 0; i<filesPresent.length; i++) {
				samplesPresent[i] = filesPresent[i].substring(0, filesPresent[i].lastIndexOf("."));
				maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(samplesPresent[i]));
            }
			samplesPresent[filesPresent.length] = refresh;
		}
		
		sampleList.setModel(new DefaultComboBoxModel<String>(samplesPresent));
		sampleList.setPreferredSize(new Dimension(maxWidth+50, 30));

	}

	public void loadMarkers() {
		Hashtable<String,String> hash;
		String[] markerNames;
		byte[] chrs;
		long time;
		int chr;

		time = new Date().getTime();

		hash = proj.getFilteredHash();
		markerSet = proj.getMarkerSet();
		fingerprint = markerSet.getFingerprint();
		markerNames = markerSet.getMarkerNames();
		numMarkers = markerNames.length;
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		dropped = new boolean[numMarkers];
		chrBoundaries = new int[27][2];
		for (int i = 0; i<chrBoundaries.length; i++) {
//			chrBoundaries[i][0] = chrBoundaries[i][1] = -1;
			chrBoundaries[i][0] = chrBoundaries[i][1] = 0;
		}
		chr = 0;
		for (int i = 0; i<numMarkers; i++) {
			dropped[i] = hash.containsKey(markerNames[i]);
			if (chrs[i]>chr) {
				if (chr!=0) {
					chrBoundaries[chr][1] = i-1;
				}
				chr = chrs[i];
				chrBoundaries[chr][0] = i;
			}
		}
		chrBoundaries[chr][1] = numMarkers-1;
		chrBoundaries[0][0] = 0;
		chrBoundaries[0][1] = numMarkers-1;

		System.out.println("Read in data for "+numMarkers+" markers in "+ext.getTimeElapsed(time));
	}

	public void loadValues() {
		Sample samp;

		samp = proj.getPartialSampleFromRandomAccessFile(sample);
		if (samp == null) {
			System.err.println("Error - sample '"+sample+"' not found in "+proj.getDir(Project.SAMPLE_DIRECTORY));
		} else if ( samp.getFingerprint()!=fingerprint) {
			System.err.println("Error - Sample "+proj.getDir(Project.SAMPLE_DIRECTORY)+sample+Sample.SAMPLE_DATA_FILE_EXTENSION+" has a different fingerprint ("+samp.getFingerprint()+") than the MarkerSet ("+fingerprint+")");
		} else {
			lrrs = samp.getLRRs();
			
			if (transformation_type > 0) {
				lrrValues = Transforms.transform(lrrs, transformation_type, transformSeparatelyByChromosome, markerSet);
			} else {
				lrrValues = lrrs;
			}
			
			bafs = samp.getBAFs();
			genotypes = samp.getAB_Genotypes();
		}
		// lrrMin = Math.floor(lrrMin);
		// lrrMax = Math.ceil(lrrMax);
	}

	public void parseLocation(String location) {
		byte oldChr;
		int[] loc;

		oldChr = chr;

		if (location == null) {
			System.err.println("Error - null location");
			return;
		}
		
		if (!location.startsWith("chr")) {
			if (track == null) {
				JOptionPane.showMessageDialog(this, "Cannot parse '"+location+"' since the gene track has either not been installed or did not load properly.", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			} else {
				loc = track.lookupPosition(location);
				if (loc[0] == -1) {
					JOptionPane.showMessageDialog(this, "'"+location+"' is not a valid gene name and is not a valid UCSC location.", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
			}
		} else {
			loc = Positions.parseUCSClocation(location);
		}
		
		if (loc == null) {
			JOptionPane.showMessageDialog(this, "'"+location+"' is not a valid UCSC location.", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}

		chr = (byte)loc[0];
		if (chr==-1) {
			chr = oldChr;
			return;
		}
		if (chr!=oldChr) {
			start = stop = -1;
			procCNVs(chr);
		}
		start = loc[1];
		stop = loc[2];
		if (start==-1||start<0) {
			start = 1;
		}
		if (stop==-1||stop>positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
		}
		
		updateGUI();
	}

	public void updateSample(String newSample) {
		boolean found;
		long time;

		found = sampleList.getSelectedItem().equals(newSample);
		if (found) {
			time = new Date().getTime();
//			System.out.println("Loading CNVs...");
//			System.out.println("took "+ext.getTimeElapsed(time));
			sample = newSample;
			time = new Date().getTime();
			System.out.print("Found "+sample+"...");
			indiPheno = sampleData.getIndiPheno(sample.toLowerCase());
			if (indiPheno == null) {
				if (!sample.equals(proj.get(Project.SAMPLE_DIRECTORY)+" directory is empty")) {
					JOptionPane.showMessageDialog(this, "Sample '"+sample+"' was not present in the SampleData file", "Error", JOptionPane.ERROR_MESSAGE);
				}
				return;
			}
			loadValues();
			procCNVs(chr);
			if (regionsList.length == 0 || !Files.exists(regionsList[regionsListIndex], jar)) {
				loadCNVsAsRegions();
			}
			updateGUI();
			System.out.println("updated in "+ext.getTimeElapsed(time));
		} else {
			for (int i = 0; i<samplesPresent.length && !found; i++) {
				if (samplesPresent[i].equals(newSample)) {
					sampleList.setSelectedIndex(i);
					found = true;
				}
	        }
			if (!found) {
				if (Files.exists(proj.getDir(Project.SAMPLE_DIRECTORY)+newSample+Sample.SAMPLE_DATA_FILE_EXTENSION, jar)) {
					createSampleList();
					updateSample(newSample);
				} else {
					JOptionPane.showMessageDialog(this, "Data was not found for the next sample in the list ("+newSample+").", "Error", JOptionPane.ERROR_MESSAGE);
				}
			}
		}
		
	}
	
	public void updateGUI() {
		if (start<=1) {
			startMarker = chrBoundaries[chr][0];
			start = 1;
		} else {
			startMarker = Array.binarySearch(positions, start, chrBoundaries[chr][0], chrBoundaries[chr][1], false);
		}

		if (stop>=positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
			stopMarker = chrBoundaries[chr][1];
		} else {
			stopMarker = Array.binarySearch(positions, stop, chrBoundaries[chr][0], chrBoundaries[chr][1], false);
		}

		if (startMarker==-1) {
			System.err.println("Error - failed to find startMarker");
			startMarker = chrBoundaries[chr][0];
		}

		if (stopMarker==-1) {
			System.err.println("Error - failed to find stopMarker");
			stopMarker = chrBoundaries[chr][1];
		}
		
		displayIndex();
		repaint();
	}

	public void loadRegions() {
		BufferedReader reader;
        Vector<String[]> v;
        String[] line;
        int ignoredLines;
		
		try {
			reader = Files.getReader(regionsList[regionsListIndex], jar, false, false);
			System.out.print("Loading regions from "+regionsList[regionsListIndex]+"...");
	        v = new Vector<String[]>();
	        ignoredLines = 0;
            while (reader.ready()) {
            	line = reader.readLine().trim().split("\t");
            	if (line.length > 1 && line[1].startsWith("chr")) {
            		v.add(line);
            	} else {
            		ignoredLines++;
            	}
            }
            System.out.println(" loaded "+v.size()+" regions");
            regions = Matrix.toStringArrays(v);
            if (ignoredLines > 1) {
            	JOptionPane.showMessageDialog(null, "Error - there were "+ignoredLines+" regions in '"+regionsList[regionsListIndex]+"' that were ignored due to improper formatting", "Error", JOptionPane.ERROR_MESSAGE);
            }
            reader.close();
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \""+regionsList[regionsListIndex]+"\" not found in data directory");
            System.exit(1);
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+regionsList[regionsListIndex]+"\"");
            System.exit(2);
        }
	}
	
	public void loadCNVsAsRegions() {
		int buffer;
		Segment[][] segs;
		int count;
		
		count = 0;
		segs = new Segment[26][];
		for (int i = 0; i<segs.length; i++) {
			procCNVs((byte)i);
			segs[i] = findUniqueRegions(cnvs);
			count += segs[i].length;
        }

		regions = new String[count][2];
		count = 0;
		for (int i = 0; i<segs.length; i++) {
			for (int j = 0; j<segs[i].length; j++) {
				regions[count][0] = sample;
				buffer = Math.max(segs[i][j].getSize() / 2, MIN_BUFFER);
				regions[count][1] = Positions.getUCSCformat(new int[] {segs[i][j].getChr(), segs[i][j].getStart()-buffer, segs[i][j].getStop()+buffer});
				count++;
            }
		}
	}

	public static Segment[] findUniqueRegions(CNVariant[][] cnvs) { // haven't actually coded the collapsing of the segments yet
		Segment[] segs;
		int count;
		
		count = 0;
		for (int i = 0; i<cnvs.length; i++) {
			count += cnvs[i].length;
        }
		
		segs = new Segment[count];
		count = 0;
		for (int i = 0; i<cnvs.length; i++) {
			for (int j = 0; j<cnvs[i].length; j++) {
				segs[count] = cnvs[i][j];
				count++;
            }
		}

		return segs;
	}

	public void showRegion() {
		System.out.println("regionIndex="+regionIndex+"\t"+"regions.length="+regions.length);
		parseLocation(regions[regionIndex][1]);
		if (regions[regionIndex].length > 2) {
			commentLabel.setText("region #"+(regionIndex+1)+":  "+ regions[regionIndex][2]);
		}
		
		if (!regions[regionIndex][0].equals(sample)) {
			updateSample(regions[regionIndex][0]);
		}
	}

//	public static Vector<String[]> loadFileToVec(String filename, boolean demo, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent) {
//		BufferedReader reader = null;
//		Vector<String[]> v = new Vector<String[]>();
//		String[] trav;
//		String[] line;
//
//		try {
//			if (demo) {
//				reader = new BufferedReader(new InputStreamReader(ClassLoader.getSystemResourceAsStream(filename)));
//			} else {
//				reader = new BufferedReader(new FileReader(filename));
//			}
//			if (ignoreFirstLine) {
//				reader.readLine();
//			}
//			while (reader.ready()) {
//				trav = line = reader.readLine().trim().split("\t", -1);
//				if (cols!=null) {
//					trav = new String[cols.length];
//					for (int i = 0; i<cols.length; i++) {
//						trav[i] = line[cols[i]];
//					}
//				}
//				if (onlyIfAbsent) {
//					if (!v.contains(trav)) {
//						v.add(trav);
//					}
//				} else {
//					v.add(trav);
//				}
//			}
//			reader.close();
//		} catch (FileNotFoundException fnfe) {
//			System.err.println("Error: file \""+filename+"\" not found in current directory");
//			System.exit(1);
//		} catch (Exception e) {
//			System.err.println("Error reading file \""+filename+"\"");
//			e.printStackTrace();
//			System.exit(2);
//		}
//
//		return v;
//	}

	public void displayIndex() {
		if (start<0) {
			start = 1;
		}
		if (stop>positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
		}

		navigationField.setText(chr==0?"all":"chr"+chr+":"+ext.addCommas(start)+"-"+ext.addCommas(stop));
	}

//	public static CNVariant[] getCNVlist(Hashtable<String,Hashtable<String,CNVariant[]>> cnvHashes, String[] filenames, String sample) {
//		Hashtable<String,CNVariant[]> hash;
//		Vector<CNVariant> v;
//		CNVariant[] set;
//		
//		v = new Vector<CNVariant>();
//		
//		for (int i = 0; i<filenames.length; i++) {
//			hash = cnvHashes.get(filenames[i]);
//			if (hash != null) {
//				set = hash.get(sample);
//				if (set != null) {
//					for (int j = 0; j<set.length; j++) {
//						set[j].setSource(i);
//						v.add(set[j]);
//                    }
//				}
//			}
//        }
//		
//		return CNVariant.toArray(v);
//	}
	
	public void procCNVs(byte chr) {
		cnvs = new CNVariant[cnvLabels.length][];
		for (int i = 0; i<cnvLabels.length; i++) {
			cnvs[i] = indiPheno.getCNVs(i, chr);
			if (cnvs[i] == null) {
				cnvs[i] = new CNVariant[0];
			}
			System.out.println("Proccessed "+cnvs[i].length+" "+cnvLabels[i]+" CNVs for chr"+chr);
        }
	}
	
//	public static CNVariant[] loadCNVfiles(Project proj, String[] filenames) {
//		BufferedReader reader;
//		Vector<CNVariant> v = null;
//		String[] line;
//		String id = null;
//		boolean jar;
//
//		jar = proj.getJarStatus();
//		
//		v = new Vector<CNVariant>();
//		for (int i = 0; i<filenames.length; i++) {
//			try {
//				reader = Files.getReader(filenames[i], jar, true, false);
//				if (reader!=null) {
//					reader.mark(1000);
//					line = reader.readLine().trim().split("[\\s]+");
//					if (!line[2].toLowerCase().equals("chr")&&ext.chromosomeNumber(line[2])!=-1) {
//						reader.reset();
//					}
//					while (reader.ready()) {
//						line = reader.readLine().trim().split("[\\s]+");
//						if (line[1].equals(id)) {
//							v.add(new CNVariant(line, i));
//						}
//					}
//					reader.close();
//				}
//			} catch (IOException ioe) {
//				System.err.println("Error reading file \""+filenames[i]+"\"");
//				ioe.printStackTrace();
//			}
//
//		}
//
//		return CNVariant.sortCNVs(CNVariant.toArray(v));
//	}

	
	public static void main(String[] args) {
		Project proj;
		boolean jar;
		
		jar = !(args.length>0&&args[0].equals("-notJar"));
		proj = new Project(Project.DEFAULT_PROJECT, jar);
		new Trailer(proj, DEFAULT_SAMPLE, proj.getFilenames(Project.CNV_FILENAMES), DEFAULT_LOCATION);
	}
}
