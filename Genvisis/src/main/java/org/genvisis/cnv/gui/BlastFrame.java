package org.genvisis.cnv.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerSeqAnnotation;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.util.CNVHelper;
import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.StrandOps;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;
import net.miginfocom.swing.MigLayout;

public class BlastFrame extends JFrame implements WindowFocusListener {

	public static class BlastUtils {

		public static int countAlignment(BlastAnnotation annotation) {
			int align = 0;
			for (CigarElement ce : annotation.getCigar().getCigarElements()) {
				switch (ce.getOperator()) {
					case D:
						break;
					case EQ:
						align += ce.getLength();
						break;
					case H:
						break;
					case I:
						break;
					case M:
						break;
					case N:
						break;
					case P:
						break;
					case S:
						break;
					case X:
						break;
					default:
						break;
				}

			}
			return align;
		}

		public static ArrayList<BlastAnnotation> filterAnnotations(Project proj,
																															 List<BlastAnnotation> annotations,
																															 int alignFilter) {
			ArrayList<BlastAnnotation> filteredList = new ArrayList<BlastAnnotation>();
			for (BlastAnnotation annot : annotations) {
				if (countAlignment(annot) > alignFilter) {
					filteredList.add(annot);
				}
			}
			return filteredList;
		}

		public static Segment getSegmentForAnnotation(MarkerSeqAnnotation ref, BlastAnnotation annot) {
			Segment origSeg = annot.getRefLoc();
			Cigar cig = annot.getCigar();
			CigarElement cigStart = cig.getCigarElement(0);
			CigarElement cigStop = cig.getCigarElement(cig.numCigarElements() - 1);

			int start = origSeg.getStart();
			if (cigStart.getOperator() == CigarOperator.X) {
				start -= cigStart.getLength();
			}
			int stop = origSeg.getStop();
			if (cigStop.getOperator() == CigarOperator.X) {
				stop += cigStop.getLength();
			}

			return new Segment(origSeg.getChr(), start, stop);
		}

		// public static String getAB(Project proj, int chr, int loc) {
		// String abCode = null;
		// String[] samples = proj.getSamples();
		// MarkerSet markerSet = proj.getMarkerSet(); // TODO should be cached
		// // TODO alter 'loc' to reflect position of next allele after PROBE_SEQ
		// int[] mkrPoss = markerSet.getPositionsByChr()[chr]; // TODO should be cached
		// int pos = Array.binarySearch(mkrPoss, loc, true); // TODO exact match?? should be okay with
		// next marker, or must be exact next base?
		// if (pos != -1) {
		// int mkrInd = markerSet.getIndicesByChr()[chr][pos];
		// HashSet<String> alleleSet1 = new HashSet<String>();
		// HashSet<String> alleleSet2 = new HashSet<String>();
		// for (String sample : samples) {
		// Sample samp = proj.getPartialSampleFromRandomAccessFile(sample, false, false, false, false,
		// true); // TODO should be cached?
		// byte[] geno = samp.getForwardGenotypes();
		// if (geno != null) {
		// String alPair = Sample.ALLELE_PAIRS[geno[mkrInd]];
		// alleleSet1.add("" + alPair.charAt(0));
		// alleleSet2.add("" + alPair.charAt(1));
		// }
		// }
		// boolean a1_A = alleleSet1.contains("A");
		// boolean a1_C = alleleSet1.contains("C");
		// boolean a1_T = alleleSet1.contains("T");
		// boolean a1_G = alleleSet1.contains("G");
		// boolean a2_A = alleleSet1.contains("A");
		// boolean a2_C = alleleSet1.contains("C");
		// boolean a2_T = alleleSet1.contains("T");
		// boolean a2_G = alleleSet1.contains("G");
		//
		// boolean a1_AT = a1_A || a1_T;
		// boolean a1_CG = a1_C || a1_G;
		// boolean a2_AT = a2_A || a2_T;
		// boolean a2_CG = a2_C || a2_G;
		// if (a1_AT && !a1_CG && !a2_AT && a2_CG) {
		// abCode = "AB";
		// } else if (!a1_AT && a1_CG && a2_AT && !a2_CG) {
		// abCode = "BA";
		// }
		// // TODO verify AB encoding
		// } else {
		// // TODO couldn't find position of next base - should log?
		// abCode = "";
		// }
		// return abCode;
		// }

	}

	private static final long serialVersionUID = 1L;

	private final Project proj;
	private MarkerSeqAnnotation referenceAnnotation;
	private ArrayList<BlastAnnotation> annotations;
	private ReferenceGenome referenceGenome;
	private final JLabel refLabel;
	private JLabel strandLbl;
	private JLabel probeLengthLbl;
	private JLabel locationLbl;
	private JLabel probeDescLbl;
	private final ReferenceLabel probeLbl;

	private JLabel probeShiftLbl;
	private JLabel nextBaseLbl;
	private final JPanel contentPane;
	private final JCheckBox chckbxExpandBlastResults;
	private final JButton btnClose;
	private final JSpinner spnFontSize;
	private final JScrollPane scrollPane;
	private final JPanel blastPanel;
	private int rowCnt = 0;
	private final JCheckBox chckbxSortByLocation;
	private final JSeparator separator_1;
	private final JCheckBox chckbxPinToFront;
	HashMap<BlastLabel, JLabel> alignCntLblMap = new HashMap<BlastLabel, JLabel>();
	private ArrayList<BlastLabel> bLabels;
	private ArrayList<BlastLabel> sorted;
	private ArrayList<JLabel> otherLabels;
	private final JSpinner spinnerAlignmentLength;
	private final JLabel lblAlignmentLength;
	private final JSeparator separator_2;
	private MarkerBlastAnnotation blastResult;
	public int currentAlignFilter;
	private final JSplitPane splitPane;
	private JPanel lowerPanel;

	private JLabel abLbl;
	String BLAST_COL_DEF = "[250px][20px][20px][grow][20px:20px:20px][120px:120px:120px]";

	/**
	 * Create the frame.
	 */
	public BlastFrame(Project proj, ChangeListener alignLengthListener) {
		super("Genvisis - BlastViewer - " + proj.PROJECT_NAME.getValue());
		this.proj = proj;
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		setMinimumSize(new Dimension(100, 100));
		UITools.setSize(this, new Dimension(1150, 800));
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);

		JPanel panel = new JPanel();
		panel.setBorder(null);
		contentPane.add(panel, BorderLayout.SOUTH);
		panel.setLayout(new MigLayout("", "[133px][][][][20px][20px][][][grow][][]", "[23px]"));

		chckbxExpandBlastResults = new JCheckBox(new AbstractAction() {
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				BlastLabel.expanded = !BlastLabel.expanded;
				BlastFrame.this.repaint();
			}
		});
		chckbxExpandBlastResults.setText("Expand BLAST Results");
		chckbxExpandBlastResults.setSelected(BlastLabel.expanded);
		chckbxExpandBlastResults.setVerticalAlignment(SwingConstants.BOTTOM);
		panel.add(chckbxExpandBlastResults, "cell 0 0,alignx left,aligny top");

		separator_1 = new JSeparator();
		separator_1.setOrientation(SwingConstants.VERTICAL);
		panel.add(separator_1, "cell 1 0,growy");

		chckbxSortByLocation = new JCheckBox(new AbstractAction() {
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				BlastFrame.this.clearLabels();
				if (((javax.swing.JCheckBox) e.getSource()).isSelected()) {
					for (BlastLabel bl : sorted) {
						BlastFrame.this.addBlastLabel(bl);
					}
				} else {
					for (BlastLabel bl : bLabels) {
						BlastFrame.this.addBlastLabel(bl);
					}
				}
				BlastFrame.this.validate();
				BlastFrame.this.repaint();
			}
		});
		chckbxSortByLocation.setText("Sort By Location");
		chckbxSortByLocation.setVerticalAlignment(SwingConstants.BOTTOM);
		panel.add(chckbxSortByLocation, "cell 2 0");

		JSeparator separator = new JSeparator();
		separator.setOrientation(SwingConstants.VERTICAL);
		panel.add(separator, "cell 4 0,alignx left,growy");

		spnFontSize = new JSpinner(new SpinnerNumberModel(15, 5, 25, 1));
		spnFontSize.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				BlastLabel.setFontSize((Integer) spnFontSize.getModel().getValue());
				refLabel.setFont(BlastLabel.LBL_FONT);
				probeLbl.setFont(BlastLabel.LBL_FONT);
				locationLbl.setFont(BlastLabel.LBL_FONT);
				probeDescLbl.setFont(BlastLabel.LBL_FONT);
				strandLbl.setFont(BlastLabel.LBL_FONT);
				probeLengthLbl.setFont(BlastLabel.LBL_FONT);
				nextBaseLbl.setFont(BlastLabel.LBL_FONT);
				abLbl.setFont(BlastLabel.LBL_FONT);
				probeShiftLbl.setFont(BlastLabel.LBL_FONT);
				for (BlastLabel lbl : bLabels) {
					lbl.setFont(BlastLabel.LBL_FONT);
				}
				for (JLabel lbl : otherLabels) {
					lbl.setFont(BlastLabel.LBL_FONT);
				}
				BlastFrame.this.repaint();
			}
		});
		panel.add(spnFontSize, "flowx,cell 5 0,alignx left,aligny center");

		JLabel lblFontSize = new JLabel("Font Size");
		panel.add(lblFontSize, "cell 5 0");

		separator_2 = new JSeparator();
		separator_2.setOrientation(SwingConstants.VERTICAL);
		panel.add(separator_2, "cell 6 0,growy");

		int pLen = proj.ARRAY_TYPE.getValue().getProbeLength();
		double pRat = proj.BLAST_PROPORTION_MATCH_FILTER.getValue();
		int pRatLen = (int) (pLen * pRat);
		spinnerAlignmentLength = new JSpinner(new SpinnerNumberModel(pRatLen, 1, pLen, 1));
		spinnerAlignmentLength.addChangeListener(alignLengthListener);
		panel.add(spinnerAlignmentLength, "flowx,cell 7 0");

		chckbxPinToFront = new JCheckBox("Pin to Front");
		panel.add(chckbxPinToFront, "cell 9 0");

		btnClose = new JButton();
		btnClose.setAction(new AbstractAction() {
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent e) {
				BlastFrame.this.setVisible(false);
				// BlastFrame.this.dispose(); // removed for state
			}
		});
		btnClose.setText("Close");
		panel.add(btnClose, "cell 10 0");

		lblAlignmentLength = new JLabel("Alignment Length");
		panel.add(lblAlignmentLength, "cell 7 0");

		splitPane = new JSplitPane();
		splitPane.setOrientation(JSplitPane.VERTICAL_SPLIT);
		contentPane.add(splitPane, BorderLayout.CENTER);

		blastPanel = new JPanel();
		blastPanel.setLayout(new MigLayout("", BLAST_COL_DEF, ""));
		scrollPane = new JScrollPane(blastPanel);
		scrollPane.getVerticalScrollBar().setUnitIncrement(16);
		splitPane.setLeftComponent(scrollPane);
		refLabel = new ReferenceLabel();
		probeLbl = new ReferenceLabel();
		JPanel hdrPanel = getHeaderPanel();
		scrollPane.setColumnHeaderView(hdrPanel);

		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);

		lowerPanel = new JPanel();
		splitPane.setRightComponent(lowerPanel);

		addWindowFocusListener(this);
		pack();
	}

	public void addBlastLabel(BlastLabel lbl) {
		JLabel locLbl = new JLabel();
		locLbl.setText(lbl.fullSegment.getUCSClocation());
		Font lblFont = BlastLabel.LBL_FONT;// Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
		locLbl.setFont(lblFont);
		otherLabels.add(locLbl);
		// locLbl.setToolTipText(lbl.myAnnotation.getCigar().toString());
		blastPanel.add(locLbl, "cell 0 " + rowCnt);
		JLabel strandLbl = new JLabel();
		strandLbl.setText(CNVHelper.decode(lbl.getStrand()));
		strandLbl.setFont(lblFont);
		otherLabels.add(strandLbl);
		alignCntLblMap.put(lbl, strandLbl);
		blastPanel.add(strandLbl, "cell 1 " + rowCnt);
		JLabel probeLengthLbl = new JLabel();
		probeLengthLbl.setText(lbl.getAlignment() + "");
		probeLengthLbl.setFont(lblFont);
		otherLabels.add(probeLengthLbl);
		blastPanel.add(probeLengthLbl, "cell 2 " + rowCnt);
		blastPanel.add(lbl, "grow, cell 3 " + rowCnt);

		JLabel nextBaseLbl = new JLabel();
		nextBaseLbl.setFont(lblFont);
		Segment seg = lbl.fullSegment;
		int len = proj.ARRAY_TYPE.getValue().getProbeLength();
		int start = referenceAnnotation.goLeft() ? seg.getStart() - len : seg.getStart();
		int stop = referenceAnnotation.goLeft() ? seg.getStart() : seg.getStart() + len;
		// locationLbl.setText(seg.getChromosomeUCSC() + ":" + start + "-" + stop);
		String[] gen = referenceGenome.getSequenceFor(new Segment(seg.getChr(), start, stop));
		nextBaseLbl.setText(lbl.positiveStrand ? gen[gen.length - 1] : gen[0]);
		otherLabels.add(nextBaseLbl);
		blastPanel.add(nextBaseLbl, "alignx left, cell 4 " + rowCnt);

		JLabel abLbl = new JLabel();
		String ab = "";
		if (referenceAnnotation != null) {
			switch (lbl.getAB()) {
				case A:
					ab = StrandOps.flipsIfNeeded(referenceAnnotation.getA().getBaseString(),
																			 referenceAnnotation.getStrand(), false);
					break;
				case B:
					ab = StrandOps.flipsIfNeeded(referenceAnnotation.getB().getBaseString(),
																			 referenceAnnotation.getStrand(), false);
					break;
				case BOTH:
					// ab = StrandOps.flipsIfNeeded(referenceAnnotation.getA().getBaseString(),
					// referenceAnnotation.getStrand(), false) + "/" +
					// StrandOps.flipsIfNeeded(referenceAnnotation.getB().getBaseString(),
					// referenceAnnotation.getStrand(), false);
					break;
			}
		}
		abLbl.setText(ab);
		abLbl.setFont(lblFont);
		abLbl.setHorizontalAlignment(SwingConstants.CENTER);
		abLbl.setHorizontalTextPosition(SwingConstants.CENTER);
		otherLabels.add(abLbl);
		blastPanel.add(abLbl, "alignx center, cell 5 " + rowCnt);
		rowCnt++;
	}

	public void clearLabels() {
		blastPanel.removeAll();
		otherLabels.clear();
	}

	public void setAnnotations(MarkerBlastAnnotation blastResult, ReferenceGenome refGen) {
		referenceAnnotation = blastResult.getMarkerSeqAnnotation();
		if (blastResult.hasPerfectMatch(proj.getLog())) {
			referenceAnnotation.setStrand(blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH,
																																	proj.getLog())
																							 .get(0).getStrand());
			referenceAnnotation.setSegment(blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH,
																																	 proj.getLog())
																								.get(0).getRefLoc());
		}
		this.blastResult = blastResult;
		referenceGenome = refGen;
		double filter = proj.BLAST_PROPORTION_MATCH_FILTER.getValue();
		int probe = proj.ARRAY_TYPE.getValue().getProbeLength();
		int alignFilter = (int) (filter * probe);
		updateAnnotations(alignFilter);
		updateLabels();
	}

	public void setSecondPanel(JPanel img) {
		lowerPanel = img;
		lowerPanel.setBorder(BorderFactory.createLineBorder(Color.GRAY, 1));
		splitPane.setRightComponent(lowerPanel);
		splitPane.setDividerLocation(0.5);
		splitPane.setOneTouchExpandable(true);
		validate();
		repaint();
	}

	public void updateAnnotations(int alignFilter) {
		currentAlignFilter = alignFilter;
		annotations = BlastFrame.BlastUtils.filterAnnotations(proj,
																													blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS,
																																												proj.getLog()),
																													alignFilter);
	}

	public void updateLabels() {
		BlastLabel.spaces = new TreeSet<Integer>();
		BlastLabel.spaceSets = new TreeMap<Integer, Integer>();
		bLabels = new ArrayList<BlastLabel>();
		sorted = new ArrayList<BlastLabel>();
		otherLabels = new ArrayList<JLabel>();
		for (BlastAnnotation annot : annotations) {
			BlastLabel lbl = new BlastLabel(proj, referenceAnnotation, annot, referenceGenome);
			// if (lbl.strandFlipped) continue;
			bLabels.add(lbl);
			sorted.add(lbl);
		}
		Collections.sort(sorted, new Comparator<BlastLabel>() {
			@Override
			public int compare(BlastLabel arg0, BlastLabel arg1) {
				int res1 = (new Byte(arg0.fullSegment.getChr())).compareTo(new Byte(arg1.fullSegment.getChr()));
				if (res1 != 0) {
					return res1;
				}
				res1 = (new Integer(arg0.fullSegment.getStart())).compareTo(new Integer(arg1.fullSegment.getStart()));
				if (res1 != 0) {
					return res1;
				}
				res1 = (new Integer(arg0.fullSegment.getStop())).compareTo(new Integer(arg1.fullSegment.getStop()));
				return res1;
			}
		});
		List<BlastAnnotation> onT = blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT,
																															proj.getLog());
		probeLbl.setAnnotation(onT.size() > 0 ? onT.get(0) : null);

		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				BlastFrame.this.clearLabels();
				for (BlastLabel lbl : bLabels) {
					BlastFrame.this.addBlastLabel(lbl);
				}
				boolean posStrand = referenceAnnotation.getStrand() == Strand.POSITIVE;
				String refSeq = referenceAnnotation.getSequence();
				refLabel.setText(!posStrand ? new StringBuilder(refSeq).reverse().toString() : refSeq);
				strandLbl.setText(CNVHelper.decode(referenceAnnotation.getStrand()));
				int len = proj.ARRAY_TYPE.getValue().getProbeLength();
				probeLengthLbl.setText(len + "");
				List<BlastAnnotation> onT = blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT,
																																	proj.getLog());
				String lbl = "";
				// if (onT.size() > 0) {
				// if (onT.size() > 1) {
				// System.out.println("Warning - multiple On-Target Non-Perfect Matches Found. Displaying
				// information for first result.");
				// }
				// BlastAnnotation onTAnnot = onT.get(0);
				// int preCount = onTAnnot.getCigar().getCigarElement(0).getOperator() == CigarOperator.X ?
				// onTAnnot.getCigar().getCigarElement(0).getLength() : 0;
				// int postCount =
				// onTAnnot.getCigar().getCigarElement(onTAnnot.getCigar().numCigarElements() -
				// 1).getOperator() == CigarOperator.X ?
				// onTAnnot.getCigar().getCigarElement(onTAnnot.getCigar().numCigarElements() -
				// 1).getLength() : 0;
				// if (preCount > 0 && postCount > 0) {
				// lbl = preCount + "-> | <-" + postCount;
				// } else if (preCount > 0) {
				// lbl = preCount + "->";
				// } else if (postCount > 0) {
				// lbl = "<-" + postCount;
				// }
				// }
				probeShiftLbl.setText(lbl);
				{
					boolean hasMatch = blastResult.hasPerfectMatch(proj.getLog());
					onT.size();
					// if (hasMatch) {
					probeDescLbl.setText("<reference sequence:>");
					Segment seg = referenceAnnotation.getSeg();
					int start = hasMatch ? seg.getStart()
															 : referenceAnnotation.goLeft() ? seg.getStart() - len
																															: seg.getStart();
					int stop = hasMatch ? seg.getStop() : referenceAnnotation.goLeft() ? seg.getStart()
																																						 : seg.getStart() + len;
					if (posStrand) {
						stop = stop + 1;
					} else {
						start = start - 1;
					}
					locationLbl.setText(seg.getChromosomeUCSC() + ":" + start + "-" + stop);

					String[] gen = referenceGenome.getSequenceFor(new Segment(seg.getChr(), start, stop));
					String[] act = ArrayUtils.subArray(gen, posStrand ? 0 : 1,
																						 posStrand ? gen.length - 1 : gen.length);

					nextBaseLbl.setText(posStrand ? gen[gen.length - 1] : gen[0]);
					String seq = ArrayUtils.toStr(act, "");
					probeLbl.setText(posStrand ? seq : BlastLabel.flipBases(seq));
					// } else if (hasOnT) {
					// probeDescLbl.setText("<on-target mismatch:>");
					// Segment seg = BlastUtils.getSegmentForAnnotation(referenceAnnotation, onT.get(0));
					// int start = referenceAnnotation.goLeft() ? seg.getStart() - len : seg.getStart();
					// int stop = referenceAnnotation.goLeft() ? seg.getStart() : seg.getStart() + len;
					// locationLbl.setText(seg.getChromosomeUCSC() + ":" + start + "-" + stop);
					// String[] gen = referenceGenome.getSequenceFor(new Segment(seg.getChr(), start, stop));
					// String[] act = Array.subArray(gen, posStrand ? 0 : 1, posStrand ? gen.length - 1 :
					// gen.length);
					// nextBaseLbl.setText(posStrand ? gen[gen.length - 1] : gen[0]);
					// String seq = Array.toStr(act, "");
					// if (!referenceAnnotation.goLeft()) {
					// seq = new StringBuilder(seq).reverse().toString();
					// }
					// probeLbl.setText(posStrand ? seq : BlastLabel.flipBases(seq));
					// }
				}
				String a = StrandOps.flipsIfNeeded(referenceAnnotation.getA().getBaseString(),
																					 referenceAnnotation.getStrand(), false);
				String b = StrandOps.flipsIfNeeded(referenceAnnotation.getB().getBaseString(),
																					 referenceAnnotation.getStrand(), false);
				String abLblStr = "(A) " + a + " | (B) " + b;
				System.out.println(abLblStr);
				abLbl.setText(abLblStr);
				BlastFrame.this.validate();
				BlastFrame.this.repaint();
			}
		});
	}

	@Override
	public void windowGainedFocus(WindowEvent e) {}

	@Override
	public void windowLostFocus(WindowEvent e) {
		if (e.getNewState() != WindowEvent.WINDOW_CLOSED) {
			if (chckbxPinToFront != null && chckbxPinToFront.isSelected()) {
				setAlwaysOnTop(true);
			} else {
				setAlwaysOnTop(false);
			}
		}

	}

	private JPanel getHeaderPanel() {
		JPanel hdrPanel = new JPanel(new MigLayout("", BLAST_COL_DEF, "[][]"));
		hdrPanel.setBorder(null);
		locationLbl = new JLabel();
		Font lblFont = BlastLabel.LBL_FONT;// Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
		locationLbl.setFont(lblFont);
		hdrPanel.add(locationLbl, "cell 0 0");
		probeDescLbl = new JLabel("<reference sequence:>");
		probeDescLbl.setFont(lblFont);
		hdrPanel.add(probeDescLbl, "cell 0 1");
		strandLbl = new JLabel(referenceAnnotation == null ? ""
																											 : CNVHelper.decode(referenceAnnotation.getStrand()));
		strandLbl.setFont(lblFont);
		hdrPanel.add(strandLbl, "cell 1 0");
		probeLengthLbl = new JLabel(proj.ARRAY_TYPE.getValue().getProbeLength() + "");
		probeLengthLbl.setFont(lblFont);
		hdrPanel.add(probeLengthLbl, "cell 2 0");
		hdrPanel.add(refLabel, "grow, cell 3 0");
		nextBaseLbl = new JLabel();
		nextBaseLbl.setFont(lblFont);
		hdrPanel.add(nextBaseLbl, "alignx left, cell 4 0");
		String abLblStr = referenceAnnotation == null ? ""
																									: "(A) "
																										+ StrandOps.flipsIfNeeded(referenceAnnotation.getA()
																																																 .getBaseString(),
																																							referenceAnnotation.getStrand(),
																																							false)
																										+ "- (B) "
																										+ StrandOps.flipsIfNeeded(referenceAnnotation.getB()
																																																 .getBaseString(),
																																							referenceAnnotation.getStrand(),
																																							false);
		abLbl = new JLabel(abLblStr);
		abLbl.setFont(lblFont);
		hdrPanel.add(abLbl, "alignx center, cell 5 0");
		probeShiftLbl = new JLabel();
		probeShiftLbl.setFont(lblFont);
		hdrPanel.add(probeShiftLbl, "cell 1 1, span 2");
		hdrPanel.add(probeLbl, "grow, cell 3 1");
		return hdrPanel;
	}

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			@Override
			public void run() {
				try {
					// FIXME this is guaranteed to cause an NPE.
					// Should convert to a standard constructor that
					// handles default project creation.
					BlastFrame frame = new BlastFrame(null, null);
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

}
