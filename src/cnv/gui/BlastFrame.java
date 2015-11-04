package cnv.gui;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;

import java.awt.BorderLayout;
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
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import common.Array;
import seq.manage.ReferenceGenome;
import cnv.annotation.MarkerBlastAnnotation;
import cnv.annotation.MarkerSeqAnnotation;
import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.filesys.Project;
import filesys.Segment;
import net.miginfocom.swing.MigLayout;

public class BlastFrame extends JFrame implements WindowFocusListener {

    public static class BlastUtils {
        
        public static ArrayList<BlastAnnotation> filterAnnotations(Project proj, List<BlastAnnotation> annotations, int alignFilter) {
            ArrayList<BlastAnnotation> filteredList = new ArrayList<BlastAnnotation>();
            for (BlastAnnotation annot : annotations) {
                if (countAlignment(annot) > alignFilter) {
                    filteredList.add(annot);
                }
            }
            return filteredList;
        }
        
        public static int countAlignment(BlastAnnotation annotation) {
            int align = 0;
            for (CigarElement ce : annotation.getCigar().getCigarElements()) {
                if (ce.getOperator() == CigarOperator.EQ) {
                    align += ce.getLength();
                }
            }
            return align;
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
        
    }

    private Project proj;
    private MarkerSeqAnnotation referenceAnnotation;
    private ArrayList<BlastAnnotation> annotations;
    private ReferenceGenome referenceGenome;
    private JLabel refLabel;
    private JLabel strandLbl;
    private JLabel probeLengthLbl;
    private JLabel locationLbl;
    private JLabel probeDescLbl;
    private JLabel probeLbl;
    
    private static final long serialVersionUID = 1L;
    private JPanel contentPane;
    private JCheckBox chckbxExpandBlastResults;
    private JButton btnClose;
    private JSpinner spnFontSize;
    private JScrollPane scrollPane;
    private JPanel blastPanel;
    private int rowCnt = 0;
    private JCheckBox chckbxSortByLocation;
    private JSeparator separator_1;
    private JCheckBox chckbxPinToFront;
    HashMap<BlastLabel, JLabel> alignCntLblMap = new HashMap<BlastLabel, JLabel>();
    private ArrayList<BlastLabel> bLabels;
    private ArrayList<BlastLabel> sorted;
    private JSpinner spinnerAlignmentLength;
    private JLabel lblAlignmentLength;
    private JSeparator separator_2;
    private MarkerBlastAnnotation blastResult;
    public int currentAlignFilter;
    private String nextBase;
    
    public void addBlastLabel(BlastLabel lbl) {
        JLabel locLbl = new JLabel();
        locLbl.setText(lbl.fullSegment.getUCSClocation());
        Font lblFont = Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
        locLbl.setFont(lblFont);
        this.blastPanel.add(locLbl, "cell 0 " + rowCnt);
        JLabel strandLbl = new JLabel();
        strandLbl.setText(lbl.getStrand().getEncoding());
        strandLbl.setFont(lblFont);
        alignCntLblMap.put(lbl, strandLbl);
        this.blastPanel.add(strandLbl, "cell 1 " + rowCnt);
        JLabel probeLengthLbl = new JLabel();
        probeLengthLbl.setText(lbl.getAlignment() + "");
        probeLengthLbl.setFont(lblFont);
        this.blastPanel.add(probeLengthLbl, "cell 2 " + rowCnt);
        this.blastPanel.add(lbl, "grow, cell 3 " + rowCnt);
        rowCnt++;
    }
    public void clearLabels() {
        this.blastPanel.removeAll();
    }
    
    public void windowGainedFocus(WindowEvent e) {}
    public void windowLostFocus(WindowEvent e) {
        if(e.getNewState() != WindowEvent.WINDOW_CLOSED) {
            if (chckbxPinToFront != null && chckbxPinToFront.isSelected()){
                setAlwaysOnTop(true);
            } else {
                setAlwaysOnTop(false);
            }
        }

    }
    
    /**
     * Create the frame.
     */
    public BlastFrame(Project proj) {
        this.proj = proj;
        setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        setBounds(100, 100, 1000, 800);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        contentPane.setLayout(new BorderLayout(0, 0));
        setContentPane(contentPane);

        blastPanel = new JPanel(new MigLayout("", "[200px][20px][20px][grow]", ""));
        scrollPane = new JScrollPane(blastPanel);
        
        contentPane.add(scrollPane, BorderLayout.CENTER);
        
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
                if (((javax.swing.JCheckBox)e.getSource()).isSelected()) {
                    for (BlastLabel bl : sorted) {
                        BlastFrame.this.addBlastLabel(bl);
                    }
                } else {
                    for (BlastLabel bl : bLabels) {
                        BlastFrame.this.addBlastLabel(bl);
                    }
                }
                BlastFrame.this.revalidate();
                BlastFrame.this.repaint();
            }
        });
        chckbxSortByLocation.setText("Sort By Location");
        chckbxSortByLocation.setVerticalAlignment(SwingConstants.BOTTOM);
        panel.add(chckbxSortByLocation, "cell 2 0");
        
        JSeparator separator = new JSeparator();
        separator.setOrientation(SwingConstants.VERTICAL);
        panel.add(separator, "cell 4 0,alignx left,growy");
        
        spnFontSize = new JSpinner(new SpinnerNumberModel(17, 5, 25, 1));
        spnFontSize.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                BlastLabel.setFontSize((Integer) spnFontSize.getModel().getValue());
                for (BlastLabel lbl : bLabels) {
                    lbl.updateFont();
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
        spinnerAlignmentLength.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                int value = ((SpinnerNumberModel) spinnerAlignmentLength.getModel()).getNumber().intValue();
                updateAnnotations(value - 1);
                updateLabels();
                BlastFrame.this.repaint();
            }
        });
        panel.add(spinnerAlignmentLength, "flowx,cell 7 0");
        
        chckbxPinToFront = new JCheckBox("Pin to Front");
        panel.add(chckbxPinToFront, "cell 9 0");
        
        btnClose = new JButton();
        btnClose.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                BlastFrame.this.setVisible(false);
//                BlastFrame.this.dispose(); // removed for state
            }
        });
        btnClose.setText("Close");
        panel.add(btnClose, "cell 10 0");
        
        lblAlignmentLength = new JLabel("Alignment Length");
        panel.add(lblAlignmentLength, "cell 7 0");
        

        refLabel = new ReferenceLabel();
        probeLbl = new ReferenceLabel();
//        indexLabel.setFont(BlastLabel.LBL_FONT);
//        indexLabel.setText(referenceAnnotation.getSequence());
        BlastLabel.setRefLabel(refLabel);
//        JPanel panel = new JPanel(new GridLayout(2, 1));
//        panel.add(indexLabel);
//        panel.add(refLabel);
//        jsp.setColumnHeaderView(panel);
        JPanel hdrPanel = getHeaderPanel();
        scrollPane.setColumnHeaderView(hdrPanel);
        
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
//        frame.pack();
        
        addWindowFocusListener(this);
    }
    
    public void setAnnotations(MarkerBlastAnnotation blastResult, ReferenceGenome refGen) {
        this.referenceAnnotation = blastResult.getMarkerSeqAnnotation();
        this.blastResult = blastResult;
        this.referenceGenome = refGen;
        double filter = proj.BLAST_PROPORTION_MATCH_FILTER.getValue();
        int probe = proj.ARRAY_TYPE.getValue().getProbeLength();
        int alignFilter = (int) (filter * probe);
        this.updateAnnotations(alignFilter);
        this.updateLabels();
    }
    
    private void updateAnnotations(int alignFilter) {
        currentAlignFilter = alignFilter;
        this.annotations = BlastFrame.BlastUtils.filterAnnotations(proj, blastResult.getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS, proj.getLog()), alignFilter);
    }
    
    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    BlastFrame frame = new BlastFrame(null);
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
    
    public void updateLabels() {
        BlastLabel.spaces = new TreeSet<Integer>();
        BlastLabel.spaceSets = new TreeMap<Integer, Integer>();
        bLabels = new ArrayList<BlastLabel>();
        sorted = new ArrayList<BlastLabel>();
        for (BlastAnnotation annot : annotations) {
            BlastLabel lbl = new BlastLabel(referenceAnnotation, annot, referenceGenome);
    //            if (lbl.strandFlipped) continue;
            bLabels.add(lbl);
            sorted.add(lbl);
        }
        Collections.sort(sorted, new Comparator<BlastLabel>() {
            @Override
            public int compare(BlastLabel arg0, BlastLabel arg1) {
                int res1 = (new Byte(arg0.fullSegment.getChr())).compareTo(new Byte(arg1.fullSegment.getChr()));
                if (res1 != 0) return res1;
                res1 = (new Integer(arg0.fullSegment.getStart())).compareTo(new Integer(arg1.fullSegment.getStart()));
                if (res1 != 0) return res1;
                res1 = (new Integer(arg0.fullSegment.getStop())).compareTo(new Integer(arg1.fullSegment.getStop()));
                return res1;
            }
        });
        
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                BlastFrame.this.clearLabels();
                for (BlastLabel lbl : bLabels) {
                    BlastFrame.this.addBlastLabel(lbl);
                }
                boolean posStrand = BlastFrame.this.referenceAnnotation.getStrand() == Strand.POSITIVE;
                refLabel.setText(!posStrand ? new StringBuilder(BlastFrame.this.referenceAnnotation.getSequence()).reverse().toString() : BlastFrame.this.referenceAnnotation.getSequence());
                strandLbl.setText(BlastFrame.this.referenceAnnotation.getStrand().getEncoding());
                int len = proj.ARRAY_TYPE.getValue().getProbeLength();
                probeLengthLbl.setText(len + "");
                Segment seg = BlastFrame.this.referenceAnnotation.getSeg();
                int start = (posStrand ? seg.getStart() - len : seg.getStart());
                int stop = (posStrand ?  seg.getStart() : seg.getStart() + len);
                locationLbl.setText(seg.getChromosomeUCSC() + ":" + start + "-" + stop);
                String[] gen = referenceGenome.getSequenceFor(new Segment(seg.getChr(), start, stop));
                String[] act = Array.subArray(gen, posStrand ? 0 : 1, posStrand ? gen.length - 1 : gen.length);
                nextBase = posStrand ? gen[gen.length - 1] : gen[0];
                probeLbl.setText(Array.toStr(act, ""));
                BlastFrame.this.revalidate();
                BlastFrame.this.repaint();
            }
        });        
    }

    private JPanel getHeaderPanel() {
        JPanel hdrPanel = new JPanel(new MigLayout("", "[200px][20px][20px][grow]", "[][]")); // TODO
        hdrPanel.setBorder(null);
        locationLbl = new JLabel();
        Font lblFont = Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
        locationLbl.setFont(lblFont);
        hdrPanel.add(locationLbl, "cell 0 0");
        probeDescLbl = new JLabel("<reference genome read:>");
        probeDescLbl.setFont(lblFont);
        hdrPanel.add(probeDescLbl, "cell 0 1");
        strandLbl = new JLabel(referenceAnnotation == null ? "" : referenceAnnotation.getStrand().getEncoding());
        strandLbl.setFont(lblFont);
        hdrPanel.add(strandLbl, "cell 1 0");
        probeLengthLbl = new JLabel(proj.ARRAY_TYPE.getValue().getProbeLength() + "");
        probeLengthLbl.setFont(lblFont);
        hdrPanel.add(probeLengthLbl, "cell 2 0");
        hdrPanel.add(refLabel, "grow, cell 3 0");
        hdrPanel.add(probeLbl, "grow, cell 3 1");
        return hdrPanel;
    }
    
}
