package cnv.gui;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.AbstractAction;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.ScrollPaneConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import net.miginfocom.swing.MigLayout;
import common.Array;
import seq.manage.ReferenceGenome;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.annotation.MarkerSeqAnnotation;
import cnv.filesys.Project;
import filesys.Segment;

public class BLASTVisualizer {
    
    private Project proj;
    private MarkerSeqAnnotation referenceAnnotation;
    private ArrayList<BlastAnnotation> annotations;
    private ReferenceGenome referenceGenome;
    private JLabel refLabel;
    private ArrayList<BlastLabel> bLabels;
    private ArrayList<BlastLabel> sorted;
    
    public BLASTVisualizer(Project proj, MarkerSeqAnnotation seq, ArrayList<BlastAnnotation> annots, ReferenceGenome refGen) {
        this.proj = proj;
        this.referenceAnnotation = seq;
        this.annotations = annots;
        this.referenceGenome = refGen;
    }

//  chr17:74,783,145-74,785,655
//  chr17:74,784,671-74,787,181
    
    public void run() {
        
        System.out.println();
        String direc = referenceAnnotation.getStrand() == Strand.NEGATIVE ? "-" : "+";
        System.out.println(referenceAnnotation.getSequence() + "\t]" + direc + "[\t");
        
        final BlastFrame frame = new BlastFrame(BlastLabel.expanded);
        
        bLabels = new ArrayList<BlastLabel>();
        sorted = new ArrayList<BlastLabel>();
        for (BlastAnnotation annot : annotations) {
            BlastLabel lbl = BlastAnnotationVisualizer.processAnnotation(referenceAnnotation, annot, referenceGenome);
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
        
        refLabel = new JLabel() {
            private static final long serialVersionUID = 1L;

            @Override
            protected void paintComponent(Graphics g) {
                FontMetrics fm = this.getFontMetrics(this.getFont());
                int h = fm.getAscent();
                g.setColor(Color.BLACK);
                int baseX = 0;

                String txt = this.getText();
                int charInd = 0;
                int added = 0;
                for (int c = 0; c < txt.length(); c++) {

                    if (BlastLabel.expanded) {

                        if (BlastLabel.spaceSets.containsKey(charInd) && added == 0) {
                            int len = BlastLabel.spaceSets.get(charInd);
//                            System.out.println(len);
                            for (int i = 0; i < len; i++) {
                                g.drawString("-", baseX, h);
                                baseX += BlastLabel.CHAR_PADDING; // char padding
                                baseX += fm.charWidth('-');
                                charInd++;
                            }
                            added += len;
                            c--;
                            continue;
                        }

                        g.drawString(txt.substring(c, c + 1), baseX, h);
                        baseX += BlastLabel.CHAR_PADDING; // char padding
                        baseX += fm.charWidth(txt.charAt(c));
                        if (added > 0) {
                            added--;
                        } else {
                            charInd++;
                        }

                    } else {
                        g.drawString(txt.substring(c, c + 1), baseX, h);
                        baseX += BlastLabel.CHAR_PADDING; // char padding
                        baseX += fm.charWidth(txt.charAt(c));
                    }

                    // if (BlastLabel.expanded && BlastLabel.spacescontains(charInd)) {
                    // g.drawString("-", baseX, h);
                    // baseX += BlastLabel.CHAR_PADDING; // char padding
                    // baseX += fm.charWidth('-');
                    // charInd++;
                    // c--;
                    // continue;
                    // }
                    //
                    //
                    // g.drawString(txt.substring(c, c+1), baseX, h);
                    // baseX += BlastLabel.CHAR_PADDING; // char padding
                    // baseX += fm.charWidth(txt.charAt(c));
                    // charInd++;

                }

            }
        };
//        indexLabel.setFont(BlastLabel.LBL_FONT);
//        indexLabel.setText(referenceAnnotation.getSequence());
        refLabel.setFont(BlastLabel.LBL_FONT);
        BlastLabel.setRefLabel(refLabel);
        refLabel.setText(referenceAnnotation.getSequence());
//        JPanel panel = new JPanel(new GridLayout(2, 1));
//        panel.add(indexLabel);
//        panel.add(refLabel);
//        jsp.setColumnHeaderView(panel);
        JPanel hdrPanel = new JPanel(new MigLayout("", "[200px][20px][grow]", ""));
        hdrPanel.setBorder(null);
        JLabel locLbl = new JLabel();
        locLbl.setText("<reference location>");
        Font lblFont = Font.decode(Font.MONOSPACED).deriveFont(Font.PLAIN, 12);
        locLbl.setFont(lblFont);
        hdrPanel.add(locLbl, "cell 0 0");
        JLabel strandLbl = new JLabel();
        strandLbl.setText("]" + referenceAnnotation.getStrand().getEncoding() + "[");
        strandLbl.setFont(lblFont);
        hdrPanel.add(strandLbl, "cell 1 0");
        hdrPanel.add(refLabel, "grow, cell 2 0");
        frame.getScrollPane().setColumnHeaderView(hdrPanel);
        
        for (BlastLabel lbl : bLabels) {
//            if (lbl.strandFlipped) {
                frame.addBlastLabel(lbl);
//            }
        }
        
        frame.setSpinnerAction(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                BlastLabel.setFontSize((Integer) frame.getSpinner().getModel().getValue());
                for (BlastLabel lbl : bLabels) {
                    lbl.updateFont();
                }
                frame.repaint();
            }
        });
        frame.setExpansionCheckBoxAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                BlastLabel.expanded = !BlastLabel.expanded;
                frame.repaint();
            }
        });
        frame.setSortCheckBoxAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                frame.clearLabels();
                if (((javax.swing.JCheckBox)e.getSource()).isSelected()) {
                    for (BlastLabel bl : sorted) {
                        frame.addBlastLabel(bl);
                    }
                } else {
                    for (BlastLabel bl : bLabels) {
                        frame.addBlastLabel(bl);
                    }
                }
                frame.revalidate();
                frame.repaint();
            }
        });
        
        frame.getScrollPane().setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        frame.getScrollPane().setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
//        frame.add(jsp);
        frame.pack();
        frame.setBounds(100, 100, 900, 800);
        frame.setVisible(true);
        
    }
    
}

class BlastAnnotationVisualizer {
    
    public static BlastLabel processAnnotation(MarkerSeqAnnotation ref, BlastAnnotation annotation, ReferenceGenome refGen) {
        
        if (refGen != null) {
            Segment seg = getSegmentForAnnotation(ref, annotation);
            String[] seqArr = refGen.getSequenceFor(seg);
            if (seqArr != null) {
                String seq = Array.toStr(seqArr, "");
                if (ref.getStrand() != annotation.getStrand()) { 
                    seq = (new StringBuilder(seq)).reverse().toString();
                }
                String direc = annotation.getStrand() == Strand.NEGATIVE ? "-" : "+";
//                if (annotation.getStrand() == ref.getStrand()) {
//                    System.out.println(seq  + "\t]" + direc + "[\t" + annotation.getCigar().toString() + "\t" + seg.getUCSClocation());
//                }
            }
        }
        
        BlastLabel label = new BlastLabel(ref, annotation, refGen);
        return label;
    }
    
    protected static Segment getSegmentForAnnotation(MarkerSeqAnnotation ref, BlastAnnotation annot) {
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



class BlastLabel extends JLabel {
    
    private static final long serialVersionUID = 1L;

    public static void setFontSize(int size) {
        LBL_FONT = Font.decode(Font.MONOSPACED).deriveFont((float) size);
        if (refLbl != null) {
            refLbl.setFont(LBL_FONT);
        }
    }
    
    public static void setRefLabel(JLabel lbl) {
        refLbl = lbl;
    }
    
    public static Font LBL_FONT = Font.decode(Font.MONOSPACED).deriveFont(17f);
    public static final int CHAR_PADDING = 2;
    private static JLabel refLbl = null;
    static boolean expanded = true; // static to affect all
    static TreeSet<Integer> spaces = new TreeSet<Integer>();
    static TreeMap<Integer, Integer> spaceSets = new TreeMap<Integer, Integer>();
    
    private BlastAnnotation myAnnotation;
    private MarkerSeqAnnotation refSeq;
    Segment fullSegment;
    private String seq;
    private TreeSet<Integer> mySpaces = new TreeSet<Integer>();
    private TreeMap<Integer, Integer> mySpaceSets = new TreeMap<Integer, Integer>();
    boolean strandFlipped = false;
    ArrayList<CigarSeq> seqParts = new ArrayList<CigarSeq>();
    
    Strand getStrand() {
        return myAnnotation.getStrand();
    }
    
    public BlastLabel(MarkerSeqAnnotation ref, BlastAnnotation annot, ReferenceGenome refGen) {
        super();
        this.refSeq = ref;
        this.myAnnotation = annot;
        this.strandFlipped = ref.getStrand() != annot.getStrand();
        this.fullSegment = BlastAnnotationVisualizer.getSegmentForAnnotation(ref, annot);
        if (refGen != null) {
            String[] seqArr = refGen.getSequenceFor(this.fullSegment);
            if (seqArr != null) {
                this.seq = Array.toStr(seqArr, "");
            } else {
                // set to probe seq, with alterations
            }
        } else {
            // set to probe seq, with alterations
        }
        this.setFont(LBL_FONT);
        parse();
        setText(getSeqPartsAsString());
    }
    
    void updateFont() {
        this.setFont(LBL_FONT);
    }
    
    private String getSeqPartsAsString() {
        StringBuilder sb = new StringBuilder();
        for (CigarSeq cs : seqParts) {
            sb.append(cs.elemSeq);
        }
        return sb.toString();
    }
    
    private void parse() {
        Cigar cig = myAnnotation.getCigar();
        int index = 0;
        int strandInd = 0;
        for (CigarElement ciggie : cig.getCigarElements()) {
            if (ciggie.getOperator().consumesReferenceBases()) {
                if (/*strandFlipped && */!ciggie.getOperator().consumesReadBases()) {
                    Integer key = strandFlipped ? refSeq.getSequence().length() - strandInd + 1 : strandInd;
                    if (spaceSets.containsKey(key)) {
                        spaceSets.put(key, Math.max(spaceSets.get(key), ciggie.getLength()));
                    } else {
                        spaceSets.put(key, ciggie.getLength());
                    }
                    mySpaceSets.put(key, ciggie.getLength());
                    for (int i = strandInd; i < strandInd + ciggie.getLength(); i++) {
                        int ind = strandFlipped ? refSeq.getSequence().length() - i + 1 : i;
                        spaces.add(ind);
                        mySpaces.add(ind);
//                        System.out.println(ind + " :: " + cig.toString());
                    }
                }
                int stop = index + ciggie.getLength();
                String cigSeq = seq.substring(index, stop);
                seqParts.add(new CigarSeq(ciggie, cigSeq, index));
                index = stop;
            } else {
                seqParts.add(new CigarSeq(ciggie, Array.toStr(Array.stringArray(ciggie.getLength(), "."), ""), index));
            }
            strandInd += ciggie.getLength();
        }
    }
    
//    chr17:74,783,145-74,785,655
//    chr17:74,784,671-74,787,181
    
    @Override
    protected void paintComponent(Graphics g) {
        FontMetrics fm = this.getFontMetrics(this.getFont()); 
        int h = fm.getAscent();
        g.setColor(Color.BLACK);
        int baseX = 0;
        int charInd = 0;
        int mySpacesCnt = 0;
        int addedSpaces = 0;
        
        for (int i = strandFlipped ? seqParts.size() - 1 : 0; strandFlipped ? i >= 0 : i < seqParts.size(); i += strandFlipped ? -1 : 1) {
            CigarSeq cs = seqParts.get(i);
            boolean diff = cs.elem.getOperator() != CigarOperator.EQ;
            boolean read = cs.elem.getOperator().consumesReadBases();
            boolean ref = cs.elem.getOperator().consumesReferenceBases();
            if (diff) {
                if (read && ref) {
                    g.setColor(Color.RED.darker());
                } else if (read) {
                    g.setColor(Color.GRAY);
                } else if (ref) {
                    g.setColor(Color.BLUE.darker());
                }
            } else {
                g.setColor(Color.GREEN.darker());
            }
            boolean collapsed = !expanded && diff && !read && ref;
            if (collapsed) {
                g.drawString("]", baseX - (3 * CHAR_PADDING), h);
                g.drawString("[", baseX - (3 * CHAR_PADDING), h);
            } else {
                boolean strike = diff && read && !ref;
                int tempX = baseX;
                for (int c = strandFlipped ? cs.elemSeq.length() - 1 : 0; strandFlipped ? c >= 0 : c < cs.elemSeq.length(); c += strandFlipped ? -1 : 1) {
//                    if (expanded && !mySpaces.contains(charInd - addedSpaces) && (spaces.contains(charInd - mySpacesCnt) && !mySpaces.contains(charInd - mySpacesCnt))) {
//                        Color col = g.getColor();
//                        g.setColor(Color.BLACK);
//                        
//                        g.drawString("-", baseX, h);
//                        baseX += CHAR_PADDING;
//                        baseX += fm.charWidth('-');
//                        tempX = baseX;
//                        
//                        charInd++;
//                        addedSpaces++;
//                        c += strandFlipped ? 1 : -1;
//                        
//                        g.setColor(col);
//                        continue;
//                    } else if (expanded && mySpaces.contains(charInd - addedSpaces)) {
//                        mySpacesCnt++;
//                    }
                    
//                    if (expanded) {
//                        
//                        if (spaceSets.containsKey(Integer.valueOf(charInd - addedSpaces)) && !mySpaceSets.containsKey(Integer.valueOf(charInd - addedSpaces))) {
//                            Color col = g.getColor();
//                            g.setColor(Color.BLACK);
//                            int cnt = spaceSets.get(Integer.valueOf(charInd - addedSpaces));
//                            for (int k = 0; k < cnt; k++) {
//                                g.drawString("-", baseX, h);
//                                baseX += CHAR_PADDING;
//                                baseX += fm.charWidth('-');
//                                tempX = baseX;
//    
////                                charInd++;
//                                addedSpaces++;
////                                c += strandFlipped ? 1 : -1;
//                            }
//                            g.setColor(col);
//                            charInd += cnt;
//                        } else if (mySpaceSets.containsKey(Integer.valueOf(charInd - addedSpaces))) {
//                            mySpacesCnt += mySpaceSets.get(Integer.valueOf(charInd - addedSpaces));
//                        }
//                    }
//                    
                    
                    
                    g.drawString(cs.elemSeq.substring(c, c+1), baseX, h);
                    baseX += CHAR_PADDING; // char padding
                    baseX += fm.charWidth(cs.elemSeq.charAt(c));

                    charInd++;
                }
                if (strike) {
                    g.setColor(Color.GRAY.darker());
                    g.drawLine(tempX, 2 * (h / 3), baseX - CHAR_PADDING, 2 * (h / 3));
                }
            }
        }
        
    }
    
}

class CigarSeq {
    int elemInd;
    CigarElement elem;
    String elemSeq;

    public CigarSeq(CigarElement elem, String seq, int index) {
        if (elem.getLength() != seq.length()) {
            throw new RuntimeException("ERROR - Sequence {" + seq + "} does not match the length of the given CigarElement {" + elem.getLength() + elem.getOperator() + "}");
        }
        this.elem = elem;
        this.elemSeq = seq;
    }
    
    public String getDisplayString() {
        return "";
    }
    
}











