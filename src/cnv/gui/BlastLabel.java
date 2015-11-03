package cnv.gui;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.JLabel;

import seq.manage.ReferenceGenome;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.annotation.MarkerSeqAnnotation;
import common.Array;
import common.Fonts;
import filesys.Segment;

class ReferenceLabel extends JLabel {
    
    private static final long serialVersionUID = 1L;
    
    public ReferenceLabel() {
        setFont(BlastLabel.LBL_FONT);
    }
    
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
}

public class BlastLabel extends JLabel {

    private static final long serialVersionUID = 1L;
    private static final Font BASE_FONT = (Fonts.SOURCE_CODE_PRO_LIGHT == null ? Font.decode(Font.MONOSPACED) : Fonts.SOURCE_CODE_PRO_LIGHT);
    
    public static void setFontSize(int size) {
        LBL_FONT = BASE_FONT.deriveFont((float) size);
        if (refLbl != null) {
            refLbl.setFont(LBL_FONT);
        }
    }
    
    public static void setRefLabel(JLabel lbl) {
        refLbl = lbl;
    }
    
    public static Font LBL_FONT = BASE_FONT.deriveFont(17f);
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
    private int alignmentCount;
    
    int getAlignment() {
        return alignmentCount;
    }
    
    Strand getStrand() {
        return myAnnotation.getStrand();
    }
    
    public BlastLabel(MarkerSeqAnnotation ref, BlastAnnotation annot, ReferenceGenome refGen) {
        super();
        this.refSeq = ref;
        this.myAnnotation = annot;
        this.strandFlipped = ref.getStrand() != annot.getStrand();
        this.fullSegment = BlastFrame.BlastUtils.getSegmentForAnnotation(ref, annot);
        this.alignmentCount = BlastFrame.BlastUtils.countAlignment(annot);
        if (refGen != null) {
            String[] seqArr = refGen.getSequenceFor(this.fullSegment);
            if (seqArr != null) {
                this.seq = Array.toStr(seqArr, "");
            } else {
                // TODO set to probe seq, with alterations, otherwise causes NPE
            }
        } else {
            // TODO set to probe seq, with alterations, otherwise causes NPE
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
//    chr17:4,914,805-5,526,347 -- #37 of 4656
    
    @Override
    protected void paintComponent(Graphics g) {
        FontMetrics fm = this.getFontMetrics(this.getFont()); 
        int h = fm.getAscent();
        g.setColor(Color.BLACK);
        int baseX = 0;
        int charInd = 0;
        int mySpacesCnt = 0;
        int addedSpaces = 0;
        int buildup = 0;
        
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
//                    if (expanded) {
//                        if (mySpaces.contains(charInd - addedSpaces)) {
//                            mySpacesCnt++;
//                        } else {
//                            if (spaces.contains(charInd - mySpacesCnt) && !mySpaces.contains(charInd - mySpacesCnt)) {
//                                Color col = g.getColor();
//                                g.setColor(Color.BLACK);
//
//                                g.drawString("-", baseX, h);
//                                baseX += CHAR_PADDING;
//                                baseX += fm.charWidth('-');
//                                tempX = baseX;
//
//                                charInd++;
//                                addedSpaces++;
//                                buildup++;
//                                c += strandFlipped ? 1 : -1;
//
//                                g.setColor(col);
//                                continue;
//                            }
//                            if (buildup > 0) {
//                                for (int j = 0; j < buildup; j++) {
//
//                                    g.drawString(cs.elemSeq.substring(c, c+1), baseX, h);
//                                    baseX += CHAR_PADDING; // char padding
//                                    baseX += fm.charWidth(cs.elemSeq.charAt(c));
//
//                                    charInd++;
//                                }
//                            }
//                        }
//                    }
                    
                    
                    if (expanded && !mySpaces.contains(charInd - addedSpaces) && (spaces.contains(charInd - mySpacesCnt) && !mySpaces.contains(charInd - mySpacesCnt))) {
                        Color col = g.getColor();
                        g.setColor(Color.BLACK);
                        
                        g.drawString("-", baseX, h);
                        baseX += CHAR_PADDING;
                        baseX += fm.charWidth('-');
                        tempX = baseX;
                        
                        charInd++;
                        addedSpaces++;
                        c += strandFlipped ? 1 : -1;
                        
                        g.setColor(col);
                        continue;
                    } else if (expanded && mySpaces.contains(charInd - addedSpaces)) {
                        mySpacesCnt++;
                    }
                    
                    g.drawString(cs.elemSeq.substring(c, c+1), baseX, h);
                    baseX += CHAR_PADDING; // char padding
                    baseX += fm.charWidth(cs.elemSeq.charAt(c));

                    charInd++;
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











