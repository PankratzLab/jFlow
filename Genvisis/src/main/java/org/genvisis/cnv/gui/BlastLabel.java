package org.genvisis.cnv.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.JLabel;

import org.genvisis.cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.MarkerSeqAnnotation;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Fonts;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;

class ReferenceLabel extends JLabel {

  private static final long serialVersionUID = 1L;

  BlastAnnotation annot;

  public ReferenceLabel() {
    setFont(BlastLabel.LBL_FONT);
  }

  public void setAnnotation(BlastAnnotation annot) {
    this.annot = annot;
  }

  @Override
  protected void paintComponent(Graphics g) {
    FontMetrics fm = getFontMetrics(getFont());
    int h = fm.getAscent();
    g.setColor(Color.BLACK);
    int baseX = 0;

    String txt = getText();
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

    }

  }
}


public class BlastLabel extends JLabel {

  private static final long serialVersionUID = 1L;
  private static final Font BASE_FONT =
                                      (Fonts.SOURCE_CODE_PRO_REGULAR == null ? Font.decode(Font.MONOSPACED)
                                                                             : Fonts.SOURCE_CODE_PRO_REGULAR);

  public static void setFontSize(int size) {
    LBL_FONT = BASE_FONT.deriveFont((float) size);
  }

  public static Font LBL_FONT = BASE_FONT.deriveFont(15f);
  public static final int CHAR_PADDING = 2;
  static boolean expanded = false; // static to affect all
  static TreeSet<Integer> spaces = new TreeSet<Integer>();
  static TreeMap<Integer, Integer> spaceSets = new TreeMap<Integer, Integer>();

  protected BlastAnnotation myAnnotation;
  Segment fullSegment;
  private String seq;
  private final TreeSet<Integer> mySpaces = new TreeSet<Integer>();
  private final TreeMap<Integer, Integer> mySpaceSets = new TreeMap<Integer, Integer>();
  boolean positiveStrand = false;
  boolean oppositeStrand = false;
  boolean reverseSequence = false;
  ArrayList<CigarSeq> seqParts = new ArrayList<CigarSeq>();
  private final int alignmentCount;

  int getAlignment() {
    return alignmentCount;
  }

  Strand getStrand() {
    return myAnnotation.getStrand();
  }

  PROBE_TAG getAB() {
    return myAnnotation.getTag();
  }

  public BlastLabel(Project proj, MarkerSeqAnnotation ref, BlastAnnotation annot,
                    ReferenceGenome refGen) {
    super();
    myAnnotation = annot;
    positiveStrand = annot.getStrand() == Strand.POSITIVE;
    oppositeStrand = ref.getStrand() != annot.getStrand();
    reverseSequence = oppositeStrand/*
                                     * && ((ref.getTopBotRef() == TOP_BOT.PLUS &&
                                     * ref.getTopBotProbe() == TOP_BOT.PLUS) || (ref.getTopBotRef()
                                     * == TOP_BOT.BOT && ref.getTopBotProbe() == TOP_BOT.BOT) ||
                                     * (ref.getTopBotRef() == TOP_BOT.TOP && ref.getTopBotProbe() ==
                                     * TOP_BOT.TOP) || (ref.getTopBotRef() == TOP_BOT.MINUS &&
                                     * ref.getTopBotProbe() == TOP_BOT.MINUS))
                                     */;
    fullSegment = BlastFrame.BlastUtils.getSegmentForAnnotation(ref, annot);
    alignmentCount = BlastFrame.BlastUtils.countAlignment(annot);
    if (refGen != null) {
      String[] seqArr = refGen.getSequenceFor(fullSegment);
      if (seqArr != null) {
        seq = Array.toStr(seqArr, "");
      } else {
        // TODO set to probe seq, with alterations, otherwise causes NPE
        seq = ref.getSequence();
      }
    } else {
      seq = ref.getSequence();
      // TODO set to probe seq, with alterations, otherwise causes NPE
    }
    setFont(LBL_FONT);
    if (!positiveStrand) {
      seq = flipBases(seq);
    }
    parse();
    String seq = getSeqPartsAsString();
    setText(seq);
  }

  public static String flipBases(String seq) {
    StringBuilder sb = new StringBuilder(seq);
    for (int i = 0; i < sb.length(); i++) {
      if (sb.charAt(i) == 'A') {
        sb.setCharAt(i, 'T');
      } else if (sb.charAt(i) == 'T') {
        sb.setCharAt(i, 'A');
      } else if (sb.charAt(i) == 'C') {
        sb.setCharAt(i, 'G');
      } else if (sb.charAt(i) == 'G') {
        sb.setCharAt(i, 'C');
      }
    }
    return sb.toString();
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
        if (/* flipSequence && */!ciggie.getOperator().consumesReadBases()) {
          Integer key =
                      /* reverseSequence ? refSeq.getSequence().length() - strandInd + 1 : */ strandInd;
          if (spaceSets.containsKey(key)) {
            spaceSets.put(key, Math.max(spaceSets.get(key), ciggie.getLength()));
          } else {
            spaceSets.put(key, ciggie.getLength());
          }
          mySpaceSets.put(key, ciggie.getLength());
          for (int i = strandInd; i < strandInd + ciggie.getLength(); i++) {
            int ind = /* reverseSequence ? refSeq.getSequence().length() - i + 1 : */ i;
            spaces.add(ind);
            mySpaces.add(ind);
            // System.out.println(ind + " :: " + cig.toString());
          }
        }
        int stop = index + ciggie.getLength();
        String cigSeq = seq.substring(index, stop);
        seqParts.add(new CigarSeq(ciggie, cigSeq, index));
        index = stop;
      } else {
        seqParts.add(new CigarSeq(ciggie,
                                  Array.toStr(Array.stringArray(ciggie.getLength(), "."), ""),
                                  index));
      }
      strandInd += ciggie.getLength();
    }
  }

  // chr17:74,783,145-74,785,655
  // chr17:74,784,671-74,787,181
  // chr17:4,914,805-5,526,347 -- #37 of 4656
  // chr7:104,057,182-104,954,272 -- #63
  // chr19:33,876,905-33,881,272

  @Override
  protected void paintComponent(Graphics g) {
    FontMetrics fm = getFontMetrics(getFont());
    int h = fm.getAscent();
    g.setColor(Color.BLACK);
    int baseX = 0;
    int charInd = 0;
    int mySpacesCnt = 0;
    int addedSpaces = 0;
    for (int i =
               reverseSequence ? seqParts.size() - 1
                               : 0; reverseSequence ? i >= 0
                                                    : i < seqParts.size(); i +=
                                                                             reverseSequence ? -1
                                                                                             : 1) {
      CigarSeq cs = seqParts.get(i);
      boolean diff = cs.elem.getOperator() != CigarOperator.EQ;
      boolean read = cs.elem.getOperator().consumesReadBases();
      boolean ref = cs.elem.getOperator().consumesReferenceBases();
      if (diff) {
        if (read && ref) {
          // g.setColor(Color.RED.darker());
          g.setColor(Color.RED);
        } else if (read) {
          g.setColor(Color.GRAY);
        } else if (ref) {
          // g.setColor(Color.BLUE.darker());
          g.setColor(Color.BLUE.brighter());
        }
      } else {
        // g.setColor(Color.GREEN.darker());
        g.setColor(Color.BLACK);
      }
      boolean collapsed = !expanded && diff && !read && ref;
      if (collapsed) {
        g.drawString("]", baseX - (3 * CHAR_PADDING), h);
        g.drawString("[", baseX - (3 * CHAR_PADDING), h);
      } else {
        boolean strike = diff && read && !ref;
        int tempX = baseX;
        for (int c =
                   reverseSequence ? cs.elemSeq.length() - 1
                                   : 0; reverseSequence ? c >= 0
                                                        : c < cs.elemSeq.length(); c +=
                                                                                     reverseSequence ? -1
                                                                                                     : 1) {
          // if (expanded) {
          // if (mySpaces.contains(charInd - addedSpaces)) {
          // mySpacesCnt++;
          // } else {
          // if (spaces.contains(charInd - mySpacesCnt) && !mySpaces.contains(charInd -
          // mySpacesCnt)) {
          // Color col = g.getColor();
          // g.setColor(Color.BLACK);
          //
          // g.drawString("-", baseX, h);
          // baseX += CHAR_PADDING;
          // baseX += fm.charWidth('-');
          // tempX = baseX;
          //
          // charInd++;
          // addedSpaces++;
          // buildup++;
          // c += strandFlipped ? 1 : -1;
          //
          // g.setColor(col);
          // continue;
          // }
          // if (buildup > 0) {
          // for (int j = 0; j < buildup; j++) {
          //
          // g.drawString(cs.elemSeq.substring(c, c+1), baseX, h);
          // baseX += CHAR_PADDING; // char padding
          // baseX += fm.charWidth(cs.elemSeq.charAt(c));
          //
          // charInd++;
          // }
          // }
          // }
          // }


          if (expanded && !mySpaces.contains(charInd - addedSpaces)
              && (spaces.contains(charInd - mySpacesCnt)
                  && !mySpaces.contains(charInd - mySpacesCnt))) {
            Color col = g.getColor();
            g.setColor(Color.BLACK);

            g.drawString("-", baseX, h);
            baseX += CHAR_PADDING;
            baseX += fm.charWidth('-');
            tempX = baseX;

            charInd++;
            addedSpaces++;
            c += reverseSequence ? 1 : -1;

            g.setColor(col);
            continue;
          } else if (expanded && mySpaces.contains(charInd - addedSpaces)) {
            mySpacesCnt++;
          }

          g.drawString(cs.elemSeq.substring(c, c + 1), baseX, h);
          baseX += CHAR_PADDING; // char padding
          baseX += fm.charWidth(cs.elemSeq.charAt(c));

          charInd++;
          // if (expanded) {
          //
          // if (spaceSets.containsKey(Integer.valueOf(charInd - addedSpaces)) &&
          // !mySpaceSets.containsKey(Integer.valueOf(charInd - addedSpaces))) {
          // Color col = g.getColor();
          // g.setColor(Color.BLACK);
          // int cnt = spaceSets.get(Integer.valueOf(charInd - addedSpaces));
          // for (int k = 0; k < cnt; k++) {
          // g.drawString("-", baseX, h);
          // baseX += CHAR_PADDING;
          // baseX += fm.charWidth('-');
          // tempX = baseX;
          //
          //// charInd++;
          // addedSpaces++;
          //// c += strandFlipped ? 1 : -1;
          // }
          // g.setColor(col);
          // charInd += cnt;
          // } else if (mySpaceSets.containsKey(Integer.valueOf(charInd - addedSpaces))) {
          // mySpacesCnt += mySpaceSets.get(Integer.valueOf(charInd - addedSpaces));
          // }
          // }
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
      throw new RuntimeException("ERROR - Sequence {" + seq
                                 + "} does not match the length of the given CigarElement {"
                                 + elem.getLength() + elem.getOperator() + "}");
    }
    this.elem = elem;
    elemSeq = seq;
  }

  public String getDisplayString() {
    return "";
  }

}

