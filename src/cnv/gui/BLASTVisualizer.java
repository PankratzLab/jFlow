package cnv.gui;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.annotation.Strand;

import java.awt.Font;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JLabel;

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
    
    public BLASTVisualizer(Project proj, MarkerSeqAnnotation seq, ArrayList<BlastAnnotation> annots, ReferenceGenome refGen) {
        this.proj = proj;
        this.referenceAnnotation = seq;
        this.annotations = annots;
        this.referenceGenome = refGen;
    }
    
//    chr17:74,783,145-74,785,655
    public void run() {
        
        System.out.println();
        String direc = referenceAnnotation.getStrand() == Strand.NEGATIVE ? "-" : "+";
        System.out.println(referenceAnnotation.getSequence() + "\t]" + direc + "[\t");
        
        for (BlastAnnotation annot : annotations) {
            BlastAnnotationVisualizer.processAnnotation(referenceAnnotation, annot, referenceGenome);
        }
        
        
    }
    
}

class BlastAnnotationVisualizer {
    static int count = 5;
    public static void processAnnotation(MarkerSeqAnnotation ref, BlastAnnotation annotation, ReferenceGenome refGen) {
        
        if (refGen != null) {
            Segment seg = getSegmentForAnnotation(ref, annotation);
            String[] seqArr = refGen.getSequenceFor(seg);
            if (seqArr != null) {
                String seq = Array.toStr(seqArr, "");
                if (ref.getStrand() != annotation.getStrand()) { 
                    seq = (new StringBuilder(seq)).reverse().toString();
                }
                String direc = annotation.getStrand() == Strand.NEGATIVE ? "-" : "+";
                System.out.println(seq  + "\t]" + direc + "[\t" + annotation.getCigar().toString() + "\t" + seg.getUCSClocation());
            }
        }
        
        if (count-- > 0) {
            JFrame frame = new JFrame();
            BlastLabel label = new BlastLabel(ref, annotation, refGen);
            frame.add(label);
            frame.pack();
            frame.setVisible(true);
        }
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
    
    private BlastAnnotation myAnnotation;
    private MarkerSeqAnnotation refSeq;
    private Segment fullSegment;
    private String seq;
    boolean strandFlipped = false;
    ArrayList<CigarSeq> seqParts = new ArrayList<CigarSeq>();
    
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
        this.setFont(Font.decode(Font.MONOSPACED));
        parse();
        setText(getSeqPartsAsString());
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
        for (CigarElement ciggie : cig.getCigarElements()) {
            if (ciggie.getOperator().consumesReferenceBases()) {
                int stop = index + ciggie.getLength();
                String cigSeq = seq.substring(index, stop);
                seqParts.add(new CigarSeq(ciggie, cigSeq));
                index = stop;
            } else {
                seqParts.add(new CigarSeq(ciggie, Array.toStr(Array.stringArray(ciggie.getLength(), "."), "")));
            }
        }
    }
    
    
    
}

class CigarSeq {
    CigarElement elem;
    String elemSeq;

    public CigarSeq(CigarElement elem, String seq) {
        if (elem.getLength() != seq.length()) {
            throw new RuntimeException("ERROR - Sequence {" + seq + "} does not match the length of the given CigarElement {" + elem.getLength() + elem.getOperator() + "}");
        }
        this.elem = elem;
        this.elemSeq = seq;
    }
    
}











