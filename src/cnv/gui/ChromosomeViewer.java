package cnv.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.util.Vector;

import javax.swing.JPanel;

import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

public class ChromosomeViewer extends JPanel {
    public static final long serialVersionUID = 8L;

    public static final int  WIDTH_BUFFER     = 25;

    private int              chromosome       = 1;
    private int              startPosition    = 0;
    private int              stopPosition     = 0;
    private GeneTrack        track            = null; // Represents where the genes are
    float                    scalingFactor;

    public ChromosomeViewer(int chr, int start, int stop, GeneTrack gt) {
        chromosome = chr;
        startPosition = start;
        stopPosition = stop;
        track = gt;

        setMaximumSize(new Dimension(getWidth(), 20));
        setMinimumSize(new Dimension(getWidth(), 20));
    }

    public void setChr(int chr) {
        chromosome = chr;
    }

    public void setStart(int start) {
        startPosition = start;
    }

    public void setStop(int stop) {
        stopPosition = stop;
    }

    public void setTrack(GeneTrack track) {
        this.track = track;
    }

    public void updateView(int chr, int start, int stop) {
        chromosome = chr;
        startPosition = start;
        stopPosition = stop;

        int window = stopPosition - startPosition;
        scalingFactor = (float) getWidth() / window;
        // scalingFactor = 1;

        repaint();
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        GeneData[] genes;
        int[][] exons;
        Vector<Segment> v = new Vector<Segment>();
        Segment[] segs;
        int width, begin, end;
        int height = 10;

        // System.out.println("Start = " + startPosition + "  Stop = " + stopPosition + " Difference is " + (stopPosition - startPosition));
        if (track == null) {
            System.out.println("Track is null");
        }

        try {
            if (track != null && stopPosition - startPosition < 10000000) {
                genes = track.getBetween(chromosome, startPosition, stopPosition, 30);
                g.setColor(Color.BLACK);
                for (int i = 0; i < genes.length; i++) {
                    // begin = (int) (getX(genes[i].getStart()) * scalingFactor);
                    // end = (int) (getX(genes[i].getStop()) * scalingFactor);
                    begin = (int) (getX(genes[i].getStart()) * scalingFactor);
                    end = (int) (getX(genes[i].getStop()) * scalingFactor);
                    width = end - begin;
                    if (width < 1) {
                        width = 1;
                    }
                    g.drawRoundRect(begin, 0, width, height, 2, 2);
                    v.add(new Segment(begin, end));
                    exons = genes[i].getExonBoundaries();
                    for (int j = 0; j < exons.length; j++) {
                        begin = getX(exons[j][0]);
                        end = getX(exons[j][1]);
                        if (j == 0 || j == exons.length - 1) {
                            g.fillRoundRect(begin, 0, width, height, 2, 2);
                        } else {
                            g.fillRect(begin, 0, width, height);
                        }

                    }
                }
                Segment.mergeOverlapsAndSort(v);
                segs = Segment.toArray(v);
                g.setFont(new Font("Arial", 0, 14));

                for (int i = 0; i < genes.length; i++) {
                    begin = getX(genes[i].getStart());
                    width = g.getFontMetrics(g.getFont()).stringWidth(genes[i].getGeneName());
                    if (!Segment.overlapsAny(new Segment(begin - width - 5, begin - 1), segs)) {
                        g.drawString(genes[i].getGeneName(), begin, (height * 2) + 2);
                    }
                }

            }
        } catch (Exception ex) {
            System.out.println("Missing data in gene track");
        }
    }

    public int getX(int pos) {
        return (int) ((double) (pos - startPosition) / (double) (stopPosition - startPosition) * (double) (getWidth() - 2 * WIDTH_BUFFER)) + WIDTH_BUFFER;
    }
}
