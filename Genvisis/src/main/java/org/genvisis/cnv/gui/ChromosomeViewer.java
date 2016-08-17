package org.genvisis.cnv.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import javax.swing.JPanel;

import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.Segment;

public class ChromosomeViewer extends JPanel implements MouseMotionListener {
  public static final long serialVersionUID = 8L;

  public static final int WIDTH_BUFFER = 25;

  private int chromosome = 1;
  private int startPosition = 0;
  private int stopPosition = 0;
  private GeneTrack track = null; // Represents where the genes are
  float scalingFactor;
  private final ArrayList<GeneRectangle> geneRects;
  private GeneData[] oldGenes;
  private final int height = 10;
  private final int chromosomeRowStart = 0;
  private final int geneNameRowStart = (height * 2) + 2;
  private final int rulerRowStart = geneNameRowStart + 2;

  public ChromosomeViewer(int chr, int start, int stop, GeneTrack gt) {
    chromosome = chr;
    startPosition = start;
    stopPosition = stop;
    track = gt;
    geneRects = new ArrayList<GeneRectangle>();
    oldGenes = new GeneData[0];
    addMouseMotionListener(this);

    setMaximumSize(new Dimension(getWidth(), 45));
    setMinimumSize(new Dimension(getWidth(), 45));
  }

  /**
   * Calculate the position given the x coordinate
   * 
   * @param x
   * @return
   */
  public int getPos(int x) {
    return (int) (x / scalingFactor) + startPosition;
  }

  /**
   * Calculate the X coordinate of the given base position
   * 
   * @param pos
   * @return
   */
  public int getX(int pos) {
    return (int) ((pos - startPosition) * scalingFactor);
  }

  @Override
  public void mouseDragged(MouseEvent e) {
    // Nothing to do on mouse drag
  }

  /**
   * Creates a tooltip indicating which genes are included in the position being moused over
   */
  @Override
  public void mouseMoved(MouseEvent e) {
    // Create a mouseover tooltip with the name of any genes within this range
    setToolTipText(null);
    ArrayList<String> genes = new ArrayList<String>();

    // Work through the list of genes in this view and see if we're inside their bounding boxes
    for (GeneRectangle geneRect : geneRects) {
      if (geneRect.contains(e.getPoint())) {
        String name = geneRect.getName();
        if (!genes.contains(name)) {
          genes.add(name);
        }
      }
    }

    // Sometimes multiple genes are marked in the same location, show all of them
    String toolTipText = "";
    for (int i = 0; i < genes.size(); i++) {
      if (i == genes.size() - 1) {
        toolTipText += genes.get(i);
      } else {
        toolTipText += genes.get(i) + ", ";
      }
    }

    // Set the tooltip based on what we're over
    if (toolTipText.equals("")) {
      if (e.getY() > rulerRowStart) {
        // Set the tooltip to the current position if we're inside the ruler
        setToolTipText(Integer.toString(getPos(e.getX())));
      }
    } else {
      // Set the tooltip to the current gene
      setToolTipText(toolTipText);
    }
  }

  @Override
  public void paintComponent(Graphics g) {
    super.paintComponent(g);
    GeneData[] genes;
    int[][] exons;
    Vector<Segment> v = new Vector<Segment>();
    int width, begin, end;
    int lastNameEnd = 0;

    if (track == null) {
      System.out.println("Track is null");
    }

    try {
      if (track != null) {
        genes = track.getBetween(chromosome, startPosition, stopPosition, 30);
        if (!Arrays.equals(genes, oldGenes)) {
          // We've changed the list of genes, so update the rectangles list as we create them
          geneRects.clear();
          oldGenes = genes;
        }

        // Draw the row containing the chromosome view
        g.setColor(Color.BLACK);
        for (GeneData gene : genes) {
          begin = getX(gene.getStart());
          end = getX(gene.getStop());
          width = end - begin;
          if (width < 1) {
            width = 1;
          }

          // Maintain a list of the rectangles so we can check them on mouseover
          GeneRectangle rectangle =
              new GeneRectangle(gene.getGeneName(), new Rectangle(begin, 0, width, height));
          geneRects.add(rectangle);

          g.drawRoundRect(begin, chromosomeRowStart, width, height, 2, 2);
          v.add(new Segment(begin, end));
          exons = gene.getExonBoundaries();
          for (int j = 0; j < exons.length; j++) {

            int exBegin = getX(exons[j][0]);
            int exEnd = getX(exons[j][1]);
            int exWidth = exEnd - exBegin;
            if (j == 0 || j == exons.length - 1) {
              g.fillRoundRect(exBegin, chromosomeRowStart, exWidth, height, 2, 2);
            } else {
              g.fillRect(exBegin, chromosomeRowStart, exWidth, height);
            }

          }
        }
        Segment.mergeOverlapsAndSort(v);

        // Draw the row containing the gene names
        g.setFont(new Font("Arial", 0, 14));
        for (GeneData gene : genes) {
          begin = getX(gene.getStart());
          width = g.getFontMetrics(g.getFont()).stringWidth(gene.getGeneName());

          // Keep track of the end of the name we rendered so we don't overwrite them
          if (begin > lastNameEnd) {
            g.drawString(gene.getGeneName(), begin, geneNameRowStart);
            lastNameEnd = begin + width;
          }
        }

        // Draw the row containing the ruler
        g.setFont(new Font("Arial", 0, 12));
        int numTicks = 10;
        int tickBottom = rulerRowStart + height;
        int tickScale = (stopPosition - startPosition) / numTicks;

        g.fillRect(getX(startPosition), rulerRowStart, getWidth(), 2);

        for (int i = 0; i <= numTicks; i++) {
          int tickPos = startPosition + (i * tickScale);

          int posWidth = 0;
          switch (i) {
            case 0:
              // The first one can be drawn with the position left-justified with full-height ticks
              g.fillRect(getX(tickPos), rulerRowStart, 2, height);
              g.drawString(Integer.toString(tickPos), getX(tickPos), tickBottom + height);
              break;
            case 5:
              // The middle one should be drawn with the position centered with full-height ticks
              g.fillRect(getX(tickPos), rulerRowStart, 2, height);
              posWidth = g.getFontMetrics(g.getFont()).stringWidth(Integer.toString(tickPos));
              g.drawString(Integer.toString(tickPos), getX(tickPos) - (posWidth / 2),
                           tickBottom + height);
              break;
            case 10:
              // The last one should be drawn with the position right-justified with full-height
              // ticks
              g.fillRect(getX(tickPos) - 2, rulerRowStart, 2, height);
              posWidth = g.getFontMetrics(g.getFont()).stringWidth(Integer.toString(tickPos));
              g.drawString(Integer.toString(tickPos), getX(stopPosition) - posWidth,
                           tickBottom + height);
              break;
            default:
              // Make the intermediate ticks half-height
              g.fillRect(getX(tickPos), rulerRowStart, 2, (height / 2));
          }
        }
      }
    } catch (Exception ex) {
      System.out.println("Missing data in gene track");
    }
  }

  /**
   * Moves the chromosome viewer to the specified position and updates the view
   * 
   * @param chr
   * @param start
   * @param stop
   */
  public void updateView(int chr, int start, int stop) {
    chromosome = chr;
    startPosition = start;
    stopPosition = stop;

    int window = stopPosition - startPosition;
    scalingFactor = (float) getWidth() / window;

    repaint();
  }
}


/**
 * A class to keep track of bounding boxes for genes
 *
 * @author Michael Vieths
 *
 */
class GeneRectangle {
  private final String geneName;
  private final Rectangle geneRect;

  public GeneRectangle(String gene, Rectangle rect) {
    geneName = gene;
    geneRect = rect;
  }

  /**
   * Returns whether the given point falls in this gene
   * 
   * @param point
   * @return
   */
  public boolean contains(Point point) {
    return geneRect.contains(point);
  }

  /**
   * Returns the name of this gene
   * 
   * @return
   */
  public String getName() {
    return geneName;
  }
}
