package org.genvisis.one.ben.imagetag;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.HashMap;
import javax.imageio.ImageIO;
import javax.swing.UIManager;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Images;

public class AnnotatedImage {

  public static class Annotation {

    private static final HashMap<String, Character> mnemonicMap = new HashMap<>();
    final String annotation;
    final char mnemonic;

    public Annotation(String ann, String mnem) {
      this.annotation = ann;
      this.mnemonic = mnem.charAt(0);
      mnemonicMap.put(this.annotation, this.mnemonic);
    }

    public Annotation(String ann) {
      this.annotation = ann;
      if (mnemonicMap.containsKey(this.annotation)) {
        this.mnemonic = mnemonicMap.get(this.annotation);
      } else {
        int ind = 0;
        Character c = this.annotation.toUpperCase().charAt(ind);
        while (mnemonicMap.containsValue(c)) {
          ind++;
          c = this.annotation.toUpperCase().charAt(ind);
        }
        if (ind == this.annotation.length()) {
          System.err.println("Error - all possible mnemonic characters already used for annotation {"
                             + this.annotation + "}.  Using alphanumerics instead.");
          String alphanum = "abcdefghijklmnopqrstuvwxyz0123456789";
          ind = 0;
          while (mnemonicMap.containsValue(alphanum.charAt(ind))) {
            ind++;
            if (ind == alphanum.length()) {
              System.err.println("Error - ran out of alphanumeric mnemonics too!");
              break;
            }
          }
          if (ind < alphanum.length()) {
            this.mnemonic = alphanum.charAt(ind);
            mnemonicMap.put(this.annotation, this.mnemonic);
          } else {
            this.mnemonic = '-';
          }
        } else {
          mnemonicMap.put(this.annotation, c);
          this.mnemonic = c;
        }
      }
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((annotation == null) ? 0 : annotation.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      Annotation other = (Annotation) obj;
      if (annotation == null) {
        if (other.annotation != null) return false;
      } else if (!annotation.equals(other.annotation)) return false;
      return true;
    }

  }

  private String name;
  private String imageFile;
  private SoftReference<BufferedImage> image;
  private ArrayList<AnnotatedImage.Annotation> annots;
  private final boolean isRoot;
  private boolean missing = false;

  private static final String DELIM = "|";

  public String exportToString() {
    StringBuilder sb = new StringBuilder();
    sb.append(name).append(DELIM).append(imageFile == null ? "" : imageFile);
    for (AnnotatedImage.Annotation a : annots) {
      sb.append(DELIM).append(a.annotation);
    }
    return sb.toString();
  }

  public AnnotatedImage(String ident, boolean isRoot) {
    this.name = ident;
    this.isRoot = isRoot;
    annots = new ArrayList<>();
  }

  public void setMissing(boolean miss) {
    missing = miss;
  }

  @Override
  public String toString() {
    return name + (missing ? " - MISSING" : "");
  }

  public BufferedImage getImage() {
    if (image == null) {
      createImage();
    }
    BufferedImage bi;
    while ((bi = image.get()) == null) {
      createImage();
    }
    return bi;
  }

  private void createImage() {
    if (image == null) {
      if (imageFile != null) {
        if (!imageFile.contains(";")) {
          try {
            image = new SoftReference<>(ImageIO.read(new File(imageFile)));
          } catch (IOException e) {
            e.printStackTrace();
            image = new SoftReference<>(createIOExceptionImage(e));
          }
        } else {
          String[] images = imageFile.split(";");
          image = new SoftReference<>(Images.stitchImages(images, Color.WHITE, false, false));
        }
      } else {
        if (imageFile == null) {
          image = new SoftReference<>(createNoFileImage());
        } else if (!Files.exists(imageFile)) {
          image = new SoftReference<>(createMissingFileImage(imageFile));
        }
      }
    }
  }

  private static BufferedImage createImage(String msg, String msg2) {
    BufferedImage bi = new BufferedImage(800, 600, BufferedImage.TYPE_INT_RGB);
    Graphics g = bi.getGraphics();
    g.setColor(UIManager.getColor("Panel.background"));
    g.fillRect(0, 0, bi.getWidth(), bi.getHeight());
    g.setColor(Color.BLACK);
    g.setFont(g.getFont().deriveFont(fontSize));
    FontMetrics fm = g.getFontMetrics();
    g.drawString(msg, (bi.getWidth() / 2) - fm.stringWidth(msg) / 2,
                 bi.getHeight() / 2 - fm.getHeight());
    if (msg2 != null) {
      g.drawString(msg2, (bi.getWidth() / 2) - fm.stringWidth(msg2) / 2,
                   bi.getHeight() / 2 + ((int) (fm.getHeight() * 1.5)));
    }
    return bi;
  }

  static float fontSize = 16f;

  private BufferedImage createIOExceptionImage(IOException e) {
    return createImage("Exception when loading image:", e.getMessage());
  }

  public static BufferedImage createReadyImage() {
    return createImage("Load Image Directory to Begin.", null);
  }

  private BufferedImage createNoFileImage() {
    return createImage("No image file found!", null);
  }

  private BufferedImage createMissingFileImage(String file) {
    return createImage("File missing:", file);
  }

  public boolean isRoot() {
    return isRoot;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public void setImageFile(String image) {
    this.imageFile = image;
  }

  public String getImageFile() {
    return this.imageFile;
  }

  public ArrayList<AnnotatedImage.Annotation> getAnnotations() {
    return this.annots;
  }

}
