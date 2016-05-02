package cnv.plots;

import java.awt.geom.Ellipse2D;
//import java.awt.geom.Path2D;
import java.awt.geom.Path2D;

public class GenericPath {
    
    public Path2D myPath;
    private byte color;
    private byte fillColor;
    private byte layer;
    private boolean fill;
    private boolean editable;
    public double[][] foci;
    
    public GenericPath(Path2D path, byte color, byte fillColor, byte layer, boolean fill, boolean editable) {
        this.myPath = path;
        this.color = color;
        this.fillColor = fillColor;
        this.layer = layer;
        this.fill = fill;
        this.editable = editable;
    }

//    public Ellipse2D getPath() {
//        return myShape;
//    }
    
    public byte getColor() {
        return color;
    }

    public byte getFillColor() {
        return fillColor;
    }

    public byte getLayer() {
        return layer;
    }

    public boolean getFill() {
        return fill;
    }

    public boolean getEditable() {
        return editable;
    }
    
    
    
}
