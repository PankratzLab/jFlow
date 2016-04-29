package cnv.plots;

import java.awt.geom.Path2D;

public class GenericPath {
    
    private Path2D myShape;
    private byte color;
    private byte fillColor;
    private byte layer;
    private boolean fill;
    private boolean editable;
    
    
    public GenericPath(Path2D myShape, byte color, byte fillColor, byte layer, boolean fill, boolean editable) {
        this.myShape = myShape;
        this.color = color;
        this.fillColor = fillColor;
        this.layer = layer;
        this.fill = fill;
    }

    public Path2D getPath() {
        return myShape;
    }
    
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
