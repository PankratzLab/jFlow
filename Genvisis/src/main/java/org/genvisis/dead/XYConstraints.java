package org.genvisis.dead;
// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 4/19/2007 9:12:19 AM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) 
// Source File Name:   XYConstraints.java


import java.io.Serializable;

public class XYConstraints
    implements Cloneable, Serializable
{
	static final long serialVersionUID = 1L;

    public XYConstraints()
    {
        this(0, 0, 0, 0);
    }

    public XYConstraints(int i, int j, int k, int l)
    {
        x = i;
        y = j;
        width = k;
        height = l;
    }

    public int getX()
    {
        return x;
    }

    public void setX(int i)
    {
        x = i;
    }

    public int getY()
    {
        return y;
    }

    public void setY(int i)
    {
        y = i;
    }

    public int getWidth()
    {
        return width;
    }

    public void setWidth(int i)
    {
        width = i;
    }

    public int getHeight()
    {
        return height;
    }

    public void setHeight(int i)
    {
        height = i;
    }

    public int hashCode()
    {
        return x ^ y * 37 ^ width * 43 ^ height * 47;
    }

    public boolean equals(Object obj)
    {
        if(obj instanceof XYConstraints)
        {
            XYConstraints xyconstraints = (XYConstraints)obj;
            return xyconstraints.x == x && xyconstraints.y == y && xyconstraints.width == width && xyconstraints.height == height;
        } else
        {
            return false;
        }
    }

    public Object clone()
    {
        return new XYConstraints(x, y, width, height);
    }

    public String toString()
    {
        return "XYConstraints[" + x + "," + y + "," + width + "," + height + "]";
    }

    int x;
    int y;
    int width;
    int height;
}