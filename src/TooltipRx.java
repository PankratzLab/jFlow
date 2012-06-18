//  <applet code="TooltipRx" width="400" height="400"></applet>
import java.awt.*;
import javax.swing.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import javax.swing.JApplet;
import javax.swing.JPanel;
import java.util.*;
 
public class TooltipRx extends JApplet
{
    @Override
    public void init()
    {
        Container con = getContentPane();
        con.setLayout(new BorderLayout());
        resize(new Dimension(300, 100));
        JScrollPane scroll = new JScrollPane();
        getContentPane().add(scroll, BorderLayout.CENTER);
//        setBackground(Color.WHITE);
        scroll.setViewportView(new ImagePanel());
    }
 
    private class ImagePanel extends JPanel implements MouseListener,
                                                       MouseMotionListener
    {
        ArrayList<Rectangle> rect_list;
        ArrayList <String> rect_descr;
        JPanel toolTip;
        JLabel label;
 
        public ImagePanel()
        {
            prepareData();
            addMouseListener(this);
            addMouseMotionListener(this);
 
            label = new JLabel();
            label.setHorizontalAlignment(JLabel.CENTER);
            label.setBorder(UIManager.getBorder("ToolTip.border"));
            toolTip = new JPanel(new BorderLayout());
            toolTip.setOpaque(true);
            toolTip.setBackground(UIManager.getColor("ToolTip.background"));
            toolTip.add(label);
            setLayout(null);
            add(toolTip);
        }
 
        public void prepareData()
        {
            int x=40;
            int y=50;
            int w=15;
            int ht=10;
            rect_list = new ArrayList<Rectangle>();
            rect_descr = new ArrayList <String> ();
 
            for(int i=0;i<5;i++)
            {
                rect_list.add(new Rectangle(x,y,w,ht));
                y +=20;
                x+=25;
            }
 
            for(int j=0;j<10;j++)
            {
                rect_descr.add("Base"+j);
            }
        }
 
        @Override
        protected void paintComponent(Graphics g)
        {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D)g;
            for(Rectangle r : rect_list)
            {
                g2.draw(r);
            }
        }
 
        public void setTooltip(int index, Point p)
        {
            label.setText(rect_descr.get(index));
            Dimension d = label.getPreferredSize();
            toolTip.setBounds(p.x, p.y, d.width, d.height);
            toolTip.setVisible(true);
            toolTip.repaint();
        }
 
        public boolean isToolTipShowing()
        {
            return toolTip.isVisible();
        }
 
        @Override
        public Dimension getPreferredSize()
        {
            return new Dimension(500,500);
        }
 
        public void hideToolTip()
        {
            toolTip.setVisible(false);
            toolTip.repaint();
        }
 
        public void mouseClicked(MouseEvent e) {}
        public void mousePressed(MouseEvent e) {}
        public void mouseReleased(MouseEvent e) {}
        public void mouseEntered(MouseEvent e) {}
        public void mouseExited(MouseEvent e) {}
        public void mouseDragged(MouseEvent e) {}
 
        public void mouseMoved(MouseEvent e)
        {
            Point p = e.getPoint();
            boolean traversing = false;
 
            Rectangle[] r =
                rect_list.toArray(new Rectangle[rect_list.size()]);
            for(int j = 0; j < r.length; j++)
            {
                Rectangle r1 = r[j];
 
                if(r1.contains(p))
                {
                    if(!isToolTipShowing())
                    {
                        // Offset toolTip from r1.
                        p.x = r1.x + r1.width + 10;
                        p.y = r1.y - 5;
                        setTooltip(j, p);
                    }
                    traversing = true;
                    break;
                }
            }
            if(!traversing && isToolTipShowing())
                hideToolTip();
        }
    }
}