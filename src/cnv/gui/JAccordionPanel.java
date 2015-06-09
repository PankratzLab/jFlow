package cnv.gui;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import net.miginfocom.swing.MigLayout;

import javax.swing.JButton;

import java.awt.Insets;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class JAccordionPanel extends JPanel {

    public static void main(String[] args) {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        frame.getContentPane().add(new JAccordionPanel());
        frame.setVisible(true);
        frame.pack();
    }
    
    
    private static final String UP = "/\\";
    private static final String DOWN = "\\/";

    private JButton expandoButton;
    public JPanel topPanel;
    public JPanel contentPanel;
    private JPanel panel;
    
    /**
     * Create the panel.
     */
    public JAccordionPanel() {
        setLayout(new MigLayout("hidemode 3", "[grow]", "[grow]"));
        
        topPanel = new JPanel();
        topPanel.setBorder(null);
        add(topPanel, "north");
        topPanel.setLayout(new MigLayout("", "[grow]", "[]"));
        
        panel = new JPanel();
        topPanel.add(panel, "east");
        panel.setLayout(new MigLayout("", "[]", "[]"));
        
        expandoButton = new JButton(DOWN);
        panel.add(expandoButton, "cell 0 0");
        expandoButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (contentPanel.isVisible()) {
                    shrink();
                } else {
                    expand();
                }
            }
        });
        expandoButton.setMargin(new Insets(0, 2, 0, 2));
        
        contentPanel = new JPanel();
        add(contentPanel, "cell 0 0,grow");
        contentPanel.setLayout(new MigLayout("", "[]", "[]"));
    }
    
    public void expand() {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                contentPanel.setVisible(true);
                expandoButton.setText(UP);
                JAccordionPanel.this.invalidate();
                JAccordionPanel.this.repaint();
            }
        });
    }
    
    public void shrink() {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                contentPanel.setVisible(false);
                expandoButton.setText(DOWN);
                JAccordionPanel.this.invalidate();
                JAccordionPanel.this.repaint();
            }
        });
    }
    
    
}
