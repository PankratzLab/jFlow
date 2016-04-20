package cnv.gui;

import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import net.miginfocom.swing.MigLayout;

import javax.swing.JButton;

import java.awt.Insets;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class JAccordionPanel extends JPanel {

    /**
     * Auto-generated svUID
     */
    private static final long serialVersionUID = 1L;

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
        
        expandoButton = new JButton(UP);
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
        expandoButton.setMargin(new Insets(0, 5, 0, 5));
        
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
