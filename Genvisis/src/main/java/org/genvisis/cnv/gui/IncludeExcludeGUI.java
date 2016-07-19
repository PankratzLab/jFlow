package cnv.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Window;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JScrollPane;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;

import java.awt.Font;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.ArrayList;

public class IncludeExcludeGUI extends JDialog {

    private final JPanel contentPanel = new JPanel();
    private JPanel panel;
    private volatile int close = JOptionPane.CLOSED_OPTION;
    private ArrayList<JCheckBox> chks = new ArrayList<JCheckBox>();
    
    /**
     * Create the dialog.
     */
    public IncludeExcludeGUI(Window owner, String[] opts, boolean[] preSel) {
        super(owner);
        setModal(true);
        setBounds(100, 100, 450, 300);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new BorderLayout(0, 0));
        {
            JScrollPane scrollPane = new JScrollPane();
            scrollPane.setViewportBorder(null);
            scrollPane.setBorder(null);
            contentPanel.add(scrollPane);
            {
                panel = new JPanel();
                scrollPane.setViewportView(panel);
                StringBuilder rows = new StringBuilder();
                for (int i = 0; i < opts.length; i++) {
                    rows.append("[]");
                }
                panel.setLayout(new MigLayout("", "[grow,fill]", rows.toString()));
            }
        }
        {
            JLabel lblSelectItemsTo = new JLabel("<html><u>Select Items to Include:</u></html>");
            lblSelectItemsTo.setFont(new Font("Arial", Font.BOLD, 14));
            contentPanel.add(lblSelectItemsTo, BorderLayout.NORTH);
        }
        int rowInd = 0;
        for (String opt : opts) {
            JCheckBox chk = new JCheckBox(opt);
            chk.setSelected(preSel[rowInd]);
            chks.add(chk);
            panel.add(chk, "cell 0 " + rowInd);
            rowInd++;
        }
        
        {
            JPanel buttonPane = new JPanel();
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            {
                JButton okButton = new JButton("OK");
                okButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        close = JOptionPane.OK_OPTION;
                        setVisible(false);
                    }
                });
                okButton.setActionCommand("OK");
                buttonPane.add(okButton);
                getRootPane().setDefaultButton(okButton);
            }
            {
                JButton cancelButton = new JButton("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        close = JOptionPane.CANCEL_OPTION;
                        setVisible(false);
                        dispose();
                    }
                });
                cancelButton.setActionCommand("Cancel");
                buttonPane.add(cancelButton);
            }
        }
    }
    
    public int getCloseCode() {
        return close;
    }
    
    public boolean[] getSelected() {
        boolean[] sel = new boolean[chks.size()];
        for (int i = 0; i < sel.length; i++) {
            sel[i] = chks.get(i).isSelected();
        }
        return sel;
    }

}
