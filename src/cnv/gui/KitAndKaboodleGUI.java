package cnv.gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;

import javax.swing.JLabel;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.HashMap;

import javax.swing.JSeparator;
import javax.swing.JCheckBox;

import common.Array;

public class KitAndKaboodleGUI extends JDialog {

    private static final long serialVersionUID = 1L;

    private final JPanel contentPanel = new JPanel();
    
    String[][] OPTIONS = new String[][]{
            {"Create Marker Positions (if not already exists)",""},
            {"Parse Illumina Sample Files",""},
            {"Transpose Data into Marker-dominant Files",""},
            {"Extract Sample Data to lrrsd.xln File",""},
            {"Run Sex Checks",""},
            {"Create/Run PLINK Files",""}
    };
    volatile boolean[] selected = Array.booleanArray(OPTIONS.length, true);

    private HashMap<String, JCheckBox> checkBoxes = new HashMap<String, JCheckBox>();
    
    volatile boolean cancelled = false;

    /**
     * Create the dialog.
     */
    public KitAndKaboodleGUI() {
        setBounds(100, 100, 300, 400);
        getContentPane().setLayout(new BorderLayout());
        contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
        getContentPane().add(contentPanel, BorderLayout.CENTER);
        contentPanel.setLayout(new MigLayout("", "[grow]", "[][][]"));
        {
            JLabel lblKitAndKaboodle = new JLabel("Kit and Kaboodle Steps:");
            lblKitAndKaboodle.setFont(new Font("Arial", Font.BOLD, 16));
            contentPanel.add(lblKitAndKaboodle, "cell 0 0,alignx center");
        }
        {
            JSeparator separator = new JSeparator();
            contentPanel.add(separator, "cell 0 1,growx");
        }
        {
            for (int i = 0; i < OPTIONS.length; i++) {
                final int index = i;
                final JCheckBox chckbx = new JCheckBox();
                chckbx.setAction(new AbstractAction() {
                    static final long serialVersionUID = 1L;
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        selected[index] = chckbx.isSelected();
                    }
                });
                chckbx.setFont(chckbx.getFont().deriveFont(Font.PLAIN));
                chckbx.setSelected(selected[i]);
                chckbx.setToolTipText(OPTIONS[i][1]);
                chckbx.setText((i+1) + ": " + OPTIONS[i][0]);
                contentPanel.add(chckbx, "cell 0 "+ (i+2) + ",alignx left");
                checkBoxes.put(OPTIONS[i][0], chckbx);
            }
        }
        {
            JPanel buttonPane = new JPanel();
            buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
            getContentPane().add(buttonPane, BorderLayout.SOUTH);
            JButton okButton = new JButton("OK");
            okButton.setActionCommand("OK");
            buttonPane.add(okButton);
            getRootPane().setDefaultButton(okButton);
            JButton cancelButton = new JButton("Cancel");
            cancelButton.setActionCommand("Cancel");
            buttonPane.add(cancelButton);
            
            AbstractAction listener = new AbstractAction() {
                private static final long serialVersionUID = 1L;

                @Override
                public void actionPerformed(ActionEvent e) {
                    if (e.getActionCommand().equals("Cancel")) {
                        cancelled = true;
                    }
                    KitAndKaboodleGUI.this.setVisible(false);
                }
            };
            okButton.addActionListener(listener);
            cancelButton.addActionListener(listener);
            
        }
        pack();
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                cancelled = true;
                super.windowClosing(e);
            }
        });
    }
    
    public boolean[] getSelectedOptions() {
        return selected;
    }

    public boolean getCancelled() {
        return cancelled;
    }
    
}
