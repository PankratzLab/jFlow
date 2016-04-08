package one.ben.fcs;

import java.awt.Font;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.Format;
import java.text.NumberFormat;

import javax.swing.JPanel;

import net.miginfocom.swing.MigLayout;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JLabel;
import javax.swing.JComboBox;
import javax.swing.SwingUtilities;
import javax.swing.event.ListDataListener;

import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;

import javax.swing.JTextField;
import javax.swing.JFormattedTextField;
import javax.swing.JCheckBox;
import javax.swing.SwingConstants;

public class FCSPlotControlPanel extends JPanel {

    private FCSPlot plot;
    
    private JComboBox<PLOT_TYPE> cbType;
    private JComboBox<String> cbYData;
    private JComboBox<String> cbXData;
    private JComboBox<AXIS_SCALE> cbYScale;
    private JComboBox<AXIS_SCALE> cbXScale;
    
    private Font lblFont = new Font("Arial", 0, 12);
    private JFormattedTextField yBndsMin;
    private JFormattedTextField yBndsMax;
    private JFormattedTextField xBndsMin;
    private JFormattedTextField xBndsMax;
    
    /**
     * Create the panel.
     */
    public FCSPlotControlPanel(FCSPlot plot) {
        this.plot = plot;
        
        setLayout(new MigLayout("", "[][][grow][]", "[][][][][][][][][][][][][][]"));
        
        JLabel lblPlotType = new JLabel("Plot Type:");
        lblPlotType.setFont(lblFont);
        add(lblPlotType, "cell 0 3,alignx trailing");
        
        cbType = new JComboBox();//<PLOT_TYPE>(PLOT_TYPE.values());
        cbType.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setPlotType((PLOT_TYPE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        add(cbType, "cell 1 3 3 1,growx");
        
        JLabel lblYaxisData = new JLabel("Y-Axis Data:");
        lblYaxisData.setFont(lblFont);
        add(lblYaxisData, "cell 0 5,alignx trailing");
        
        cbYData = new JComboBox<String>();
        cbYData.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setYDataName(arg0.getItem().toString());
                    plot.updateGUI();
                }
            }
        });
        cbYData.setMaximumRowCount(15);
        add(cbYData, "cell 1 5 3 1,growx");
        
        JLabel lblScale = new JLabel("Scale:");
        lblScale.setFont(lblFont);
        add(lblScale, "cell 0 6 2 1,alignx trailing");
        
        chckbxShowMedianY = new JCheckBox("Show Median");
        chckbxShowMedianY.setHorizontalAlignment(SwingConstants.TRAILING);
        chckbxShowMedianY.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setMedianVisible(show, true);
                plot.updateGUI();
            }
        });
        add(chckbxShowMedianY, "cell 0 7 3 1,alignx trailing");
        
        chckbxShowSdY = new JCheckBox("Show SD");
        chckbxShowSdY.setHorizontalAlignment(SwingConstants.TRAILING);
        chckbxShowSdY.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setSDVisible(show, true);
                plot.updateGUI();
            }
        });
        add(chckbxShowSdY, "cell 3 7,alignx trailing");
        
        JLabel lblYbounds = new JLabel("Y-Bounds:");
        lblYbounds.setFont(lblFont);
        add(lblYbounds, "cell 0 8 2 1,alignx trailing");
        
        JLabel lblXaxisData = new JLabel("X-Axis Data:");
        lblXaxisData.setFont(lblFont);
        add(lblXaxisData, "cell 0 10,alignx trailing");
        
        cbXData = new JComboBox<String>();
        cbXData.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setXDataName(arg0.getItem().toString());
                    plot.updateGUI();
                }
            }
        });
        cbXData.setMaximumRowCount(15);
        add(cbXData, "cell 1 10 3 1,growx");
        
        cbYScale = new JComboBox();//<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
        cbYScale.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setYScale((AXIS_SCALE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        add(cbYScale, "cell 2 6 2 1,growx");
        
        cbXScale = new JComboBox();//<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
        cbXScale.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setXScale((AXIS_SCALE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        
        JLabel lblScale_1 = new JLabel("Scale:");
        lblScale_1.setFont(lblFont);
        add(lblScale_1, "cell 0 11 2 1,alignx trailing");
        add(cbXScale, "cell 2 11 2 1,growx");
        
        Format numberFormat = NumberFormat.getNumberInstance();
        yBndsMin = new JFormattedTextField(numberFormat);
        yBndsMin.addPropertyChangeListener("value", pcl);
        yBndsMin.setColumns(10);
        yBndsMin.setValue(0);
        add(yBndsMin, "flowx,cell 2 8 2 1");
        
        yBndsMax = new JFormattedTextField(numberFormat);
        yBndsMax.addPropertyChangeListener("value", pcl);
        yBndsMax.setColumns(10);
        add(yBndsMax, "cell 2 8 2 1");
        
        chckbxShowMedianX = new JCheckBox("Show Median");
        chckbxShowMedianX.setHorizontalAlignment(SwingConstants.TRAILING);
        chckbxShowMedianX.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setMedianVisible(show, false);
                plot.updateGUI();
            }
        });
        add(chckbxShowMedianX, "cell 0 12 3 1,alignx right");
        
        chckbxShowSdX = new JCheckBox("Show SD");
        chckbxShowSdX.setHorizontalAlignment(SwingConstants.TRAILING);
        chckbxShowSdX.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setSDVisible(show, false);
                plot.updateGUI();
            }
        });
        add(chckbxShowSdX, "cell 3 12,alignx right");
        
        JLabel lblXbounds = new JLabel("X-Bounds:");
        lblXbounds.setFont(lblFont);
        add(lblXbounds, "cell 0 13 2 1,alignx trailing");
        
        xBndsMin = new JFormattedTextField(numberFormat);
        xBndsMin.addPropertyChangeListener("value", pcl);
        xBndsMin.setColumns(10);
        xBndsMin.setValue(0);
        add(xBndsMin, "flowx,cell 2 13 2 1");
        
        xBndsMax = new JFormattedTextField(numberFormat);
        xBndsMax.addPropertyChangeListener("value", pcl);
        xBndsMax.setColumns(10);
        add(xBndsMax, "cell 2 13 2 1");
        
    }
    
    PropertyChangeListener pcl = new PropertyChangeListener() {
        @Override
        public void propertyChange(PropertyChangeEvent evt) {
            if (evt == null || plot == null || progSet) return;
            JFormattedTextField jftf = (JFormattedTextField) evt.getSource();
            String prop = evt.getPropertyName();
            if (jftf == xBndsMin) {
                prop = AbstractPanel2.X_MIN;
            } else if (jftf == xBndsMax) {
                prop = AbstractPanel2.X_MAX;
            } else if (jftf == yBndsMin) {
                prop = AbstractPanel2.Y_MIN;
            } else if (jftf == yBndsMax) {
                prop = AbstractPanel2.Y_MAX;
            }
            Object oldV = evt.getOldValue();
            Object newV = evt.getNewValue();
            Double oldV2 = oldV == null ? null : ((Number) oldV).doubleValue(); 
            Double newV2 = newV == null ? null : ((Number) newV).doubleValue();
            /*FCSPlotControlPanel.this.*/firePropertyChange(prop, oldV2, newV2);
//            plot.firePropertyChange(prop, oldV2, newV2);
//            PropertyChangeEvent pce = new PropertyChangeEvent(FCSPlotControlPanel.this, prop, oldV2, newV2);
//            plot.propertyChange(pce);
            
        }
    };
    
    public void setPlotType(PLOT_TYPE typ) {
        cbType.setSelectedItem(typ);
    }
    
    public void setScale(AXIS_SCALE scl, boolean x) {
        (x ? cbXScale : cbYScale).setSelectedItem(scl);
    }
    
    public void setColumns(String[] dataNames, boolean x, int selected) {
        (x ? cbXData : cbYData).setModel(new DefaultComboBoxModel<String>(dataNames));
        (x ? cbXData : cbYData).setSelectedIndex(selected);
        (x ? cbXData : cbYData).repaint();
    }
    
    volatile boolean progSet = false;
    private JCheckBox chckbxShowMedianX;
    private JCheckBox chckbxShowSdX;
    private JCheckBox chckbxShowMedianY;
    private JCheckBox chckbxShowSdY;

    private void resetProgSet() {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                progSet = false;
            }
        });
    }
    
    public void setXMin(double xMin) {
        progSet = true;
        xBndsMin.setValue(xMin);
        resetProgSet();
    }
    public void setXMax(double xMax) {
        progSet = true;
        xBndsMax.setValue(xMax);
        resetProgSet();
    }
    public void setYMin(double yMin) {
        progSet = true;
        yBndsMin.setValue(yMin);
        resetProgSet();
    }
    public void setYMax(double yMax) {
        progSet = true;
        yBndsMax.setValue(yMax);
        resetProgSet();
    }
    
}
