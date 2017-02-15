package org.genvisis.cnv;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.QueueControl;
import org.genvisis.common.QueueControl.JobQueue;

public class QueueProperties {
	
	private static final String PROPERTIES_FILE = "queues.properties";
	private static Properties props = null;
	private static final String INIT_FAILURE = "Failed to create PBS queue properties: " + PROPERTIES_FILE
      + ". Genvisis will continue, but generated PBS files should be reviewed before submission."	;
	
	private static Logger log;
	private static List<JobQueue> qList;
	
	private enum QueueKeys {
		Q_NAME() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getName(); }
		},
		Q_WALL_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinWalltime(); }
		},
		Q_WALL_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxWalltime(); }
		},
		Q_WALL_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultWalltime(); }
		},
		Q_MEM_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinMem(); }
		},
		Q_MEM_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxMem(); }
		},
		Q_MEM_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultMem(); }
		},
		Q_PROC_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinProc(); }
		},
		Q_PROC_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxProc(); }
		},
		Q_PROC_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultProc(); }
		},
		Q_NODE_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinNodeCnt(); }
		},
		Q_NODE_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxNodeCnt(); }
		},
		Q_NODE_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultNodeCnt(); }
		};
		
		public abstract String getValue(JobQueue q);
	}
	
	
	
	
	
	
	
	public static synchronized void voidQueues() {
		// clear queue container object
		
		if (Files.exists(PROPERTIES_FILE)) {
			boolean del = (new File(PROPERTIES_FILE)).delete();
			if (!del) {
				log.reportError("Could not delete PBS queue properties file.  Please delete this file manually and re-initialize.");
			}
		}
	}
	
	private static synchronized void init() {
		if (props == null) {
			props = new Properties();
			try {
				File propFile = new File(PROPERTIES_FILE);
				if (!propFile.exists() && !propFile.createNewFile()) {
					System.err.println(INIT_FAILURE);
					System.exit(-1);
				}
				InputStream is = new FileInputStream(propFile);
				props.load(is);
				is.close();
				
				List<JobQueue> qList = QueueControl.parseAllowedQueues(log);
				
				
				save();
			} catch (Exception e) {
				System.err.println(INIT_FAILURE);
				System.exit(-1);
			}
		}
	}
	

	/**
	 * Write the current properties to disk
	 */
	private static synchronized void save() {
		PrintWriter out;
		
		if (qList != null) {
			try {
				out = Files.getAppropriateWriter(PROPERTIES_FILE);
				QueueKeys[] keys = QueueKeys.values();
				
				for (JobQueue q : qList) {
					out.println(QueueKeys.Q_NAME.getValue(q));
					for (int i = 1; i < keys.length; i++) {
						out.println(keys[i].getValue(q));
					}
					out.println();
				}
				
				out.close();
			} catch (Exception e) {
				System.err.println("Failed to save PBS queue properties: " + PROPERTIES_FILE);
			}
		}
	}

	
}
