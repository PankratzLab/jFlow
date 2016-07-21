//package one.JL;
//
//import java.io.BufferedReader;
//import java.io.File;
//import java.io.FileNotFoundException;
//import java.io.IOException;
//
//import jxl.Workbook;
//import jxl.WorkbookSettings;
//import jxl.write.Label;
//import jxl.write.WritableSheet;
//import jxl.write.WritableWorkbook;
//import jxl.write.WriteException;
//import jxl.write.biff.RowsExceededException;
//import common.Files;
//import common.Logger;
//
//public class ExcelWriter {
//	private String[] files;
//	private String[] tabTitles;
//	private Logger log;
//
//	public ExcelWriter(String[] files, String[] tabTitles, Logger log) {
//		super();
//		this.files = files;
//		this.tabTitles = tabTitles;
//		this.log = log;
//	}
//
//	public void write(String output) {
//		try {
//			 
//			WritableWorkbook wb = Workbook.createWorkbook(new File(output));
//			WritableSheet[] sheets = new WritableSheet[tabTitles.length];
//			for (int i = 0; i < sheets.length; i++) {
//				sheets[i] = wb.createSheet(tabTitles[i], i);
//				
//			}
//			for (int i = 0; i < files.length; i++) {
//				try {
//					BufferedReader reader = Files.getAppropriateReader(files[i]);
//					int row = 0;
//					while (reader.ready()) {
//						String[] line = reader.readLine().trim().split("\t");
//						for (int j = 0; j < line.length; j++) {
//							Label label = new Label(j,row , line[j]);
//							try {
//								sheets[i].addCell(label);
//							} catch (RowsExceededException e) {
//								log.reportError("Too many Rows");
//								log.reportException(e);
//								// TODO Auto-generated catch block
//								e.printStackTrace();
//							} catch (WriteException e) {
//								// TODO Auto-generated catch block
//								e.printStackTrace();
//							}
//						}
//						row++;
//					}
//					reader.close();
//				} catch (FileNotFoundException fnfe) {
//					log.reportError("Error: file \"" + files[i] + "\" not found in current directory");
//					return;
//				} catch (IOException ioe) {
//					log.reportError("Error reading file \"" + files[i] + "\"");
//					return;
//				}
//
//			}
//			wb.write();
//			try {
//				wb.close();
//			} catch (WriteException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//
//		} catch (IOException e) {
//			log.reportTimeError("Could not create " + output);
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//
//	}
//
//}
