package common;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 * Mainly to consolidate files into a single excel workbook, likely more options to come...delimiters, numeric v string format, etc
 * 
 *
 */
public class ExcelConverter {

	private List<ExcelConversionParams> files;
	private String output;
	private Logger log;

	public ExcelConverter(ArrayList<String> files, String output, Logger log) {
		super();

		this.files = new ArrayList<ExcelConversionParams>();
		for (int i = 0; i < files.size(); i++) {
			this.files.add(new ExcelConversionParams(files.get(i), "\t", ext.rootOf(files.get(i))));
		}
		this.output = output;
		if (!output.endsWith(".xlsx")) {
			throw new IllegalArgumentException("output must have .xlsx extension");
		}
		this.log = log;
	}

	public ExcelConverter(List<ExcelConversionParams> files, String output, Logger log) {
		super();
		this.files = files;
		this.output = output;
		if (!output.endsWith(".xlsx")) {
			throw new IllegalArgumentException("output must have .xlsx extension");
		}
		this.log = log;
	}

	public void convert(boolean overwrite) {

		if (!Files.exists(output) || overwrite) {
			new File(ext.parseDirectoryOfFile(output)).mkdirs();
			try {
				XSSFWorkbook workbook = new XSSFWorkbook();
				for (ExcelConversionParams fileParams : files) {
					try {
						XSSFSheet sheet = workbook.createSheet(fileParams.getSheetName());
						try {
							BufferedReader reader = Files.getAppropriateReader(fileParams.getFile());
							int rownum = 0;
							while (reader.ready()) {
								String[] line = reader.readLine().trim().split(fileParams.getDelimiter());
								Row row = sheet.createRow(rownum++);
								int cellnum = 0;
								for (int i = 0; i < line.length; i++) {
									Cell cell = row.createCell(cellnum++);
									try {
										if (Double.isNaN(Double.parseDouble(line[i]))) {
											cell.setCellValue(line[i]);

										} else {
											cell.setCellValue(Double.parseDouble(line[i]));
										}
									} catch (NumberFormatException nfe) {
										cell.setCellValue(line[i]);

									}

								}
							}

							reader.close();
						} catch (FileNotFoundException fnfe) {
							log.reportError("Error: file \"" + fileParams.getFile() + "\" not found in current directory");
							workbook.close();
							return;
						} catch (IOException ioe) {
							log.reportError("Error reading file \"" + fileParams.getFile() + "\"");
							workbook.close();
							return;
						}
					} catch (IllegalArgumentException ile) {
						log.reportTimeError("Offending = " + fileParams.getSheetName());
						ile.printStackTrace();
					}
				}

				FileOutputStream out = new FileOutputStream(new File(output));
				workbook.write(out);
				out.close();
				workbook.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		} else {
			log.reportFileExists(output);
		}

	}

	public static class ExcelConversionParams {
		private String file;

		private String delimiter;
		private String sheetName;

		/**
		 * @param file
		 *            the file to convert
		 * @param delimiter
		 *            to parse to cells
		 * @param sheetName
		 *            name of this sheet in the workbook
		 */
		public ExcelConversionParams(String file, String delimiter, String sheetName) {
			super();
			this.file = file;
			this.delimiter = delimiter;
			this.sheetName = sheetName;
		}

		public String getFile() {
			return file;
		}

		public String getDelimiter() {
			return delimiter;
		}

		public String getSheetName() {
			return sheetName;
		}

	}

}
