package org.genvisis.seq.analysis;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.genvisis.common.Logger;
import org.junit.Test;

/**
 * @author Travis Rogers
 */
public class BlastTableTest {

	/**
	 * Test of {@link BlastTable#createBlastTableFile()} method.
	 * 
	 * @throws IOException if a problem is encountered while cleaning up/deleting the output file
	 */
	@Test
	public void createBlastTableFileTEST() throws IOException {
		String[] searchTerms = {"chrX", "chrY"};
		String inFile = "src/test/resources/BlastTableTEST_input.vcf";
		String outFile = "src/test/resources/BlastTableTEST_actual_ouput.txt";
		String expectedOutput = "src/test/resources/BlastTableTEST_expected_output.txt";
		String token = "OFF_T_ALIGNMENTS";
		BlastTable instance = new BlastTable(token, searchTerms, inFile, outFile, new Logger());
		assertTrue(instance.createBlastTableFile());
		Path p1 = Paths.get(expectedOutput);
		Path p2 = Paths.get(outFile);
		try {
			List<String> expected = Files.readAllLines(p1, Charset.forName("utf-8"));
			List<String> actual = Files.readAllLines(p2, Charset.forName("utf-8"));
			assertEquals(expected.size(), actual.size());
			for (int i = 0; i < expected.size(); i++) {
				assertEquals(expected.get(i), actual.get(i));
			}
		} catch (IOException e) {
			fail(e.getMessage());
		} finally {
			Files.deleteIfExists(p2);
		}
	}

}
