package org.genvisis.seq.qc;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.Plink;

/**
 * Utility class for validating information in a pedigree. For each family in a pedigree, a full
 * relationship map is derived. From this information, two items can be validated:
 * <ul>
 * <li>If a {@link SexChecks} output file is provided, the pedigree-declared sex will be validated
 * against the genome-estimated sex.</li>
 * <li>If a {@link Plink} genome file is provided, the expected familial relationship will be
 * validated against inherited-by-descent thresholds.</li>
 * </ul>
 */
public class CheckPed {
	private static final String GENOME_PATH = "./plink/quality_control/genome/plink.genome";
	private static final String SEXCHECK_PATH = "./results/sexCheck.xln";
	private static final Double RELATED_THRESHOLD = 0.1;

	/**
	 * {@link Project}-based {@link #validate(String, String, String, Double, Logger, String)} entry
	 * point, using the pedigree, genome and sexCheck files from the project, as needed.
	 */
	public static void validate(Project proj, String genomePath, String sexPath, Double minRel) {
		Pedigree pedigree = proj.loadPedigree();

		if (genomePath == null || !Files.exists(genomePath)) {
			genomePath = proj.GENOME_CLUSTER_FILENAME.getValue();
		}
		if (sexPath == null || !Files.exists(sexPath)) {
			sexPath = proj.SEXCHECK_RESULTS_FILENAME.getValue();
		}

		validate(pedigree.getPedigreeFile(), genomePath, sexPath, minRel, proj.getLog(),
		         proj.RESULTS_DIRECTORY.getValue());
	}

	/**
	 * Parses the given pedigree into familial mappings and validates these putative values against a
	 * genome and sexCheck analysis. The genome and sexCheck are optional but their absence will limit
	 * the usefulness of the output.
	 *
	 * @param genome Path to Plink .genome output (optional)
	 * @param sexCheck Path to SexChecks output (optional)
	 * @param minRel from [0 - 1.0]. Pedigree-computed relationships with an expected relatedness
	 *        below this threshold will be ignored.
	 * @param outputDir Path to write validation results.
	 */
	public static void validate(String pedigree, String genome, String sexCheck, Double minRel,
	                            Logger log, String outputDir) {
		Multimap<String, Person> families = ArrayListMultimap.create();

		// Map of fid\tiid to Person data structure
		Map<String, Person> peopleByFidIid = new HashMap<String, Person>();

		try {
			Hashtable<String, String> sexDict = null;

			// Load SexCheck file
			if (sexCheck != null && new File(sexCheck).exists()) {
				sexDict = HashVec.loadFileToHashString(sexCheck, new int[] {0},
				                                       new int[] {ext.indexOfStr("Sex",
				                                                                 SexChecks.SEX_HEADER)},
				                                       false, "\t", true, false, true);
			} else {
				if (sexCheck != null) {
					log.reportError("Sex check file: " + sexCheck + " is not valid.");
				}
				log.reportError("No sex data found. Only pedigree sexes will be used.");
			}

			// STEP 1 - parse pedigree to family networks.
			// Each fid maps to a map of iids to Person instances
			BufferedReader reader = Files.getAppropriateReader(pedigree);
			while (reader.ready()) {
				// Assume (FID, IID, FA, MO, SEX, PHENO, DNA, [MZID])
				String[] iPed = reader.readLine().split("\t");
				Person p = new Person(iPed[6], iPed[0], iPed[1], iPed[2], iPed[3], iPed[4],
				                      sexDict == null ? "0" : sexDict.get(iPed[6]),
				                      iPed.length > 7 ? iPed[7] : "");
				families.put(p.familyId(), p);
				peopleByFidIid.put(p.fidiid(), p);
			}
			reader.close();

			// Create links from parent > child, needed for determining expected relatedness
			for (Person child : peopleByFidIid.values()) {
				if (peopleByFidIid.containsKey(child.motherId())) {
					peopleByFidIid.get(child.motherId()).addChild(child);
				}
				if (peopleByFidIid.containsKey(child.fatherId())) {
					peopleByFidIid.get(child.fatherId()).addChild(child);
				}
			}

			// For each FID, we'll create an IIDxIID table, with 1 row per individual in the family and a
			// column entry for each related member of the family. The cell values will be used to store
			// expected/actual information about the given relationship.
			Map<String, Table<String, String, Relation>> familyTables =
			                                                          new HashMap<String, Table<String, String, Relation>>();

			// STEP 2 - Traverse family networks, populating expected relatedness values for all
			// relationships
			// in our pedigree
			for (String fid : families.keySet()) {
				Table<String, String, Relation> t = HashBasedTable.create();
				for (Person p : families.get(fid)) {
					Set<String> visited = new HashSet<String>();
					noteVisit(visited, p, p);
					// base, relative, relatedness, distance, searchParents, visited IID pairs, table
					populateRelations(p, p, 1.0, 0, true, visited, t, peopleByFidIid, minRel);
				}
				familyTables.put(fid, t);
			}

			// STEP 3 - Parse the genome and update Relations with computed values
			if (genome != null && new File(genome).exists()) {
				BufferedReader genomeReader = Files.getAppropriateReader(genome);
				int[] indices = ext.indexFactors(
				                                 new String[] {"FID1", "IID1", "FID2", "IID2", "Z0", "Z1",
				                                               "Z2", "PI_HAT"},
				                                 Plink.CLUSTER_HEADER, false, true);
				int fid1 = indices[0];
				int iid1 = indices[1];
				int fid2 = indices[2];
				int iid2 = indices[3];
				int z0 = indices[4];
				int z1 = indices[5];
				int z2 = indices[6];
				int piHat = indices[7];

				// Skip the header
				genomeReader.readLine();
				// process each line of the genome
				while (genomeReader.ready()) {
					String[] line = genomeReader.readLine().trim().split("\\s+");
					Relatedness r = getRel(line[z0], line[z1], line[z2], line[piHat]);
					if (line[fid1].equals(line[fid2])) {
						// Same family, record relatedness
						Table<String, String, Relation> t = familyTables.get(line[fid1]);
						// See if these people actually have a genetic relationship - e.g. spouses are typically
						// not related to each other
						if (t.contains(line[iid1], line[iid2])) {
							t.get(line[iid1], line[iid2]).setRelatedness(r);
							t.get(line[iid2], line[iid1]).setRelatedness(r);
						}
					} else if (!r.equals(Relatedness.UNRELATED)) {
						// Related from different families - add extraFamilial entry to both individuals
						peopleByFidIid.get(line[fid1] + "\t" + line[iid1]).addExtraFamilial(line[fid2],
						                                                                    line[iid2], r);
						peopleByFidIid.get(line[fid2] + "\t" + line[iid2]).addExtraFamilial(line[fid1],
						                                                                    line[iid1], r);
					}
				}
				genomeReader.close();

			} else {
				if (genome != null) {
					log.reportError("Invalid genome file: " + genome);
				}
				log.reportError("CheckPed: Without genome information, only pedigree-calculated relationships will be output.");
			}

			// write the output
			writeResults(peopleByFidIid, familyTables, outputDir, log);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Helper method to write {@link #validate} results.
	 */
	private static void writeResults(Map<String, Person> peopleByFidIid,
	                                 Map<String, Table<String, String, Relation>> familyTables,
	                                 String outputDir, Logger log) {
		if (outputDir == null) {
			outputDir = "";
		}

		if (!outputDir.isEmpty() && !new File(outputDir).exists()) {
			new File(outputDir).mkdirs();
		}

		String outFile = outputDir + "checkPed.csv";

		PrintWriter writer = Files.getAppropriateWriter(outFile);

		int maxRelationships = -1;

		// Find the max number of relationships possessed by a single person to use as column count
		for (String fid : familyTables.keySet()) {
			Table<String, String, Relation> t = familyTables.get(fid);
			for (Map<String, Relation> colMap : t.rowMap().values()) {
				// We subtract one since the columns should contain an identity cell
				int size = colMap.size() - 1;
				if (size > maxRelationships) {
					maxRelationships = size;
				}
			}
		}

		// Header format:
		// sample, fid, iid, putativesex, estimated sex, sexmismatch (1/0), pedmismatches (list of iid's
		// that failed validation),
		// [for each familial relationship - iid, pedrelationship, genomecall],
		// extrafamilials (semicolon-separated fid+iid(relatedness) entries)
		StringBuilder header =
		                     new StringBuilder("Sample\tfid\tiid\tpedSex\testimatedSex\tsexMismatch\tpedMismatches\t");

		// add an iid, pedCall and genomeCall column for each relationship
		for (int i = 0; i < maxRelationships; i++) {
			header.append("rel");
			header.append(i + 1);
			header.append("- iid\t");
			header.append("rel");
			header.append(i + 1);
			header.append("- pedCall\t");
			header.append("rel");
			header.append(i + 1);
			header.append("- genomeCall\t");
		}

		header.append("extrafamilials");
		writer.println(header.toString());

		// Write output by family
		for (String fid : familyTables.keySet()) {
			Map<String, Map<String, Relation>> family = familyTables.get(fid).rowMap();
			// Each row in the family map will be a row of output
			for (String iid : family.keySet()) {
				Person p = peopleByFidIid.get(fid + "\t" + iid);
				StringBuilder linePrefix = new StringBuilder();

				// Write basic info
				appendDelimited(linePrefix, p.dna());
				appendDelimited(linePrefix, fid);
				appendDelimited(linePrefix, iid);
				appendDelimited(linePrefix, p.pedSex());
				appendDelimited(linePrefix, p.estSex());
				appendDelimited(linePrefix, p.sexMisMatch());

				// For each non-identity column for this individual, write the iid, predicted relationship
				// and estimated relationship
				StringBuilder lineSuffix = new StringBuilder();
				StringBuilder pedMismatches = new StringBuilder("");
				Map<String, Relation> colMap = family.get(iid);
				//FIXME sort on relation first
				List<String> sortedFids = new ArrayList<String>(colMap.keySet());
				Collections.sort(sortedFids);
				for (String relId : sortedFids) {
					Relation pr = colMap.get(relId);
					if (!pr.isIdentity()) {
						appendDelimited(lineSuffix, relId);
						appendDelimited(lineSuffix, pr.getPedLabel());
						appendDelimited(lineSuffix, pr.getRelatedness());
						if (!pr.matches()) {
							if (pedMismatches.length() != 0) {
								pedMismatches.append(";");
							}
							pedMismatches.append(relId);
						}
					}
				}
				// fill empty columns
				for (int i = colMap.size() - 1; i < maxRelationships; i++) {
					for (int j = 0; j < 3; j++) {
						appendDelimited(lineSuffix, "-");
					}
				}
				appendDelimited(linePrefix, pedMismatches.toString());
				writer.print(linePrefix.toString());

				appendDelimited(lineSuffix, p.getExtraFamilials());
				writer.println(lineSuffix.toString());
			}
		}

		log.report("Wrote CheckPed output to: " + outFile);
		writer.flush();
		writer.close();
	}

	/**
	 * Helper method to append target string as a column for the current line
	 */
	private static void appendDelimited(StringBuilder sb, String s) {
		sb.append(s);
		sb.append("\t");
	}

	/**
	 * Updates the relatedness value in the given table for the base row and relation column and
	 * recurses to valid relationships.
	 */
	private static void populateRelations(Person base, Person relation, double relatedness,
	                                      int distance, boolean visitParents, Set<String> history,
	                                      Table<String, String, Relation> t,
	                                      Map<String, Person> peopleByFidIid, Double minRel) {
		// Account for MZ twins
		if (!base.individualId().equals(relation.individualId()) && base.mzEquals(relation)) {
			relatedness = 1.0;
		}

		// Updated relatedness for this relationship
		Relation r = t.get(base.individualId(), relation.individualId());
		if (r == null) {
			r = new Relation(base, relation, distance);
			t.put(base.individualId(), relation.individualId(), r);
		}
		r.addRelation(relatedness);

		// Everyone visited from this individual will gain half as much relatedness
		relatedness *= 0.5;

		// Cutoff point for minRel
		if (Double.compare(relatedness, minRel) < 0) {
			relatedness = 0;
		}

		// recursive step - parents
		if (visitParents) {
			// We can only visit parents if we have a direct parental lineage
			// e.g. if this is my nephew through my aunt, I'm not related to the father
			Person father = peopleByFidIid.get(relation.fatherId());
			populateRelationsRecursive(base, relation, father, relatedness, distance - 1, true, history,
			                           t, peopleByFidIid, minRel);
			Person mother = peopleByFidIid.get(relation.motherId());
			populateRelationsRecursive(base, relation, mother, relatedness, distance - 1, true, history,
			                           t, peopleByFidIid, minRel);
		}

		// recursive step - children
		for (Person child : relation.children()) {
			// If we're related to this person we're half as related to their children
			populateRelationsRecursive(base, relation, child, relatedness, distance + 1, false, history,
			                           t, peopleByFidIid, minRel);
		}
	}

	/**
	 * Recursive helper step for {@link #populateRelations}. Tests the relation > next linkage for
	 * validity. If valid, consumes the link and populates the base > next relationship.
	 */
	private static void populateRelationsRecursive(Person base, Person relation, Person next,
	                                               double relatedness, int distance,
	                                               boolean visitParents, Set<String> history,
	                                               Table<String, String, Relation> t,
	                                               Map<String, Person> peopleByFidIid,
	                                               Double minRel) {
		if (next != null && !visited(history, relation, next)) {
			noteVisit(history, relation, next);
			populateRelations(base, next, relatedness, distance, visitParents, history, t, peopleByFidIid,
			                  minRel);
		}
	}

	/**
	 * Create bidirectional entries in the visited set to ensure this node isn't visited this run
	 * through the same relationship again.
	 */
	private static void noteVisit(Set<String> visited, Person p1, Person p2) {
		visited.add(p1.individualId() + "\t" + p2.individualId());
		visited.add(p2.individualId() + "\t" + p1.individualId());
	}

	/**
	 * @return true if the relationship between these two was already used for node traversal in this
	 *         run.
	 */
	private static boolean visited(Set<String> visited, Person p1, Person p2) {
		return visited.contains(p1.individualId() + "\t" + p2.individualId());
	}

	/**
	 * Helper method for {@link #getRel(double, double, double, double)}. Takes care of parsing strings.
	 */
	public static Relatedness getRel(String ibd0, String ibd1, String ibd2, String piHat) {
		return getRel(ext.getValidDouble(ibd0), ext.getValidDouble(ibd1), ext.getValidDouble(ibd2),
		              ext.getValidDouble(piHat));
	}

	/**
	 * Helper method to get the best match {@link Relatedness} value for a set of IBD values from a
	 * plink .genome file.
	 */
	public static Relatedness getRel(double ibd0, double ibd1, double ibd2, double piHat) {
		//NB: relies on Relatedness values being declared in order of specificity
		for (Relatedness r : Relatedness.values()) {
			if (r.matches(ibd0, ibd1, ibd2, piHat)) {
				return r;
			}
		}
		return Relatedness.UNRELATED;
	}

	/**
	 * Helper class for tracking relatedness, generational distance, and estimaed sexes between two
	 * individuals.
	 */
	public static class Relation {
		private final String basePedSex;
		private final String relPedSex;
		private final String baseEstSex;
		private final String relEstSex;
		private final int distance;
		private final boolean identity;
		private final boolean siblings;
		private Relatedness computedRelatedness;
		private Relatedness pedEstRelatedness;
		private String label = null;
		private double relation;

		public Relation(Person base, Person relation, int distance) {
			basePedSex = base.pedSex();
			relPedSex = relation.pedSex();
			baseEstSex = base.estSex();
			relEstSex = relation.estSex();
			siblings = base.fatherId().equals(relation.fatherId())
			           || base.motherId().equals(relation.motherId());
			this.relation = 0.0;
			this.distance = distance;
			identity = base == relation;
		}

		public String getRelatedness() {
			return computedRelatedness == null ? "N/A" : computedRelatedness.label();
		}

		public void setRelatedness(Relatedness r) {
			this.computedRelatedness = r;
		}

		/**
		 * Each time this node is reached in {@link #populateRelations}, this method should be called.
		 * This number, from [0-1], represents the amount of genetic material we expect to share with
		 * this person given the particular relationship path that was followed.
		 */
		public void addRelation(double v) {
			relation += v;
		}

		/**
		 * @return True if this defines a relationship between a person and themselves.
		 */
		public boolean isIdentity() {
			return identity;
		}

		/**
		 * @return A common description of the relationship of the relative to the base individual
		 */
		public String getPedLabel() {
			if (label == null) {
				estimateRelatedness();
			}
			return label;
		}

		/**
		 * @return true if the genome-estimated Relatedness for this relationship matches what was
		 *         predicted by the pedigree.
		 */
		public boolean matches() {
			if (computedRelatedness == null) {
				return true;
			}

			if (pedEstRelatedness == null) {
				estimateRelatedness();
			}

			return pedEstRelatedness.equals(computedRelatedness);
		}

		/**
		 * Given the relation, distance and sex of the relation to this base Person, estimate a
		 * relationship label and Relatedness value.
		 */
		private void estimateRelatedness() {
			// TODO may need to track both distance in steps and relative generation ?
			// not sure I properly am accounting for half-X or X-removed
			if (relation >= 0.9) {
				pedEstRelatedness = Relatedness.DUPLICATE;
				if (isIdentity()) {
					label = "Duplicate";
				} else {
					label = "Twin";
				}
			} else if (relation >= 0.49) {
				if (distance < 0) {
					label = "1".equals(relPedSex) ? "Father" : "Mother";
					pedEstRelatedness = Relatedness.PARENTCHILD;
				} else if (distance > 0) {
					label = "1".equals(relPedSex) ? "Son" : "Daughter";
					pedEstRelatedness = Relatedness.PARENTCHILD;
				} else {
					label = "1".equals(relPedSex) ? "Brother" : "Sister";
					pedEstRelatedness = Relatedness.SIBLING;
				}
			} else if (relation >= 0.24) {
				pedEstRelatedness = Relatedness.SECOND;
				if (distance <= -2) {
					label = "1".equals(relPedSex) ? "Grandfather" : "Grandmother";
				} else if (distance == -1) {
					label = "1".equals(relPedSex) ? "Uncle" : "Aunt";
				} else if (distance == 0) {
					if (siblings) {
						label = "Half-sibling";
					} else {
						label = "Cousin";
					}
				} else if (distance == 1) {
					label = "1".equals(relPedSex) ? "Nephew" : "Niece";
				} else if (distance >= 2) {
					label = "1".equals(relPedSex) ? "Grandson" : "Granddaughter";
				}
			} else if (relation >= 0.11) {
				pedEstRelatedness = Relatedness.THIRD;
				if (distance <= -3) {
					label = "1".equals(relPedSex) ? "Great-grandfather" : "Great-grandmother";
				} else if (distance == -2) {
					label = "1".equals(relPedSex) ? "Great-uncle" : "Great-aunt";
				} else if (distance == -1) {
					label = "1".equals(relPedSex) ? "Half-uncle" : "Half-aunt"; // NB: could also be cousin
					                                                            // once-removed?
				} else if (distance == 0) {
					label = "Second cousin";
				} else if (distance == 1) {
					label = "1".equals(relPedSex) ? "Half-nephew" : "Half-niece";
				} else if (distance == 2) {
					label = "1".equals(relPedSex) ? "Great-nephew" : "Great-niece";
				} else if (distance >= 3) {
					label = "1".equals(relPedSex) ? "Great-grandson" : "Great-granddaughter";
				}
			} else {
				label = "Third cousin";
				pedEstRelatedness = Relatedness.FOURTH;
			}
		}
	}

	/**
	 * Enumerations of relationship levels of interest, unifying their label and plink genome
	 * thresholds.
	 */
	public enum Relatedness {
														DUPLICATE("Dup", 0, 0, 0,
														          0.9), PARENTCHILD("PO", 0, 0.8, 0,
														                            0.49), SIBLING("Sib", 0, 0.3, 0.1,
														                                           0.35), SECOND("2nd", 0, 0.4,
														                                                         0,
														                                                         0), THIRD("3rd",
														                                                                   0,
														                                                                   0.2,
														                                                                   0,
														                                                                   0), FOURTH("4th",
														                                                                              0,
														                                                                              0.1,
														                                                                              0,
														                                                                              0), UNRELATED("UN",
														                                                                                            0,
														                                                                                            0,
														                                                                                            0,
														                                                                                            0);

		private String label;
		private double ibd0;
		private double ibd1;
		private double ibd2;
		private double piHat;

		private Relatedness(String label, double ibd0, double ibd1, double ibd2, double piHat) {
			this.label = label;
			this.ibd0 = ibd0;
			this.ibd1 = ibd1;
			this.ibd2 = ibd2;
			this.piHat = piHat;
		}

		/**
		 * @return A shorthand label suitable for printing
		 */
		public String label() {
			return label;
		}

		/**
		 * @return true if the given values all meet or exceed the required thresholds for this
		 *         Relatedness.
		 */
		public boolean matches(double z0, double z1, double z2, double propIbd) {
			return Double.compare(z0, ibd0) >= 0 && Double.compare(z1, ibd1) >= 0
			       && Double.compare(z2, ibd2) >= 0 && Double.compare(propIbd, piHat) >= 0;
		}
	}

	/**
	 * Helper data structure for storing information about an individual
	 */
	public static class Person {

		private final String dna;
		private final String familyId;
		private final String individualId;
		private final String fatherId;
		private final String motherId;
		private final List<Person> children;
		private final String pedSex;
		private final String estSex;
		private final String mzId;
		private final String fidiid;
		private final Set<ExtraFamilialRel> extraFamilial;

		public Person(String dna, String fid, String iid, String fatherId, String motherId,
		              String pedSex, String estSex, String mzId) {
			this.dna = dna;
			this.familyId = fid;
			this.individualId = iid;
			this.fatherId = parentId(fid, fatherId);
			this.motherId = parentId(fid, motherId);
			this.pedSex = pedSex;
			this.estSex = estSex;
			this.mzId = mzId;
			fidiid = fid + "\t" + iid;
			children = new ArrayList<Person>();
			extraFamilial = new TreeSet<ExtraFamilialRel>();
		}

		private String parentId(String fid, String parent) {
			if (parent == null || parent.isEmpty()) {
				return "";
			}
			return fid + "\t" + parent;
		}

		public String dna() {
			return dna;
		}

		public String familyId() {
			return familyId;
		}

		public String individualId() {
			return individualId;
		}

		public String fatherId() {
			return fatherId;
		}

		public String motherId() {
			return motherId;
		}

		public List<Person> children() {
			return Collections.unmodifiableList(children);
		}

		public void addChild(Person c) {
			children.add(c);
		}

		public String pedSex() {
			return pedSex;
		}

		public String estSex() {
			return estSex;
		}

		public String fidiid() {
			return fidiid;
		}

		public String mzId() {
			return mzId;
		}

		/**
		 * @return 1 if expected/computed sex mismatches, 0 if they agree
		 */
		public String sexMisMatch() {
			if (pedSex().equals(estSex())) {
				return "0";
			}
			return "1";
		}

		/**
		 * @return true iff this and the targeted person both have the same, non-empty MZ-sibling ID
		 */
		public boolean mzEquals(Person o) {
			return !mzId.isEmpty() && mzId.equals(o.mzId());
		}

		public void addExtraFamilial(String fid, String iid, Relatedness r) {
			extraFamilial.add(new ExtraFamilialRel(fid, iid, r));
		}

		/**
		 * @return A semi-colon-delimited list of all {@link ExtraFamilialRel}s that were found to match
		 *         this individual.
		 */
		public String getExtraFamilials() {
			StringBuilder ef = new StringBuilder("");
			boolean first = true;
			for (ExtraFamilialRel r : extraFamilial) {
				if (!first) {
					ef.append(";");
				}
				ef.append(r.ref());
				first = false;
			}

			return ef.toString();
		}
	}

	/**
	 * Helper class for tracking relationships discovered in the genome that were not defined in the
	 * pedigree.
	 */
	private static class ExtraFamilialRel implements Comparable<ExtraFamilialRel> {
		final String fid;
		final String iid;
		final Relatedness r;

		public ExtraFamilialRel(String fid, String iid, Relatedness r) {
			this.fid = fid;
			this.iid = iid;
			this.r = r;
		}

		/**
		 * @return A string formatted: "fid+iid(label)". A '+' is used to avoid column overflow in
		 *         tab-delimited output. The label comes from {@link Relatedness#label()} for the
		 *         estimated relatedness between this and the base individual.
		 */
		public String ref() {
			return new StringBuilder(fid).append("+").append(iid).append("(").append(r.label())
			                             .append(")").toString();
		}

		@Override
		public int compareTo(ExtraFamilialRel o) {
			int c = r.compareTo(o.r);
			if (c == 0) {
				c = fid.compareTo(o.fid);
			}
			if (c == 0) {
				c = iid.compareTo(o.iid);
			}
			return c;
		}
	}

	/**
	 * CLI entry point
	 */
	public static void main(String... args) {
		final String pedigree = "ped";
		final String genome = "genome";
		final String sexCheck = "sexCheck";
		final String relatedThresh = "minRel";

		CLI c = new CLI(CheckPed.class);
		c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);
		c.addArg(pedigree, "standalone pedigree file to validate", "./pedigree.dat", CLI.Arg.FILE);
		c.addArg(sexCheck, "output from " + SexCheck.class.getName() + " for validating pedigree sexes",
		         SEXCHECK_PATH);
		c.addArg(genome, "plink .genome for validating pedigree relationships", GENOME_PATH);
		c.addArgWithDefault(relatedThresh,
		                    "minimum relatedness in .genome to identify pedigree mismatch",
		                    RELATED_THRESHOLD.toString(), CLI.Arg.NUMBER);

		c.addGroup(CLI.ARG_PROJ, pedigree);

		c.parseWithExit(args);

		String genomePath = c.get(genome);
		String sexPath = c.get(sexCheck);
		Double minRel = c.getD(relatedThresh);

		if (c.has(CLI.ARG_PROJ)) {
			validate(new Project(c.get(CLI.ARG_PROJ), false), genomePath, sexPath, minRel);
		} else {
			validate(c.get(pedigree), genomePath, sexPath, minRel,
			         new Logger("checkPed.log"), "");
		}
	}
}
