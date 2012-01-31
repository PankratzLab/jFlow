package dead;

import java.io.*;
import java.util.*;

import mining.Distance;

import common.*;

public class GWorld {
	public static final String DIR = "C:/home/npankrat/auto/";
	public static final String[] TYPES = {"WILD", "SPECIAL", "TERRORIST CAMP", "WASTELAND", "CITY"};
	public static final String[] WILDS = {"GRASSLAND", "OIL FIELD", "HILLS", "MOUNTAIN", "RIVER", "PLAIN"};
	public static final String[] SPECIALS = {"TITANIUM", "GRAPHENE", "URANIUM", "DIAMONDS"};
	public static final String[] FREQUENT_ERRORS = {"REQUEST BOOST", "CityFieldsMap", "PRODUCTION BOOST"};

	public class Tile {
		String title;
		int x;
		int y;
		String region;
		int level;
		int type;
		Gen occupier;
		
		public Tile() {}
	}

	public class Gen {
		String name;
		int power;
		Alliance alliance;
		Vector<Tile> cities;
		Vector<Tile> wilds;
		Vector<String> powers;
		IntVector powerHistory;
		
		public Gen(String name) {
			this.name = name;
			this.power = -1;
			this.alliance = null;
			this.powers = new Vector<String>();
			this.cities = new Vector<Tile>();
			this.wilds = new Vector<Tile>();
			this.powerHistory = new IntVector();
		}
		
		public void determineMaxPower() {
			power = powers.size()==0?-1:Sort.putInOrder(Array.toIntArray(Array.toStringArray(powers)))[powers.size()-1];
			powerHistory.add(power);
		}
	}

	public class Alliance {
		String name;
		int power;
		Vector<String> members;
		Vector<String> powers;
		IntVector powerHistory;
		
		public Alliance(String name) {
			this.name = name;
			this.power = -1;
			this.members = new Vector<String>();
			this.powers = new Vector<String>();
			this.powerHistory = new IntVector();
		}

		public void determineMaxPower() {
			power = powers.size()==0?-1:Sort.putInOrder(Array.toIntArray(Array.toStringArray(powers)))[powers.size()-1];
			powerHistory.add(power);
		}
	}
	
	public class CoordIncer {
		private int x;
		private int y;
		private int incX;
		private int incY;
		
		public CoordIncer() {
			x = 5;
			y = 3;
			incX = 0;
			incY = 0;
		}
		
		public String getNext() {
			String coords;

			if (y > 800)
			{
				return null;
			}

			coords = (x-4+incX)+","+(y-2+incY);
			incX += 1;
			if (incX == 6) {
				incX = 0;
				incY += 1;

				if (x > 798) {
					incX = 4;
				}
			}
			if (incY == 4) {
				incY = 0;
				x += 6;

				if (x > 804)
				{
					x = 5;
					y += 4;
					incX = 0;
					incY = 0;
				}

				if (x > 798)
				{
					x = 799;
					incX = 4;
				}
			}
			
			return coords;
		}
	}

	private Hashtable<String, Tile> tiles;
	private Hashtable<String, Gen> gens;
	private Hashtable<String, Alliance> alliances;
	
	public GWorld() {
		tiles = new Hashtable<String, Tile>();
		gens = new Hashtable<String, Gen>();
		alliances = new Hashtable<String, Alliance>();
	}
	
	public GWorld(String filename) {
		this();
		
		BufferedReader reader;
//		PrintWriter writer, writer2, writer3;
		String[] lines;
		Tile tile;
		Alliance alliance;
		String[] keys;
		
		alliances.put("Specials", alliance = new Alliance("Specials"));
		for (int i = 0; i < SPECIALS.length; i++) {
			gens.put(SPECIALS[i], new Gen(SPECIALS[i]));
			alliance.members.add(SPECIALS[i]);
		}
		
		alliances.put("TERRORIST CONTROLLED", new Alliance("TERRORIST CONTROLLED"));

		try {
			reader = new BufferedReader(new FileReader(filename));
//			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_good_stuff.out"), true);
//			writer2 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"bad_stuff.out"), true);
//			writer3 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"superbad_stuff.out"), true);

			while (reader.ready()) {
				lines = getNextChunk(reader);
				
				if (lines.length > 6 && lines.length <= 12 && Array.max(ext.indexFactors(FREQUENT_ERRORS, lines, true, new Logger(), false, false)) == -1) {
//					writer.println(Array.toStr(lines, "\n")+"\n");
//					tile = parseTile(lines, writer3);
					tile = parseTile(lines);
					if (tile != null) {
						tiles.put(tile.x+","+tile.y, tile);
					} else {
						System.err.println(Array.toStr(lines, "\n"));
						System.err.println("----------------------");
					}
//				} else {
//					writer2.println(Array.toStr(lines, "\n")+"\n");
				}
			}
			reader.close();
//			writer.close();
//			writer2.close();
//			writer3.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		keys = HashVec.getKeys(gens);
		for (int i = 0; i < keys.length; i++) {
			gens.get(keys[i]).determineMaxPower();
		}

		keys = HashVec.getKeys(alliances);
		for (int i = 0; i < keys.length; i++) {
			alliances.get(keys[i]).determineMaxPower();
		}
	}
	
	public void mergeWorlds(GWorld trav) {
//		String[] tileKeys, genKeys, allianceKeys;
//		Alliance thisAlliance, travAlliance;
//		Gen thisGen, travGen;
//		Tile thisTile, travTile;
//		
//		tileKeys = HashVec.getKeys(trav.tiles);
//		genKeys = HashVec.getKeys(trav.gens);
//		allianceKeys = HashVec.getKeys(trav.alliances);
//
//		for (int i = 0; i < allianceKeys.length; i++) {
//			if (alliances.containsKey(allianceKeys[i])) {
//				travAlliance = trav.alliances.get(allianceKeys[i]);
//				thisAlliance = alliances.get(allianceKeys[i]);
//				for (int j = 0; j < travAlliance.members.size(); j++) {
//					HashVec.addIfAbsent(travAlliance.members.elementAt(j), thisAlliance.members);
//				}
//				thisAlliance.powerHistory.add(travAlliance.power);
//			} else {
//				alliances.put(allianceKeys[i], trav.alliances.get(allianceKeys[i]));
//			}
//			travAlliance = trav.alliances.get(allianceKeys[i]);
//		}
		
	}
	
	public static String[] getNextChunk(BufferedReader reader) {
		boolean done, started;
		Vector<String> v;
		String temp;
		
		v = new Vector<String>();

		try {
			done = false;
			started = false;
			while (!done) {
				if (reader.ready()) {
					temp = reader.readLine().trim();
					if (temp.equals("")) {
						if (started) {
							done = true; 
						}
					} else {
						v.add(temp);
						started = true;
					}
				} else {
					done = true;
				}
			}
		} catch (IOException ioe) {
			System.err.println("Error reading next chunk");
			System.exit(2);
		}
		
		return Array.toStringArray(v);
	}
	
	public Tile parseTile(String[] lines) {
		String genName, allianceName;
		Gen gen;
		Alliance alliance;
		String[] coords;
		Tile tile;
		boolean novel;

		tile = new Tile();
		tile.title = lines[0];
		if (ext.indexOfStr(tile.title, WILDS) >= 0) {
			tile.type = 0;
		} else if (ext.indexOfStr(tile.title, SPECIALS) >= 0) {
			tile.type = 1;
		} else if (tile.title.equals("TERRORIST CAMP")) {
			tile.type = 2;
		} else if (tile.title.equals("WASTELAND")) {
			tile.type = 3;
		} else {
			tile.type = 4;
		}
		
		if (tile.type <= 3) {
			if (!lines[1].startsWith("LEVEL ")) {
				System.err.println("Error - found a "+TYPES[tile.type]+" without a level");
				return null;
			} else {
				tile.level = Integer.parseInt(lines[1].substring(6));
				lines = Array.removeFromArray(lines, 1);
			}
		}
		// no scouting allowed if newbie or if under treaty
		if (lines[1].startsWith("SCOUT")) {
			lines = Array.removeFromArray(lines, 1);
		}
		
		//no attacking allowed if newbie or if under treaty
		if (lines[1].startsWith("ATTACK")) {
			lines = Array.removeFromArray(lines, 1);
		}

		// if belongs to your own alliance
		if (lines[1].startsWith("TRANSPORT")) {
			lines = Array.removeFromArray(lines, 1);
		}

		// if belongs to your own alliance
		if (lines[1].startsWith("REINFORCE")) {
			lines = Array.removeFromArray(lines, 1);
		}
		if (lines[1].startsWith("UNOCCUPIED")) {
			genName = null;
		} else {
			if (lines[1].startsWith("GENERAL")) {
				genName = lines[1].substring(8);
			} else if (tile.type != 2 && tile.type != 1) {
				System.err.println("Error - found a "+TYPES[tile.type]+" with an invalid occupier: "+lines[1]);
				return null;				
			} else {
				genName = lines[1]; // special occupied by an alliance
			}
		}
		if (tile.type == 4) {
			try {
				tile.level = Integer.parseInt(lines[2]);
				lines = Array.removeFromArray(lines, 2);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - found a "+TYPES[tile.type]+" that didn't have the level listed on line 5 (found '"+lines[2]+"' instead)");
				return null;
			}
		}
		
		if (genName == null) {
			if (lines[2].startsWith("POWER")) {
				lines = Array.removeFromArray(lines, 2);
			} else {
				System.err.println("Error - found a "+TYPES[tile.type]+" that should have had a power listed");
				return null;
			}
		} else {
			if (gens.containsKey(genName)) {
				gen = gens.get(genName);
			} else {
				gens.put(genName, gen = new Gen(genName));
			}

			if (tile.type == 1) {
//				make occupier the special itself and the title the alliance
				tile.occupier = gens.get(tile.title);
				tile.title = genName;

				if (alliances.containsKey(genName)) {
					alliance = alliances.get(genName);
				} else {
					alliances.put(genName, alliance = new Alliance(genName));
				}
				
				if (lines[2].startsWith("ALLIANCE POWER:")) {
					HashVec.addIfAbsent(ext.replaceAllWith(lines[2].substring(16).trim(), ",", ""), alliance.powers);
					lines = Array.removeFromArray(lines, 2);
				} else {
					System.err.println("Error - found a "+TYPES[tile.type]+" occupied by "+genName+" that didn't have an alliance power listed: "+lines[2]);
					return null;
				}
				if (lines[2].startsWith("ALLIANCE DIPLOMACY:")) {
					lines = Array.removeFromArray(lines, 2);
				} else {
					System.err.println("Error - found a "+TYPES[tile.type]+" occupied by "+genName+" that didn't have diplomacy listed");
					return null;
				}
				if (lines[2].startsWith("CONTROLLED:")) {
					lines = Array.removeFromArray(lines, 2);
				} else {
					System.err.println("Error - found a "+TYPES[tile.type]+" occupied by "+genName+" that didn't list how long it was controlled");
					return null;
				}
			} else {
				tile.occupier = gen;
				if (lines[2].startsWith("POWER")) {
					HashVec.addIfAbsent(ext.replaceAllWith(lines[2].substring(7), ",", ""), gen.powers);
					lines = Array.removeFromArray(lines, 2);
				} else {
					System.err.println("Error - found a "+TYPES[tile.type]+" occupied by "+genName+" that did not have a power listed");
//					writer3.println(Array.toStr(lines, "\n")+"\n");
					return null;
				}

				allianceName = "error!";
				if (lines[2].startsWith("ALLIANCE:")) {
					allianceName = lines[2].substring(10);
					if (alliances.containsKey(allianceName)) {
						alliance = alliances.get(allianceName);
					} else {
						alliances.put(allianceName, alliance = new Alliance(allianceName));
					}

					if (gen.alliance == null) {
						gen.alliance = alliance;
					} else if (!gen.alliance.name.equals(alliance.name)) {
						System.err.println("FYI - General "+genName+" has switched alliaces from "+gen.alliance.name+" to "+alliance.name);
						gen.alliance.members.remove(genName);
						gen.alliance = alliance;
					}
					HashVec.addIfAbsent(genName, gen.alliance.members);
					lines = Array.removeFromArray(lines, 2);
				} else if (!genName.equals("TERRORIST CONTROLLED")) {
					System.err.println("Error - found a "+TYPES[tile.type]+" occupied by "+genName+" that didn't have an alliance listed");
					return null;
				}
				if (lines[2].startsWith("ALLIANCE DIPLOMACY:")) {
					lines = Array.removeFromArray(lines, 2);
//				} else if (!allianceName.equals("None") && !genName.equals("TERRORIST CONTROLLED")) {
//					System.err.println("Error - found a "+TYPES[tile.type]+" occupied by "+genName+" that didn't have diplomacy listed");
//					writer3.println(Array.toStr(lines, "\n")+"\n");
//					return null;
				}
			}
		}
		if (lines[2].startsWith("COORDINATES:")) {
			coords = ext.replaceAllWith(lines[2].substring(13).trim(), " ", "").split("[,-]");
			tile.x = Integer.parseInt(coords[0]);
			tile.y = Integer.parseInt(coords[1]);
			tile.region = coords[2];
		} else {
			System.err.println("Error - found a "+TYPES[tile.type]+" that didn't have any coordinates!: '"+lines[2]+"'");
			return null;
		}

		if (tile.occupier != null) {
			novel = true;
			for (int i = 0; i < tile.occupier.cities.size(); i++) {
				if (tile.x == tile.occupier.cities.elementAt(i).x && tile.y == tile.occupier.cities.elementAt(i).y) {
					novel = false;
				}
			}
			for (int i = 0; i < tile.occupier.wilds.size(); i++) {
				if (tile.x == tile.occupier.wilds.elementAt(i).x && tile.y == tile.occupier.wilds.elementAt(i).y) {
					novel = false;
				}
			}
			if (novel) {
				if (tile.type == 4) {
					tile.occupier.cities.add(tile);
				} else {
					tile.occupier.wilds.add(tile);
				}
			}
		}
		
		return tile;
	}
	
	public static void clean(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String temp, trav;
		Hashtable<String, String> hash;
		Vector<String> v = new Vector<String>();
		String coord;
		
		hash = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_clean.out"));
			trav = "";
			coord = null;
			while (reader.ready()) {
				temp = reader.readLine().trim();
				if (temp.equals("")) {
					if (!trav.equals("")) {
						if (coord == null) {
							System.err.println(trav+"\n");
						} else if (hash.containsKey(coord)) {
							if (!trav.equals(hash.get(coord))) {
								v.add(coord);
							}
						} else {
							hash.put(coord, trav);
							writer.println(trav+"\n");
						}
					}
					trav = "";
					coord = null;
				} else {
					trav += temp+"\n";
					if (temp.startsWith("COORDINATES:")) {
						coord = temp;
					}
				}
				
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public String[] getStats() {
		String startCoord, endCoord;
		int count, start, end, good, total;
		String coords;
		boolean done;
		CoordIncer ci;
		
		count = 0;
		start = end = -1;
		good = total = 0;
		done = false;
		endCoord = null;
		startCoord = null;
		ci = new CoordIncer();
		while (!done) {
			count++;
			coords = ci.getNext();
			if (coords == null) {
				done = true;
			} else {
				if (tiles.containsKey(coords)) {
					if (start == -1) {
						start = count;
						startCoord = coords;
					}
					end = count;
					endCoord = coords;
					good++;
				}
			}
		}
		
		total = end - start + 1;
		
		return new String[] {startCoord, endCoord, good+"", total+"", ext.formDeci((double)good/(double)total, 3)};
	}
	
	public void generateReports(String root) {
		PrintWriter writer, writer2, writer3;
		String[] keys, peeps;
		int[] stats, order;
		Tile tile;
		Gen gen;	
		CoordIncer ci;
		String coords;
		boolean done;
				
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(root, false)+"_alliances.out"));
			writer2 = new PrintWriter(new FileWriter(ext.rootOf(root, false)+"_places.out"));
			writer3 = new PrintWriter(new FileWriter(ext.rootOf(root, false)+"_players.xln"));
			writer3.println("Player\tPower\tAlliance\tCities\tWilds");
			keys = HashVec.getKeys(alliances);
			for (int i = 0; i < keys.length; i++) {
				peeps = Array.toStringArray(alliances.get(keys[i]).members);
				if (peeps.length > 0) {
					writer.println(keys[i]);
					writer2.println(Array.toStr(Array.stringArray(keys[i].length(), "-"),""));
					writer2.println(keys[i]);
					writer2.println(Array.toStr(Array.stringArray(keys[i].length(), "-"),""));
					writer2.println();
					stats = new int[peeps.length];
					for (int j = 0; j < peeps.length; j++) {
						stats[j] = gens.get(peeps[j]).power;
					}
					order = Sort.quicksort(stats, Sort.DESCENDING);
					for (int j = 0; j < peeps.length; j++) {
						if (j < 20) {
							writer.println("  "+peeps[order[j]]+" ("+stats[order[j]]+")");
						}
						writer2.println(peeps[order[j]]+(stats[order[j]]==-1?"":" ("+stats[order[j]]+")"));
	
						gen = gens.get(peeps[order[j]]);
						for (int k = 0; k < gen.cities.size(); k++) {
							tile = gen.cities.elementAt(k);
							writer2.println(ext.formStr(tile.x+"", 3)+","+ext.formStr(tile.y+"", 3)+" "+ext.formStr(tile.region, 9)+"  "+tile.title);
						}
						writer2.println("-----");
						for (int k = 0; k < gen.wilds.size(); k++) {
							tile = gen.wilds.elementAt(k);
							writer2.println(ext.formStr(tile.x+"", 3)+","+ext.formStr(tile.y+"", 3)+" "+ext.formStr(tile.region, 9)+"  "+tile.title+" L"+tile.level);
						}
						writer3.println(peeps[order[j]]+"\t"+stats[order[j]]+"\t"+keys[i]+"\t"+gen.cities.size()+"\t"+gen.wilds.size());
						writer2.println();
					}
				
					writer.println();
				}
			}
			writer.close();
			writer2.close();
			writer3.close();
		} catch (Exception e) {
			System.err.println("Error writing to files");
			e.printStackTrace();
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(root, false)+"_new_terrs.out"));
			done = false;
			ci = new CoordIncer();
			while (!done) {
				coords = ci.getNext();
				if (coords == null) {
					done = true;
				} else {
					if (tiles.containsKey(coords)) {
						tile = tiles.get(coords);
						if (TYPES[tile.type].equals("TERRORIST CAMP")) {
							writer.println(tile.x+" "+tile.y+" "+tile.level);
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to files");
			e.printStackTrace();
		}
		
	}
	
	public void distributeTargets(int type, int level) {
		Vector<int[]> v;
		int count, limit;
		Tile tile;
		CoordIncer ci;
		String coords;
		boolean done;
		int[][] locations;
		int[][][] cityPositions;
		int[][][] orders;
		double[] dists;
		boolean[] used;
		String[][][] distributions;
		String[] gens;
		
		v = new Vector<int[]>();
		ci = new CoordIncer();
		done = false;
		while (!done) {
			coords = ci.getNext();
			if (coords == null) {
				done = true;
			} else {
				if (tiles.containsKey(coords)) {
					tile = tiles.get(coords);
					if (tile.type == type && tile.level == level) {
						v.add(new int[] {tile.x, tile.y});
					}
				}	
			}
		}
		
		locations = Matrix.toMatrix(v);
		gens = HashVec.loadFileToStringArray(DIR+"gens.txt", false, new int[] {0}, false);
		cityPositions = new int[gens.length][][];

		orders = new int[gens.length][][];
		count = 0;
		for (int i = 0; i < gens.length; i++) {
			cityPositions[i] = Matrix.toIntArrays(HashVec.loadFileToStringMatrix(DIR+gens[i]+"cities.txt", false, new int[] {0,1}, false));
			orders[i] = new int[cityPositions[i].length][];
			for (int j = 0; j < cityPositions.length; j++) {
				dists = new double[locations.length];
				for (int k = 0; k < locations.length; k++) {
					dists[k] = Distance.euclidean(cityPositions[i][j], locations[k]);
				}
				orders[i][j] = Sort.quicksort(dists);
				count++;
			}
		}

		limit = (int)Math.floor((double)locations.length /(double)count);
		distributions = new String[gens.length][][];
		used = new boolean[locations.length];		
		count = 0;
		while (count < limit) {
			for (int i = 0; i < gens.length; i++) {
				if (count == 0) {
					distributions[i] =  new String[cityPositions[i].length][limit];
				}
				for (int j = 0; j < cityPositions.length; j++) {
					for (int k = 0; k < locations.length; k++) {
						if (!used[orders[i][j][k]]) {
							distributions[i][j][count] = Array.toStr(locations[orders[i][j][k]], " ");
							System.out.println((i+1)+"\t"+(j+1)+"\t"+Array.toStr(locations[orders[i][j][k]]));
							used[orders[i][j][k]] = true;
							break;
						}
					}
				}
			}
			count++;
		}
		
		for (int i = 0; i < gens.length; i++) {
			for (int j = 0; j < cityPositions.length; j++) {
				Files.writeList(distributions[i][j], DIR+gens[i]+(j+1)+"TC"+level+"s.txt");
			}
		}
	}
	
	public void writeMissingCoords(String filename) {
		PrintWriter writer;
		CoordIncer ci;
		boolean done;
		String coords;
		
		try {
			writer = new PrintWriter(new FileWriter(filename));
			ci = new CoordIncer();
			done = false;
			while (!done) {
				coords = ci.getNext();
				if (coords == null) {
					done = true;
				} else {
					if (!tiles.containsKey(coords)) {
						writer.println(Array.toStr(coords.split(","), " "));
					}	
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}
	
	public static void stats(String dir) {
		PrintWriter writer;
		int count;
		boolean done;
		GWorld master, world;
		
		done = false;
		master = new GWorld();
		try {
			writer = new PrintWriter(new FileWriter(dir+"summary.out"));
			count = 1;
			while (!done) {
				if (Files.exists(dir+count+".txt", false)) {
					world = new GWorld(dir+count+".txt");
					System.out.println(count+"\t"+Array.toStr(world.getStats()));
					writer.println(count+"\t"+Array.toStr(world.getStats()));
					
					master = world;
					count++;
				} else {
					done = true;
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + dir+"summary.out");
			e.printStackTrace();
		}

		master.writeMissingCoords(dir+"missing_cords.out");
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Gwar.dat";

		String usage = "\n" + "dead.Gwar requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

//		String dir = "C:\\home\\npankrat\\auto\\fullMap\\";
		String dir = "C:\\home\\npankrat\\auto\\";
//		dir = "C:\\home\\npankrat\\auto\\fullMap\\dbs\\";
		
//		filename = dir+"highs!.txt";
		filename = dir+"db.txt";
//		filename = dir+"tcs.txt";
		filename = dir+"cities_db.txt";
		filename = dir+"terr410.txt";
		filename = dir+"drwar_db.txt";
		filename = dir+"db.txt";
		filename = dir+"rivera/all.txt";
		filename = dir+"allall.txt";
		filename = dir+"pain.txt";

		filename = dir+"GW/allall.txt";
		
		
//		filename = dir+"db_bloated.txt";

		filename = dir+"level3TCs.txt";
		
		try {
//			parse(filename);
//			clean(filename);
//			stats(dir);
//			new GWorld(filename).generateReports(filename);
			new GWorld(filename).distributeTargets(ext.indexOfStr("TERRORIST CAMP", TYPES), 3);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
